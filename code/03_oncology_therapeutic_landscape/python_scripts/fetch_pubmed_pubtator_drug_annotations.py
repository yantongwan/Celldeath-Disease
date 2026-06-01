#!/usr/bin/env python3
"""Fetch and merge PubMed/PubTator chemical evidence for oncology drug analysis.

This script upgrades the exploratory local drug lexicon workflow by adding:
1. PubMed EFetch XML ChemicalList annotations.
2. PubTator3 chemical entity annotations.
3. Local title/abstract lexicon fallback from the existing oncology module.
4. Sentence-level drug/death-mode relation candidates with blank manual labels.

It never invents missing counts or annotations. Network calls are opt-in via
--run or --run-all, and every response is cached.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MODULE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = MODULE_ROOT.parents[1]

DEFAULT_ARTICLE_INPUT = MODULE_ROOT / "12_source_data" / "oncology_article_records_scoped.csv"
DEFAULT_LOCAL_LEXICON = MODULE_ROOT / "05_drug_therapy_extraction" / "drug_lexicon_normalized.csv"
DEFAULT_LOCAL_MENTIONS = MODULE_ROOT / "05_drug_therapy_extraction" / "oncology_drug_mentions.csv"

OUT_DRUG = MODULE_ROOT / "05_drug_therapy_extraction"
OUT_REL = MODULE_ROOT / "09_sentence_level_relation_audit"
OUT_QC = MODULE_ROOT / "15_QC_reports"
DEFAULT_CACHE = OUT_DRUG / "pubmed_pubtator_cache"

EUTILS_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PUBTATOR_BIOCJSON = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson"


DEATH_VARIANTS = {
    "ferroptosis": r"ferroptosis|ferroptotic",
    "pyroptosis": r"pyroptosis|pyroptotic",
    "netosis": r"netosis|neutrophil extracellular trap|neutrophil extracellular traps",
    "necroptosis": r"necroptosis|necroptotic",
    "immunogenic cell death": r"immunogenic cell death|\bICD\b",
    "icd": r"immunogenic cell death|\bICD\b",
    "cuproptosis": r"cuproptosis|cuproptotic",
    "panoptosis": r"panoptosis|panoptotic",
    "disulfidptosis": r"disulfidptosis|disulfidptotic",
}

RELATION_CUES = [
    ("induction_or_activation", r"\b(induce|induced|induces|induction|trigger|triggered|promote|promoted|activate|activated|enhance|enhanced|increase|increased)\b"),
    ("sensitization", r"\b(sensitize|sensitized|sensitise|sensitised|enhance sensitivity|overcome resistance)\b"),
    ("inhibition_or_protection", r"\b(inhibit|inhibited|suppress|suppressed|block|blocked|prevent|prevented|protect|protected)\b"),
    ("resistance_escape", r"\b(resistance|resistant|escape)\b"),
    ("combination_therapy", r"\b(combination|combined|synergize|synergistic|synergy|co-treatment|cotreatment)\b"),
]

NEGATION_RE = re.compile(r"\b(not|failed to|did not|no significant|without)\b", re.I)
TUMOR_CONTEXT_RE = re.compile(r"cancer|tumou?r|neoplasm|carcinoma|glioma|melanoma|leukemia|lymphoma|sarcoma|malignan", re.I)


@dataclass
class HttpResult:
    ok: bool
    text: str
    source: str
    cache_path: Path
    error: str = ""
    status: str = ""


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_text(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def progress(msg: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def normalize_space(text: str | None) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def normalize_name(text: str | None) -> str:
    text = normalize_space(text).lower()
    text = text.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return normalize_space(text)


def boundary_regex(term: str) -> re.Pattern:
    escaped = re.escape(term)
    return re.compile(rf"(^|[^A-Za-z0-9])({escaped})([^A-Za-z0-9]|$)", re.I)


def synonym_regex(synonyms: str) -> re.Pattern:
    parts = [p.strip() for p in synonyms.split("|") if p.strip()]
    escaped = "|".join(re.escape(p) for p in parts)
    return re.compile(rf"(^|[^A-Za-z0-9])({escaped})([^A-Za-z0-9]|$)", re.I)


def batched(items: list[str], size: int) -> Iterable[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def batch_hash(prefix: str, pmids: list[str]) -> str:
    key = prefix + ":" + ",".join(pmids)
    return hashlib.sha256(key.encode("utf-8")).hexdigest()[:20]


def request_with_cache(
    url: str,
    params: dict[str, str],
    cache_path: Path,
    sleep_seconds: float,
    force: bool,
    run: bool,
    timeout: int = 60,
) -> HttpResult:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists() and not force:
        return HttpResult(True, cache_path.read_text(encoding="utf-8"), "cache", cache_path, status="cached")
    if not run:
        return HttpResult(False, "", "dry_run", cache_path, error="dry-run; network not requested", status="dry_run")
    query = urllib.parse.urlencode(params)
    request_url = f"{url}?{query}"
    try:
        time.sleep(max(sleep_seconds, 0))
        request = urllib.request.Request(request_url, headers={"User-Agent": "cell-death-atlas-drug-annotation/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
        cache_path.write_text(text, encoding="utf-8")
        return HttpResult(True, text, "network", cache_path, status="fetched")
    except urllib.error.HTTPError as exc:
        return HttpResult(False, "", "network", cache_path, error=f"HTTP {exc.code}: {exc.reason}", status="error")
    except Exception as exc:  # noqa: BLE001
        return HttpResult(False, "", "network", cache_path, error=str(exc), status="error")


def fetch_pubmed_efetch(
    pmid_batch: list[str],
    args: argparse.Namespace,
) -> HttpResult:
    params = {
        "db": "pubmed",
        "retmode": "xml",
        "id": ",".join(pmid_batch),
        "tool": args.tool,
    }
    if args.email:
        params["email"] = args.email
    if args.api_key:
        params["api_key"] = args.api_key
    cache = args.cache_dir / "pubmed_efetch_xml" / f"{batch_hash('efetch', pmid_batch)}.xml"
    return request_with_cache(EUTILS_EFETCH, params, cache, args.sleep_seconds, args.force, args.run or args.run_all)


def fetch_pubtator(
    pmid_batch: list[str],
    args: argparse.Namespace,
) -> HttpResult:
    params = {
        "pmids": ",".join(pmid_batch),
        "concepts": "chemical",
    }
    cache = args.cache_dir / "pubtator_biocjson" / f"{batch_hash('pubtator', pmid_batch)}.json"
    return request_with_cache(PUBTATOR_BIOCJSON, params, cache, args.sleep_seconds, args.force, args.run or args.run_all)


def parse_pubmed_chemical_xml(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not text.strip():
        return rows
    root = ET.fromstring(text)
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//MedlineCitation/PMID")
        pmid = pmid_el.text.strip() if pmid_el is not None and pmid_el.text else ""
        if not pmid:
            continue
        for chemical in article.findall(".//ChemicalList/Chemical"):
            name_el = chemical.find("NameOfSubstance")
            reg_el = chemical.find("RegistryNumber")
            chemical_name = normalize_space(name_el.text if name_el is not None else "")
            if not chemical_name:
                continue
            rows.append(
                {
                    "pmid": pmid,
                    "chemical_name": chemical_name,
                    "chemical_name_normalized": normalize_name(chemical_name),
                    "mesh_ui": name_el.attrib.get("UI", "") if name_el is not None else "",
                    "registry_number": normalize_space(reg_el.text if reg_el is not None else ""),
                    "source": "PubMed_MeSH_ChemicalList",
                }
            )
    return rows


def parse_pubtator_biocjson(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not text.strip():
        return rows
    data = json.loads(text)
    if isinstance(data, list):
        documents = data
    elif isinstance(data, dict):
        # PubTator3 BioC JSON responses currently use a top-level "PubTator3"
        # array. Older BioC-style responses may use "documents".
        documents = data.get("documents") or data.get("PubTator3") or []
    else:
        documents = []
    for doc in documents:
        pmid = str(doc.get("id", "")).split("|")[0]
        for passage in doc.get("passages", []):
            section = passage.get("infons", {}).get("section_type", passage.get("infons", {}).get("type", ""))
            for ann in passage.get("annotations", []):
                infons = ann.get("infons", {})
                ann_type = (infons.get("type") or infons.get("biotype") or "").lower()
                if "chemical" not in ann_type:
                    continue
                text_mention = normalize_space(ann.get("text", ""))
                if not text_mention:
                    continue
                loc = ann.get("locations", [{}])[0] if ann.get("locations") else {}
                rows.append(
                    {
                        "pmid": pmid,
                        "chemical_name": text_mention,
                        "chemical_name_normalized": normalize_name(text_mention),
                        "pubtator_identifier": infons.get("identifier", ""),
                        "pubtator_type": ann_type,
                        "mention_text": text_mention,
                        "mention_section": section,
                        "start_offset": loc.get("offset", ""),
                        "end_offset": int(loc.get("offset", 0)) + int(loc.get("length", 0)) if str(loc.get("offset", "")).isdigit() and str(loc.get("length", "")).isdigit() else "",
                        "source": "PubTator3_chemical",
                    }
                )
    return rows


def load_local_lexicon(path: Path) -> tuple[dict[str, dict[str, str]], list[tuple[re.Pattern, dict[str, str]]]]:
    rows = read_csv_rows(path)
    by_norm: dict[str, dict[str, str]] = {}
    patterns: list[tuple[re.Pattern, dict[str, str]]] = []
    for row in rows:
        norm = normalize_name(row.get("drug_normalized") or row.get("drug_display"))
        if norm:
            by_norm[norm] = row
        synonyms = row.get("synonyms", "")
        if synonyms:
            patterns.append((synonym_regex(synonyms), row))
    return by_norm, patterns


def classify_chemical(name: str, lex_by_norm: dict[str, dict[str, str]], lex_patterns: list[tuple[re.Pattern, dict[str, str]]]) -> dict[str, str]:
    norm = normalize_name(name)
    if norm in lex_by_norm:
        row = lex_by_norm[norm]
        return {
            "drug_normalized": row.get("drug_normalized", norm),
            "drug_display": row.get("drug_display", name),
            "therapy_class": row.get("therapy_class", "mapped_local_lexicon"),
            "specific_drug_yes_no": row.get("is_specific_drug_yes_no", "yes"),
            "generic_therapy_term_yes_no": row.get("is_generic_therapy_term_yes_no", "no"),
            "mapping_confidence": row.get("mapping_confidence", "high"),
            "local_lexicon_match_yes_no": "yes",
        }
    for pattern, row in lex_patterns:
        if pattern.search(name):
            return {
                "drug_normalized": row.get("drug_normalized", norm),
                "drug_display": row.get("drug_display", name),
                "therapy_class": row.get("therapy_class", "mapped_local_lexicon"),
                "specific_drug_yes_no": row.get("is_specific_drug_yes_no", "yes"),
                "generic_therapy_term_yes_no": row.get("is_generic_therapy_term_yes_no", "no"),
                "mapping_confidence": row.get("mapping_confidence", "high"),
                "local_lexicon_match_yes_no": "yes",
            }
    return {
        "drug_normalized": norm,
        "drug_display": name,
        "therapy_class": "unmapped_chemical",
        "specific_drug_yes_no": "needs_review",
        "generic_therapy_term_yes_no": "no",
        "mapping_confidence": "needs_review",
        "local_lexicon_match_yes_no": "no",
    }


def sentence_split(text: str) -> list[str]:
    text = normalize_space(text)
    if not text:
        return []
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if len(s.strip()) > 20]


def classify_relation(sentence: str) -> tuple[str, str, str]:
    for category, pattern in RELATION_CUES:
        match = re.search(pattern, sentence, flags=re.I)
        if match:
            neg = "yes" if NEGATION_RE.search(sentence) else "no"
            return category, match.group(0), neg
    return "unclear_co_mention", "", "yes" if NEGATION_RE.search(sentence) else "no"


def death_pattern(mode: str) -> re.Pattern:
    key = normalize_space(mode).lower()
    return re.compile(DEATH_VARIANTS.get(key, re.escape(key)), re.I)


def token_distance(sentence: str, drug: str, death_re: re.Pattern) -> str:
    lower = sentence.lower()
    drug_match = re.search(re.escape(drug), lower, flags=re.I)
    death_match = death_re.search(sentence)
    if not drug_match or not death_match:
        return ""
    prefix_drug = lower[: drug_match.start()]
    prefix_death = lower[: death_match.start()]
    return str(abs(len(prefix_drug.split()) - len(prefix_death.split())))


def relation_confidence(category: str, negation: str, distance: str) -> str:
    if category in {"induction_or_activation", "sensitization"} and negation == "no":
        try:
            if distance and int(distance) <= 18:
                return "high"
        except ValueError:
            pass
    if category != "unclear_co_mention" and negation == "no":
        return "medium"
    return "low"


def build_local_lexicon_fallback(local_mentions_path: Path) -> list[dict[str, str]]:
    rows = read_csv_rows(local_mentions_path)
    fallback = []
    for row in rows:
        fallback.append(
            {
                "pmid": row.get("pmid", ""),
                "drug_normalized": row.get("drug_normalized", ""),
                "drug_display": row.get("drug_display", row.get("drug_normalized", "")),
                "therapy_class": row.get("therapy_class", ""),
                "specific_drug_yes_no": row.get("is_specific_drug_yes_no", ""),
                "generic_therapy_term_yes_no": row.get("is_generic_therapy_term_yes_no", ""),
                "mapping_confidence": row.get("mapping_confidence", ""),
                "evidence_sources": "local_title_abstract_lexicon",
                "evidence_tier": "Tier3_local_lexicon_fallback",
                "mesh_ui": "",
                "registry_number": "",
                "pubtator_identifier": "",
                "mention_text": "",
                "local_lexicon_match_yes_no": "yes",
            }
        )
    return fallback


def merge_drug_evidence(
    article_rows: list[dict[str, str]],
    mesh_rows: list[dict[str, str]],
    pubtator_rows: list[dict[str, str]],
    local_rows: list[dict[str, str]],
    lex_by_norm: dict[str, dict[str, str]],
    lex_patterns: list[tuple[re.Pattern, dict[str, str]]],
) -> list[dict[str, str]]:
    article_by_pmid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in article_rows:
        article_by_pmid[str(row.get("pmid", ""))].append(row)

    merged_by_key: dict[tuple[str, str], dict[str, str]] = {}

    def add_record(pmid: str, chemical_name: str, source: str, extra: dict[str, str]) -> None:
        if not pmid or not chemical_name:
            return
        cls = classify_chemical(chemical_name, lex_by_norm, lex_patterns)
        key = (pmid, cls["drug_normalized"])
        rec = merged_by_key.setdefault(
            key,
            {
                "pmid": pmid,
                **cls,
                "evidence_sources": "",
                "evidence_tier": "",
                "mesh_ui": "",
                "registry_number": "",
                "pubtator_identifier": "",
                "mention_text": "",
            },
        )
        sources = set(filter(None, rec["evidence_sources"].split(";")))
        sources.add(source)
        rec["evidence_sources"] = ";".join(sorted(sources))
        tiers = set(filter(None, rec["evidence_tier"].split(";")))
        tiers.add("Tier1_PubMed_MeSH_ChemicalList" if source == "PubMed_MeSH_ChemicalList" else "Tier2_PubTator3_chemical")
        rec["evidence_tier"] = ";".join(sorted(tiers))
        for field in ["mesh_ui", "registry_number", "pubtator_identifier", "mention_text"]:
            val = extra.get(field, "")
            if val and val not in rec.get(field, "").split(";"):
                rec[field] = ";".join(filter(None, [rec.get(field, ""), val]))

    for row in mesh_rows:
        add_record(row.get("pmid", ""), row.get("chemical_name", ""), "PubMed_MeSH_ChemicalList", row)
    for row in pubtator_rows:
        add_record(row.get("pmid", ""), row.get("chemical_name", ""), "PubTator3_chemical", row)
    for row in local_rows:
        pmid = row.get("pmid", "")
        norm = row.get("drug_normalized", "")
        key = (pmid, norm)
        rec = merged_by_key.setdefault(key, {"pmid": pmid, **row})
        sources = set(filter(None, rec.get("evidence_sources", "").split(";")))
        sources.add("local_title_abstract_lexicon")
        rec["evidence_sources"] = ";".join(sorted(sources))
        tiers = set(filter(None, rec.get("evidence_tier", "").split(";")))
        tiers.add("Tier3_local_lexicon_fallback")
        rec["evidence_tier"] = ";".join(sorted(tiers))

    output: list[dict[str, str]] = []
    for (pmid, _drug), rec in merged_by_key.items():
        contexts = article_by_pmid.get(pmid, [{}])
        for ctx in contexts:
            out = dict(rec)
            out.update(
                {
                    "tumor_family": ctx.get("tumor_family", ""),
                    "disease_term": ctx.get("disease_term", ""),
                    "death_mode": ctx.get("death_mode", ""),
                    "year": ctx.get("publication_year_num") or ctx.get("year", ""),
                    "journal": ctx.get("journal", ""),
                    "sentence_available_yes_no": "yes" if ctx.get("title") or ctx.get("abstract") else "no",
                }
            )
            output.append(out)
    return output


def build_relation_candidates(article_rows: list[dict[str, str]], merged_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    article_by_pmid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in article_rows:
        article_by_pmid[str(row.get("pmid", ""))].append(row)

    rels: list[dict[str, str]] = []
    seen = set()
    for rec in merged_rows:
        if rec.get("generic_therapy_term_yes_no") == "yes":
            continue
        pmid = rec.get("pmid", "")
        drug_display = rec.get("drug_display", "")
        if not pmid or not drug_display:
            continue
        drug_re = boundary_regex(drug_display)
        for ctx in article_by_pmid.get(pmid, []):
            mode = ctx.get("death_mode", "")
            death_re = death_pattern(mode)
            text = normalize_space((ctx.get("title", "") + ". " + ctx.get("abstract", "")).strip())
            for sent in sentence_split(text):
                if not drug_re.search(sent) or not death_re.search(sent):
                    continue
                category, cue, neg = classify_relation(sent)
                dist = token_distance(sent, drug_display, death_re)
                conf = relation_confidence(category, neg, dist)
                key = (pmid, rec.get("drug_normalized", ""), mode, sent)
                if key in seen:
                    continue
                seen.add(key)
                rels.append(
                    {
                        "pmid": pmid,
                        "year": ctx.get("publication_year_num") or rec.get("year", ""),
                        "tumor_family": ctx.get("tumor_family", rec.get("tumor_family", "")),
                        "disease_term": ctx.get("disease_term", rec.get("disease_term", "")),
                        "death_mode": mode,
                        "drug_normalized": rec.get("drug_normalized", ""),
                        "drug_display": drug_display,
                        "sentence": sent,
                        "relation_category_rule_based": category,
                        "cue_word": cue,
                        "drug_death_distance_tokens": dist,
                        "tumor_context_present_yes_no": "yes" if TUMOR_CONTEXT_RE.search(text) else "no",
                        "negation_flag": neg,
                        "confidence_rule_based": conf,
                        "evidence_sources": rec.get("evidence_sources", ""),
                        "manual_label": "",
                        "curator_notes": "",
                    }
                )
    return rels


def manual_sample(rows: list[dict[str, str]], n: int = 300) -> list[dict[str, str]]:
    rng = random.Random(20260515)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row.get("relation_category_rule_based", ""), row.get("death_mode", ""))].append(row)
    sample: list[dict[str, str]] = []
    for group_rows in grouped.values():
        rng.shuffle(group_rows)
        sample.extend(group_rows[:3])
    rng.shuffle(sample)
    return sample[:n]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--article-input", type=Path, default=DEFAULT_ARTICLE_INPUT)
    parser.add_argument("--local-lexicon", type=Path, default=DEFAULT_LOCAL_LEXICON)
    parser.add_argument("--local-mentions", type=Path, default=DEFAULT_LOCAL_MENTIONS)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    parser.add_argument("--tool", default="cell_death_atlas_oncology_drug_annotation")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--max-pmids", type=int, default=0, help="Limit PMIDs for testing; 0 means no limit.")
    parser.add_argument("--sleep-seconds", type=float, default=0.35)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and planned batches without network calls. This is also the default if --run/--run-all/--merge-only is absent.")
    parser.add_argument("--run", action="store_true", help="Perform network calls for selected PMID set.")
    parser.add_argument("--run-all", action="store_true", help="Perform network calls for all PMIDs. Alias for --run with no max limit.")
    parser.add_argument("--merge-only", action="store_true", help="Skip network and parse cached responses only.")
    parser.add_argument("--checkpoint-every", type=int, default=10, help="Print and write checkpoint summary every N batches per source.")
    args = parser.parse_args()

    OUT_DRUG.mkdir(parents=True, exist_ok=True)
    OUT_REL.mkdir(parents=True, exist_ok=True)
    OUT_QC.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    article_rows = read_csv_rows(args.article_input)
    if not article_rows:
        write_text(OUT_QC / "pubmed_pubtator_missing_input_report.md", [
            "# Missing Input Report",
            "",
            f"Article input not found or empty: {args.article_input}",
            "Run the oncology R pipeline first to create oncology_article_records_scoped.csv.",
        ])
        return 1

    pmids = sorted({str(row.get("pmid", "")).strip() for row in article_rows if str(row.get("pmid", "")).strip()})
    if args.max_pmids and not args.run_all:
        pmids = pmids[: args.max_pmids]

    total_batches = (len(pmids) + args.batch_size - 1) // args.batch_size if args.batch_size else 0
    progress("PubMed/PubTator oncology drug annotation started")
    progress(f"PMIDs considered: {len(pmids)}")
    progress(f"Batch size: {args.batch_size}; batches per source: {total_batches}")
    progress(f"Mode: {'run-all' if args.run_all else 'run' if args.run else 'merge-only' if args.merge_only else 'dry-run'}")
    progress(f"Cache directory: {args.cache_dir}")

    lex_by_norm, lex_patterns = load_local_lexicon(args.local_lexicon)
    pmid_set = set(pmids)
    local_fallback_rows = [row for row in build_local_lexicon_fallback(args.local_mentions) if row.get("pmid", "") in pmid_set]

    query_log: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    mesh_rows: list[dict[str, str]] = []
    pubtator_rows: list[dict[str, str]] = []

    dry_run = args.dry_run or not (args.run or args.run_all or args.merge_only)
    output_suffix = ""
    if dry_run:
        output_suffix = f"_dry_run_test{len(pmids)}"
    elif args.max_pmids and not args.run_all:
        output_suffix = f"_test{len(pmids)}"

    def drug_out(stem: str) -> Path:
        return OUT_DRUG / f"{stem}{output_suffix}.csv"

    def rel_out(stem: str) -> Path:
        return OUT_REL / f"{stem}{output_suffix}.csv"
    if dry_run:
        write_text(OUT_QC / "pubmed_pubtator_dry_run.md", [
            "# PubMed/PubTator Drug Annotation Dry Run",
            "",
            f"PMIDs available: {len(pmids)}",
            f"Batch size: {args.batch_size}",
            f"Cache directory: {args.cache_dir}",
            f"Output suffix: {output_suffix}",
            "",
            "No network calls were made. Re-run with --run --email YOUR_EMAIL or --run-all --email YOUR_EMAIL.",
        ])
        progress("Dry-run mode: no network calls will be made")

    for source_name, fetcher, parser_func, ext in [
        ("PubMed_EFetch", fetch_pubmed_efetch, parse_pubmed_chemical_xml, "xml"),
        ("PubTator3", fetch_pubtator, parse_pubtator_biocjson, "json"),
    ]:
        source_parsed_total = 0
        source_error_total = 0
        progress(f"{source_name}: starting {total_batches} batches")
        for batch_index, pmid_batch in enumerate(batched(pmids, args.batch_size), start=1):
            batch_start = time.time()
            if args.merge_only:
                cache = args.cache_dir / ("pubmed_efetch_xml" if source_name == "PubMed_EFetch" else "pubtator_biocjson") / f"{batch_hash('efetch' if source_name == 'PubMed_EFetch' else 'pubtator', pmid_batch)}.{ext}"
                result = HttpResult(cache.exists(), cache.read_text(encoding="utf-8") if cache.exists() else "", "cache", cache, "" if cache.exists() else "cache missing", "cached" if cache.exists() else "missing")
            else:
                result = fetcher(pmid_batch, args)
            query_log.append(
                {
                    "source": source_name,
                    "batch_index": batch_index,
                    "pmid_count": len(pmid_batch),
                    "status": result.status,
                    "cache_path": str(result.cache_path),
                    "error": result.error,
                }
            )
            if result.ok and result.text:
                try:
                    parsed = parser_func(result.text)
                    source_parsed_total += len(parsed)
                    if source_name == "PubMed_EFetch":
                        mesh_rows.extend(parsed)
                    else:
                        pubtator_rows.extend(parsed)
                    elapsed = time.time() - batch_start
                    progress(
                        f"{source_name}: batch {batch_index}/{total_batches} "
                        f"pmids={len(pmid_batch)} status={result.status} source={result.source} "
                        f"parsed={len(parsed)} total_parsed={source_parsed_total} "
                        f"elapsed={elapsed:.1f}s cache={result.cache_path}"
                    )
                except Exception as exc:  # noqa: BLE001
                    source_error_total += 1
                    errors.append({"source": source_name, "batch_index": str(batch_index), "error": f"parse error: {exc}"})
                    progress(f"{source_name}: batch {batch_index}/{total_batches} PARSE_ERROR {exc}")
            elif result.error and result.status != "dry_run":
                source_error_total += 1
                errors.append({"source": source_name, "batch_index": str(batch_index), "error": result.error})
                progress(f"{source_name}: batch {batch_index}/{total_batches} ERROR {result.error}")
            else:
                elapsed = time.time() - batch_start
                progress(
                    f"{source_name}: batch {batch_index}/{total_batches} "
                    f"pmids={len(pmid_batch)} status={result.status} source={result.source} "
                    f"elapsed={elapsed:.1f}s cache={result.cache_path}"
                )

            if args.checkpoint_every > 0 and batch_index % args.checkpoint_every == 0:
                checkpoint = [
                    "# PubMed/PubTator Checkpoint",
                    "",
                    f"- Source: {source_name}",
                    f"- Completed batches: {batch_index}/{total_batches}",
                    f"- Parsed annotations for this source: {source_parsed_total}",
                    f"- Errors for this source: {source_error_total}",
                    f"- Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}",
                ]
                checkpoint_path = OUT_QC / f"pubmed_pubtator_checkpoint_{source_name}.md"
                write_text(checkpoint_path, checkpoint)
                progress(f"{source_name}: checkpoint written to {checkpoint_path}")

        progress(f"{source_name}: completed; parsed={source_parsed_total}; errors={source_error_total}")

    mesh_fields = ["pmid", "chemical_name", "chemical_name_normalized", "mesh_ui", "registry_number", "source"]
    pubtator_fields = ["pmid", "chemical_name", "chemical_name_normalized", "pubtator_identifier", "pubtator_type", "mention_text", "mention_section", "start_offset", "end_offset", "source"]
    write_csv_rows(drug_out("pubmed_mesh_chemical_annotations"), mesh_rows, mesh_fields)
    write_csv_rows(drug_out("pubtator_chemical_annotations"), pubtator_rows, pubtator_fields)
    write_csv_rows(drug_out("pubmed_pubtator_query_log"), query_log, ["source", "batch_index", "pmid_count", "status", "cache_path", "error"])
    write_csv_rows(drug_out("pubmed_pubtator_errors"), errors, ["source", "batch_index", "error"])

    merged = merge_drug_evidence(article_rows, mesh_rows, pubtator_rows, local_fallback_rows, lex_by_norm, lex_patterns)
    merged_fields = [
        "pmid", "year", "tumor_family", "disease_term", "death_mode", "drug_normalized", "drug_display",
        "therapy_class", "specific_drug_yes_no", "generic_therapy_term_yes_no", "mapping_confidence",
        "local_lexicon_match_yes_no", "evidence_tier", "evidence_sources", "mesh_ui", "registry_number",
        "pubtator_identifier", "mention_text", "journal", "sentence_available_yes_no",
    ]
    write_csv_rows(drug_out("oncology_drug_mentions_evidence_merged"), merged, merged_fields)

    relations = build_relation_candidates(article_rows, merged)
    rel_fields = [
        "pmid", "year", "tumor_family", "disease_term", "death_mode", "drug_normalized", "drug_display",
        "sentence", "relation_category_rule_based", "cue_word", "drug_death_distance_tokens",
        "tumor_context_present_yes_no", "negation_flag", "confidence_rule_based", "evidence_sources",
        "manual_label", "curator_notes",
    ]
    write_csv_rows(rel_out("drug_death_relation_sentences_validated_input"), relations, rel_fields)
    write_csv_rows(rel_out("drug_death_relation_manual_validation_sample_pubmed_pubtator"), manual_sample(relations), rel_fields)

    summary = [
        "# PubMed/PubTator Drug Annotation Summary",
        "",
        f"- Article input: {args.article_input}",
        f"- Unique PMIDs considered: {len(pmids)}",
        f"- PubMed MeSH chemical annotations parsed: {len(mesh_rows)}",
        f"- PubTator chemical annotations parsed: {len(pubtator_rows)}",
        f"- Local fallback mentions loaded: {len(local_fallback_rows)}",
        f"- Merged PMID-drug-context rows: {len(merged)}",
        f"- Sentence relation candidates: {len(relations)}",
        f"- Network mode: {'run' if (args.run or args.run_all) else 'merge-only' if args.merge_only else 'dry-run'}",
        f"- Output suffix: {output_suffix or 'canonical'}",
        "",
        "Interpretation boundary: merged annotations identify chemical/drug evidence in article records. They do not validate drug-induced RCD mechanisms without manual sentence review.",
    ]
    summary_name = f"pubmed_pubtator_drug_annotation_summary{output_suffix}.md"
    write_text(OUT_QC / summary_name, summary)
    progress("PubMed/PubTator oncology drug annotation completed")
    print("\n".join(summary), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
