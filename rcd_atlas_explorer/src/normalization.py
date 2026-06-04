from __future__ import annotations

import re
from typing import Any


CANONICAL_DEATH_MODES: tuple[str, ...] = (
    "Ferroptosis",
    "Pyroptosis",
    "Necroptosis",
    "NETosis",
    "Immunogenic cell death",
    "Cuproptosis",
    "PANoptosis",
    "Disulfidptosis",
)


DEATH_MODE_SYNONYMS: dict[str, str] = {
    "ferroptosis": "Ferroptosis",
    "pyroptosis": "Pyroptosis",
    "necroptosis": "Necroptosis",
    "netosis": "NETosis",
    "nets": "NETosis",
    "neutrophil extracellular traps": "NETosis",
    "icd": "Immunogenic cell death",
    "immunogenic cell death": "Immunogenic cell death",
    "cuproptosis": "Cuproptosis",
    "panoptosis": "PANoptosis",
    "disulfidptosis": "Disulfidptosis",
    "disulfidptosis like": "Disulfidptosis",
}


TRUE_VALUES = {"1", "true", "t", "yes", "y", "positive", "present"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "negative", "absent"}


def _to_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null", "<na>"}:
        return ""
    return text


def normalize_pmid(value: Any) -> str:
    """Return a PMID as a clean string, preserving empty values as empty strings."""
    text = _to_text(value)
    if not text:
        return ""
    text = re.sub(r"\.0+$", "", text)
    match = re.search(r"\d+", text)
    return match.group(0) if match else ""


def _death_key(value: Any) -> str:
    text = _to_text(value)
    text = re.sub(r"[\s_\-/]+", " ", text)
    text = re.sub(r"[^A-Za-z0-9 ]+", "", text)
    return text.strip().lower()


def normalize_death_mode(value: Any) -> str:
    key = _death_key(value)
    if not key:
        return ""
    return DEATH_MODE_SYNONYMS.get(key, _to_text(value).strip())


def is_canonical_death_mode(value: Any) -> bool:
    return normalize_death_mode(value) in CANONICAL_DEATH_MODES


def normalize_disease(value: Any) -> str:
    text = _to_text(value)
    text = re.sub(r"\s+", " ", text)
    text = text.strip(" \t\r\n;,.")
    return text


def normalize_gene_symbol(value: Any) -> str:
    text = _to_text(value)
    text = re.sub(r"\s+", "", text)
    return text.upper()


def normalize_drug_name(value: Any) -> str:
    text = _to_text(value)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" \t\r\n;,.").lower()


def as_bool(value: Any) -> bool | None:
    text = _to_text(value).lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    return None


def as_int(value: Any) -> int | None:
    text = _to_text(value)
    if not text:
        return None
    match = re.search(r"-?\d+", text)
    return int(match.group(0)) if match else None


def as_float(value: Any) -> float | None:
    text = _to_text(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def infer_oncology_flag(disease_term: Any, explicit: Any = None) -> bool | None:
    explicit_bool = as_bool(explicit)
    if explicit_bool is not None:
        return explicit_bool
    disease = normalize_disease(disease_term).lower()
    if not disease:
        return None
    oncology_tokens = (
        "neoplasm",
        "cancer",
        "tumor",
        "tumour",
        "carcinoma",
        "sarcoma",
        "leukemia",
        "leukaemia",
        "lymphoma",
        "melanoma",
        "glioma",
        "blastoma",
    )
    return any(token in disease for token in oncology_tokens)


def flag_gene_artifact(gene_symbol: Any, death_mode: Any = None) -> bool:
    symbol = normalize_gene_symbol(gene_symbol)
    mode = normalize_death_mode(death_mode)
    return symbol == "ASCL4" and (not mode or mode == "Ferroptosis")


def make_pubmed_link(pmid: str) -> str:
    clean = normalize_pmid(pmid)
    return f"https://pubmed.ncbi.nlm.nih.gov/{clean}/" if clean else ""


def make_doi_link(doi: Any) -> str:
    text = _to_text(doi)
    if not text:
        return ""
    if text.startswith("http://") or text.startswith("https://"):
        return text
    return f"https://doi.org/{text}"


def abstract_snippet(abstract: Any, max_chars: int = 280) -> str:
    text = re.sub(r"\s+", " ", _to_text(abstract))
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 1].rstrip() + "..."
