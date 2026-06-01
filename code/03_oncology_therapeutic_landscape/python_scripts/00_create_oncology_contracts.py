#!/usr/bin/env python3
"""Create evidence-gated contracts for the oncology therapeutic RCD module."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_DIR = ROOT / "00_analysis_contracts"


COMMON_BOUNDARIES = """\
Allowed interpretation:
- Literature-level oncology signal in PubMed-indexed article records.
- Tumor-specific death-mode literature profile or therapy-linked RCD vocabulary.
- Hypothesis-generating triads requiring sentence-level and manual validation before mechanistic claims.

Forbidden interpretation:
- Do not claim causal death-pathway induction from co-mention alone.
- Do not claim validated mechanism, clinical efficacy, therapeutic maturity, or true tumor vulnerability.
- Do not turn sentence-level candidates into experimental proof.
"""


CONTRACTS = {
    "oncology_module_master_contract.md": {
        "question": "Can an oncology-focused, therapy-aware module deepen the RCD atlas without overclaiming drug mechanism or clinical readiness?",
        "inputs": "Audited pair metrics, pair-PMID links, disease-system mapping, oncology term classification, PMID metadata, cleaned vocabulary tables, low-count warnings, and available title/abstract/MeSH chemical fields.",
        "evidence": "Tiered: disease/pair-level oncology evidence; article-level metadata; MeSH/title/abstract drug mentions; sentence-level relation candidates; blank manual-validation sample.",
        "confounders": "Oncology publication volume, tumor-name synonymy, generic therapy words, tool compounds, uneven death-mode ages, low-count triads, and abstract-only relation ambiguity.",
        "normalization": "Tumor-normalized shares, oncology-background selectivity, death-mode breadth, recency components, and transparent component-wise heat index.",
        "gate": "Main-text eligible only if oncology death-mode and tumor-profile gates pass and therapy relation evidence clears sentence-level thresholds. Otherwise supplementary/table-only.",
        "table": "Module-level gate summary and integration decision.",
        "figure": "Only generate figures whose module gates pass; reject or demote therapy figures if relation evidence is insufficient.",
        "placement": "Main text only for robust tumor-profile or relation-gated therapy results; otherwise supplementary diagnostics.",
    },
    "death_mode_hotness_contract.md": {
        "question": "Which RCD-associated death-mode concepts are currently hottest in oncology literature?",
        "inputs": "Oncology-filtered pair-PMID records with year, death mode, tumor family, low-count flags, and therapy/drug vocabulary indicators if available.",
        "evidence": "Article-level oncology counts, tumor-family breadth, recent growth/burst, therapy vocabulary fraction, and drug-associated PMID fraction.",
        "confounders": "Ferroptosis or other high-volume modes may dominate by raw count; newer terms may appear hot because denominator is small.",
        "normalization": "Report components separately and compute a transparent z-score heat index; test Pair_Count >=3 sensitivity when pair counts are available.",
        "gate": "At least five death modes with oncology_unique_pmids >= 50, ranking not collapsed by one tumor family >80% unless disclosed, and Pair_Count >=3 sensitivity does not invert the main message.",
        "table": "03_oncology_death_mode_hotness/oncology_death_mode_hotness.csv",
        "figure": "Figure_oncology_death_mode_hotness if gate passes.",
        "placement": "Supplementary or main-text support depending on strength; never framed as therapeutic target ranking.",
    },
    "tumor_specificity_contract.md": {
        "question": "Do tumor families show different death-mode literature profiles after oncology-background normalization?",
        "inputs": "Oncology tumor mapping, pair-PMID links, pair counts, years, low-count flags, and disease/tumor family mappings.",
        "evidence": "Tumor family x death mode pair counts, unique PMID support, tumor-normalized share, oncology-background selectivity, recent growth, and low-count flags.",
        "confounders": "Tumor family volume, broad versus specific tumor MeSH terms, over-representation of generic neoplasm records, and low-support tumor subtypes.",
        "normalization": "Compute tumor-normalized shares and log2 oncology-background selectivity with displayed count support.",
        "gate": "Displayed tumor families require unique_pmid_count >= 20 or total_pair_count >= 20; displayed cells require pair_count >=3 or explicit low-support marking.",
        "table": "04_tumor_specific_death_profiles/tumor_death_selectivity.csv",
        "figure": "Figure_tumor_specific_death_landscape if gate passes.",
        "placement": "Can support main-text oncology heterogeneity paragraph if readable and robust.",
    },
    "drug_extraction_contract.md": {
        "question": "Which drugs or therapy classes are linked to oncology RCD-associated article records?",
        "inputs": "PMID metadata with title/abstract and any MeSH chemical/keyword fields; oncology-filtered article records; curated therapy/drug lexicon from existing metadata plus explicit oncology therapy keywords.",
        "evidence": "Tier 1 MeSH/Supplementary Concept chemical terms, Tier 2 title/abstract lexicon matches, Tier 3 sentence-window relation candidates.",
        "confounders": "Generic words such as therapy/treatment/agent, abbreviations, natural products, tool compounds, and brand/generic synonyms.",
        "normalization": "Normalize drug names, keep generic therapy terms separate from specific drugs, and label each mention by evidence tier and mapping confidence.",
        "gate": "Drug extraction passes if metadata text or chemical fields exist and specific drug mentions are not dominated by a single generic term. Relation plotting requires the sentence-level gate.",
        "table": "05_drug_therapy_extraction/oncology_drug_mentions.csv",
        "figure": "No standalone figure unless downstream relation/triad gates pass.",
        "placement": "Source-data and supplementary diagnostics unless relation evidence is strong.",
    },
    "drug_death_triads_contract.md": {
        "question": "Which tumor-drug-death-mode triads have enough article and relation support to be hypothesis-generating?",
        "inputs": "Oncology drug mentions, sentence-level relation candidates, tumor mapping, death modes, years, and low-count flags.",
        "evidence": "PMID count, high-confidence relation count, relation categories, recency, specificity, and evidence-tier summary.",
        "confounders": "Tool-compound literatures, generic tumor records, same abstract mentioning multiple drugs/death modes, and unclear relation direction.",
        "normalization": "Grade triads by transparent thresholds; separate unclear co-mentions from induction/sensitization candidates.",
        "gate": "Strong/moderate triads require pmid_count >=5, high_confidence_relation_count >=2, non-generic drug name, and relation categories not entirely unclear.",
        "table": "06_drug_death_mode_triads/oncology_drug_death_triads.csv",
        "figure": "Triad figure only if enough strong/moderate triads are present and readable.",
        "placement": "Supplementary unless manually validated or highly robust.",
    },
    "same_drug_multi_death_contract.md": {
        "question": "Do the same specific drugs appear in multiple death-mode induction/sensitization contexts across tumor literature?",
        "inputs": "Relation-gated drug-death records with tumor family, death mode, therapy class, and evidence tier.",
        "evidence": "Drug death-mode breadth, tumor-family breadth, high-confidence relation count, death-mode entropy, and low-confidence fraction.",
        "confounders": "High-profile drugs and tool compounds may appear broadly because they are heavily studied, not because they induce multiple death pathways.",
        "normalization": "Compute a polydeath index from breadth, tumor breadth, relation support, and entropy; show components rather than relying on a single score.",
        "gate": "Display only specific drugs with death_mode_breadth >=2, high_confidence_relation_count >=3, and total_pmid_count >=5.",
        "table": "07_same_drug_multi_death_signals/drug_polydeath_signal_table.csv",
        "figure": "Figure_same_drug_multi_death_signals if gate passes.",
        "placement": "Supplementary unless relation evidence is strong and not dominated by tool compounds.",
    },
    "therapy_class_contract.md": {
        "question": "Do therapy classes differ in death-mode literature contexts and relation categories?",
        "inputs": "Drug mention and relation tables with therapy class mappings, death modes, tumor families, years, and confidence tiers.",
        "evidence": "Therapy class x death-mode counts, high-confidence relation counts, relation category composition, tumor breadth, selectivity, and recent growth.",
        "confounders": "Broad classes such as chemotherapy or targeted therapy can dominate; experimental tool compounds may not represent clinical therapy.",
        "normalization": "Class-death selectivity relative to oncology background; display support counts and mark low-count cells.",
        "gate": "Plot if at least three therapy classes and five death modes have enough high-confidence relation support; otherwise table-only.",
        "table": "08_therapy_class_landscape/therapy_class_death_mode_matrix.csv",
        "figure": "Figure_therapy_class_death_mode_landscape if gate passes.",
        "placement": "Supplementary relation-gated diagnostic unless strongly supported.",
    },
    "sentence_relation_audit_contract.md": {
        "question": "Can title/abstract sentences support induction/sensitization/inhibition relation candidates between drug and death-mode vocabulary?",
        "inputs": "Title/abstract text, tumor context, drug lexicon, death-mode term variants, sentence tokenizer, and relation cue lexicons.",
        "evidence": "Same-sentence or adjacent-sentence drug-death co-occurrence with relation cue, negation flag, token distance, tumor context, and confidence grade.",
        "confounders": "Abstract compression, negation, passive wording, multiple drugs per sentence, review articles, and generic therapy terms.",
        "normalization": "Require same-sentence cue when possible; separate high/medium/low confidence; generate a stratified manual validation sample with blank labels.",
        "gate": "At least 500 sentence-level candidates and at least 100 high-confidence induction_or_activation or sensitization records; no single generic term >50%; manual validation sample exists.",
        "table": "09_sentence_level_relation_audit/drug_death_relation_sentences.csv",
        "figure": "Sentence audit figure S19 only if relation gate passes; otherwise report as incomplete/exploratory.",
        "placement": "Always supplementary/audit-facing unless manual validation is completed.",
    },
}


def render_contract(name: str, spec: dict) -> str:
    return f"""# {name.replace('_', ' ').replace('.md', '').title()}

## Scientific Question
{spec['question']}

## Input Data Required
{spec['inputs']}

## Evidence Level
{spec['evidence']}

## Possible Confounders
{spec['confounders']}

## Normalization Strategy
{spec['normalization']}

## Gate Criteria
{spec['gate']}

## Output Table
{spec['table']}

## Candidate Figure
{spec['figure']}

## Manuscript Placement
{spec['placement']}

## Interpretation Boundaries
{COMMON_BOUNDARIES}
"""


def main() -> None:
    CONTRACT_DIR.mkdir(parents=True, exist_ok=True)
    for filename, spec in CONTRACTS.items():
        (CONTRACT_DIR / filename).write_text(render_contract(filename, spec), encoding="utf-8")
    print(f"contracts_written={len(CONTRACTS)}")
    print(f"contract_dir={CONTRACT_DIR}")


if __name__ == "__main__":
    main()
