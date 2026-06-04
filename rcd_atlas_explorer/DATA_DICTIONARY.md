# Data Dictionary

## Core Article Columns

- `pmid`: PubMed identifier stored as a string.
- `title`, `journal`, `publication_year`, `publication_date`, `doi`, `abstract`: article metadata.
- `publication_types`: PubMed publication-type text when available.
- `evidence_stage`: source-level stage label. This is descriptive and not a validation grade.
- `publicationtype_supported_trial`: true when publication-type text supports trial-level evidence.
- `is_oncology`, `tumor_family`: oncology-oriented flags from source data or disease-term inference.

## Atlas Association Columns

- `death_mode`: one of the eight canonical regulated cell-death modes.
- `disease_term`: normalized disease term from source tables.
- `disease_system`: optional disease-system grouping.
- `pair_count`, `unique_pmids`: literature-volume measures for a death-mode disease pair.
- `literature_hub_score`, `selectivity_score`: manuscript-derived descriptive metrics when available.

## Gene Columns

- `gene_symbol`: uppercase gene symbol after normalization.
- `source_field`, `matched_text`, `context_sentence`: text-mining provenance fields.
- `artifact_flag`: true for entries excluded from biological summaries by default.

## Drug Columns

- `normalized_drug_name`: conservative normalized drug label.
- `reference_class_label`: source drug/reference-class label where available.
- `approved_clinical_drug_flag`, `fda_approved_flag`: source flags; absence does not imply non-approval.
- `relation_cue`, `relation_sentence`, `negation_flag`, `token_distance`: sentence-level relation provenance.

## Interpretation Boundary

This application summarizes literature-level associations from PubMed-derived data. A death-mode–disease pair, gene co-mention or drug–death relation candidate does not prove biological causality, pathway activation, therapeutic efficacy or clinical validity.
