# Expected Raw Files

Place manuscript-derived CSV, TSV, or XLSX files directly in this directory. Files in
`data/raw/demo/` are only used when no non-demo raw files are found.

The builder detects logical tables from filename patterns and maps flexible source
columns through `config/schema_config.yaml`.

Recommended source families:

- Article metadata: `article_metadata*.csv`
- Death-mode disease article links: `pair_pmids*.csv` or `pair_pmid_links*.csv`
- Death-mode disease pair summaries: `death_disease_pairs*.csv` or `pair_metrics*.csv`
- Gene mention files: `gene_matches_accepted*.csv`, `gene_death_disease_pmid_links*.csv`
- Drug mention files: `oncology_drug_mentions*.csv`, `drug_death_relation_sentences*.csv`
- Optional disease dictionary: `disease_system_mapping*.csv`
- Optional synonyms: `synonyms*.csv`

Raw files are never edited by the build script.
