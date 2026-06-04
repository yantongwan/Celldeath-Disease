# Interactive atlas of regulated cell death across human disease

This repository builds a local Streamlit + DuckDB companion app for the
biomedical literature atlas project:

> Cartography of emerging regulated cell death across human disease

The app converts static bibliometric, molecular and pharmacological tables into
queryable article-level evidence views. It is designed for local use and does
not require an internet connection after the database has been built.

**Interpretation warning:** This application summarizes literature-level associations from PubMed-derived data. A death-mode–disease pair, gene co-mention or drug–death relation candidate does not prove biological causality, pathway activation, therapeutic efficacy or clinical validity.

## Project Structure

```text
rcd_atlas_explorer/
  app.py
  requirements.txt
  config/
  data/
    raw/
    processed/
    exports/
  scripts/
  src/
  tests/
```

## Required Data Files

Place source CSV, TSV or XLSX files directly under `data/raw/`. Files in
`data/raw/demo/` are loaded only when no non-demo raw files are found.

Recommended source families:

- Article metadata: `article_metadata*.csv`
- Death-mode disease PMID links: `pair_pmids*.csv`, `pair_pmid_links*.csv`
- Pair summaries: `death_disease_pairs*.csv`, `pair_metrics*.csv`
- Gene layers: `gene_matches_accepted*.csv`, `gene_death_disease_pmid_links*.csv`
- Drug layers: `oncology_drug_mentions*.csv`, `drug_death_relation_sentences*.csv`
- Optional dictionaries: `disease_system_mapping*.csv`, `synonyms*.csv`

Column mapping is configured in `config/schema_config.yaml`; update that file
instead of changing the builder when source column names differ.

## Installation

```bash
pip install -r requirements.txt
```

## Build The Database

```bash
python scripts/build_database.py
```

This writes:

- `data/processed/rcd_atlas.duckdb`
- `data/processed/validation_report.md`
- `data/processed/rcd_atlas.parquet/`
- `data/exports/*.csv`

Raw files are never modified.

## Run The App

```bash
streamlit run app.py
```

For the public GitHub/Streamlit deployment package, the large local DuckDB
snapshot is not required. If `data/processed/rcd_atlas.duckdb` is absent but
`data/processed/rcd_atlas.parquet/` is present, the app builds a temporary
runtime DuckDB database from the bundled Parquet tables.

## App Tabs

1. `Atlas overview`: PMID counts, death-mode breadth and summary exports.
2. `Death mode to diseases`: death-mode disease summaries and article evidence.
3. `Gene + disease`: gene, disease and death-mode queries with artefact handling.
4. `Drug / disease / death mode`: drug-disease and death-mode-disease candidate relation queries.
5. `Article explorer`: global search across article metadata and association layers.
6. `Data dictionary`: loaded files, row counts, validation warnings and definitions.

## Data Interpretation

The app reports literature-level associations. Co-occurrence, text-mined gene
mentions and drug-death relation candidates are hypothesis-generating. All
drug-death-disease relationships should be read as literature-supported
candidate associations, not manually validated mechanisms.

`ASCL4` is flagged as a likely `ACSL4` orthographic/entity-normalization
artefact in ferroptosis contexts and is excluded from biological summaries by
default. It remains visible when the UI toggle `Show flagged artefacts` is
enabled.

Negated drug relation sentences are retained in the database for transparency
and excluded by default in the UI.

## Adding New Data

1. Place new raw files in `data/raw/`.
2. Add filename patterns or aliases in `config/schema_config.yaml` if needed.
3. Rebuild with `python scripts/build_database.py`.
4. Check `data/processed/validation_report.md`.

## Local Server Or Streamlit Cloud

For local network use:

```bash
streamlit run app.py --server.address 0.0.0.0
```

For Streamlit Community Cloud, commit the repository with a built database only
if your data-sharing policy permits it and the file-size policy allows it. The
recommended public deployment path for this repository is to commit the bundled
Parquet snapshot and set the app entry point to:

```text
rcd_atlas_explorer/app.py
```

## Troubleshooting

- `Database not found`: run `python scripts/build_database.py`.
- `No rows loaded`: verify that source files are in `data/raw/`, not only in a
  sibling manuscript folder.
- `Missing columns`: update `config/schema_config.yaml` aliases.
- Slow first query: DuckDB is filtering on disk; subsequent Streamlit reruns use
  cached query results.
- Missing titles or journals: the builder does not fabricate metadata. It keeps
  PMIDs and marks missing metadata in the validation report.
