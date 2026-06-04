# Eight-Mode Regulated Cell Death Disease Atlas

This repository provides analysis and plotting code for a PubMed/PubTator-based literature atlas of disease, gene, oncology and therapeutic contexts across eight selected regulated cell death modes.

The outputs quantify literature-level associations. Co-mentions, PubMed/PubTator annotations, sentence-window relation candidates and composite literature scores are bibliometric and hypothesis-generating signals rather than evidence of causal mechanisms, therapeutic activity or clinical efficacy.

## Death Modes

- ferroptosis
- pyroptosis
- necroptosis
- NETosis
- immunogenic cell death
- cuproptosis
- PANoptosis
- disulfidptosis

## Repository Layout

```text
code/
  01_selected_death_modes_v5/            Core PubMed atlas and figure regeneration scripts
  02_gene_v3/                            Gene/protein co-mention and VOSviewer extension scripts
  03_oncology_therapeutic_landscape/     Oncology, drug and therapeutic landscape scripts

rcd_atlas_explorer/
  app.py                                 Streamlit companion website entry point
  data/processed/rcd_atlas.parquet/      Compact public data snapshot for app deployment
  README.md                              App-specific local and deployment instructions

table/00_manifest/
  selected_death_modes_manifest.csv      Final eight-mode screening and pair-atlas summary
  source_trace_README.md                 Source-table and raw-data traceability notes

DATA_AVAILABILITY_STATEMENT.md           Data and code availability statement
FIGURE_METRIC_AXIS_DEFINITIONS.md        Figure metrics and axis definitions
NATURE_METHODS_ALL_PNGS.md               Nature-style methods summary for selected PNG figures
requirements.txt                         Python package requirements
R_PACKAGES.md                            R package requirements
```

## Data Scope

The figure package covers publication years 2000-2025 and the eight selected death-mode concepts listed above. The public manifest is:

```text
table/00_manifest/selected_death_modes_manifest.csv
```

The manifest includes PubMed term-screen counts, pair-atlas unique PMID counts and first pair-atlas publication years for the final eight concepts.

This repository tracks code, methods documentation, metric definitions and a lightweight manifest. Large raw PMID tables, PubMed/PubTator context tables, figure source-data directories and editable figure assets are excluded from Git history because of file-size and database-snapshot constraints.

## Interactive Atlas Website

The Streamlit companion website is included in:

```text
rcd_atlas_explorer/app.py
```

For local preview:

```bash
cd rcd_atlas_explorer
streamlit run app.py
```

For Streamlit Community Cloud deployment, use this repository and set the main
file path to:

```text
rcd_atlas_explorer/app.py
```

The app can run from the bundled Parquet snapshot in
`rcd_atlas_explorer/data/processed/rcd_atlas.parquet/`; it creates a temporary
DuckDB database at runtime when `rcd_atlas.duckdb` is not present. The full raw
CSV corpus and local DuckDB snapshot are intentionally excluded from ordinary
Git history.

## Methods And Metrics

Detailed methodological descriptions are provided in:

- `NATURE_METHODS_ALL_PNGS.md`
- `FIGURE_METRIC_AXIS_DEFINITIONS.md`
- `table/00_manifest/source_trace_README.md`

## Environment

Python package requirements are listed in `requirements.txt`.

R package requirements are listed in `R_PACKAGES.md`.
