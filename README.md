# Eight-Mode Regulated Cell Death Disease Atlas

This repository contains the analysis and plotting code for a PubMed/PubTator-based literature atlas of disease, gene, oncology and therapeutic contexts across eight selected regulated cell death modes:

- ferroptosis
- pyroptosis
- necroptosis
- NETosis
- immunogenic cell death
- cuproptosis
- PANoptosis
- disulfidptosis

The analyses quantify literature-level associations. Co-mentions, PubMed/PubTator annotations, sentence-window relation candidates and composite literature scores should be interpreted as bibliometric and hypothesis-generating signals, not as proof of causal mechanisms, therapeutic activity or clinical efficacy.

## Repository Layout

```text
code/
  01_selected_death_modes_v5/            Core PubMed atlas and figure regeneration scripts
  02_gene_v3/                            Gene/protein co-mention and VOSviewer extension scripts
  03_oncology_therapeutic_landscape/     Oncology, drug and therapeutic landscape scripts

table/00_manifest/
  selected_death_modes_manifest.csv      Final eight-mode screening and pair-atlas summary
  source_trace_README.md                 Source-table and raw-data traceability notes

DATA_AVAILABILITY_STATEMENT.md           Data and code availability draft
FIGURE_METRIC_AXIS_DEFINITIONS.md        Figure metrics and axis definitions
NATURE_METHODS_ALL_PNGS.md               Nature-style methods summary for selected PNG figures
GITHUB_RELEASE_AUDIT.md                  Release audit notes
GITHUB_UPLOAD_COMMANDS.md                Upload command template
requirements.txt                         Python package requirements
R_PACKAGES.md                            R package requirements
```

Large raw PMID tables, PubMed/PubTator context tables, figure source-data directories and editable figure assets are intentionally excluded from ordinary Git history. They should be released through a DOI-backed archive such as Zenodo, figshare or OSF.

## Data Scope

The final figure package uses publication years 2000-2025 and the eight selected death-mode concepts listed above. The public manifest is:

```text
table/00_manifest/selected_death_modes_manifest.csv
```

The manifest includes PubMed term-screen counts, pair-atlas unique PMID counts and first pair-atlas publication years for the final eight concepts.

## Environment

Python scripts were written for Python 3 and use common scientific packages listed in `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

R scripts require packages listed in `R_PACKAGES.md`. Bioconductor packages are used for the gene/network extension module.

## Path Configuration

The public code avoids hard-coded local disk paths. Set these environment variables when reproducing the full analysis from archived data:

```bash
export CELLDEATH_ATLAS_ROOT=/path/to/archived/analysis-root
export CELLDEATH_ATLAS_PACKAGE_DIR=/path/to/selected_death_modes_2000_2025
export CELLDEATH_FIGURE_ROOT=/path/to/CellDeath_Figure
```

Optional variables used by some upstream-regeneration scripts:

```bash
export CELLDEATH_ATLAS_PACKAGE=selected_death_modes_2000_2025
export CELLDEATH_ATLAS_V4_DIR=/path/to/v4_regenerated_source_package
export CELLDEATH_PAIR_METRICS_BASE=/path/to/pair_metrics_2025_cutoff.csv
export NCBI_EMAIL=your_email@example.com
export NCBI_API_KEY=your_ncbi_api_key
export ENTREZ_EMAIL=your_email@example.com
export ENTREZ_API_KEY=your_ncbi_api_key
```

Do not commit real API keys, private email settings or `.env` files.

## Typical Use

Run a module after setting the data paths:

```bash
python code/01_selected_death_modes_v5/regenerate_figure2B_corrected_top15_disease_hubs.py
Rscript code/02_gene_v3/run_all_gene_v3_figures.R "$CELLDEATH_ATLAS_ROOT"
Rscript code/03_oncology_therapeutic_landscape/R_scripts/run_all_oncology_therapeutic_landscape.R
```

Some scripts regenerate complete upstream packages, while others redraw a specific figure panel. Check the script name and header before running.

## Data Availability

The intended public release model is:

- GitHub: code, lightweight manifest, methods documentation and metric definitions.
- DOI archive: raw PMID lists, PMID-level context tables, PubMed/PubTator intermediate outputs, large drug-reference tables and final figure source-data bundles.
- Manuscript source data: final tables required for direct figure inspection, as requested by the journal.

Suggested manuscript wording is provided in `DATA_AVAILABILITY_STATEMENT.md`.

## Citation

If using this code before publication, cite the GitHub repository and the corresponding DOI archive once available.
