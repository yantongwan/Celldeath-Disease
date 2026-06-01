# GitHub Release Audit

Date: 2026-06-01

Scope: `CellDeath_Figure`

## Release Summary

- The folder is not currently initialized as a Git repository.
- Code inventory: 44 Python/R source files, 12,125 source lines.
- Main code modules:
  - `code/01_*_v5`: PubMed atlas construction and main/supplementary figure regeneration.
  - `code/02_gene_v3`: gene-level and comorbidity figure/table workflows.
  - `code/03_oncology_therapeutic_landscape`: oncology, drug, PubMed/PubTator, and therapeutic landscape workflows.
- Updated public scope: eight selected regulated cell death modes: ferroptosis, pyroptosis, necroptosis, NETosis, immunogenic cell death, cuproptosis, PANoptosis, and disulfidptosis.

## Findings Before Public Upload

1. Large files should not be committed to ordinary GitHub.
   - Several raw or intermediate CSV files exceed GitHub's 100 MB file limit.
   - Recommended handling: deposit raw and large intermediate source-data tables in Zenodo, figshare, OSF, or Git LFS, then cite the DOI in Data Availability.

2. macOS resource-fork files are present.
   - 39 `._*` files were detected.
   - These should be excluded by `.gitignore` and optionally deleted before the first commit.

3. No hard-coded private credentials were detected in the inspected code.
   - PubMed/NCBI access is configured through environment variables such as `NCBI_EMAIL`, `NCBI_API_KEY`, `ENTREZ_EMAIL`, and `ENTREZ_API_KEY`.
   - Keep `.env` files excluded and do not commit real API keys or personal email addresses.

4. Absolute local paths were replaced in the release code.
   - Public scripts now use configurable roots such as `CELLDEATH_ATLAS_ROOT`, `CELLDEATH_ATLAS_PACKAGE_DIR`, and `CELLDEATH_FIGURE_ROOT`.

5. Some internal file and folder labels still reflect an older screening branch.
   - This conflicts with the updated eight-mode public screening language.
   - Recommended handling: keep the public README, Methods, Data Availability, and figure documentation aligned to the eight selected modes. If time permits, rename public release folders and script stems to neutral names in a clean release branch.

6. Runtime environments are not yet documented.
   - Add a short `README.md` plus either `requirements.txt` / `environment.yml` for Python and a package list/session information block for R.

## Recommended Repository Contents

Commit to GitHub:

- `code/`
- `.gitignore`
- `README.md` or release README draft
- `DATA_AVAILABILITY_STATEMENT.md`
- `GITHUB_UPLOAD_COMMANDS.md`
- `NATURE_METHODS_ALL_PNGS.md`
- `FIGURE_METRIC_AXIS_DEFINITIONS.md`, if this is the latest metric-axis documentation
- `table/00_manifest/`
- selected small figure source tables needed for peer review, if they are not redundant with the DOI archive

Archive externally with a DOI:

- `raw/`
- large intermediate CSV files
- large drug-reference lookup tables
- editable Illustrator/TIFF/ZIP figure bundles, unless the journal specifically requires them in a source-data package

## Minimum Pre-Release Fixes

1. Keep the eight selected death modes explicit in all public text.
2. Exclude raw and very large intermediate data from ordinary Git history.
3. Remove or ignore macOS resource-fork files.
4. Document required Python/R packages and environment variables.
5. Add a DOI-backed source-data archive before manuscript submission or final public release.
