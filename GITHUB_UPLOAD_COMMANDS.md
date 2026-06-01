# GitHub Upload Commands

Run these commands from the project root:

```bash
cd CellDeath_Figure
```

## 1. Inspect Files Before Committing

```bash
git status --short
find . -name '._*' -type f -print
find . -type f -size +90M -print
```

Optional cleanup for macOS resource-fork files:

```bash
find . -name '._*' -type f -delete
```

## 2. Initialize Git

```bash
git init
git branch -M main
```

## 3. Add Only Code and Lightweight Documentation

Recommended first public commit:

```bash
git add .gitignore
git add README.md
git add requirements.txt
git add R_PACKAGES.md
git add code
git add DATA_AVAILABILITY_STATEMENT.md
git add GITHUB_RELEASE_AUDIT.md
git add GITHUB_UPLOAD_COMMANDS.md
git add NATURE_METHODS_ALL_PNGS.md
git add FIGURE_METRIC_AXIS_DEFINITIONS.md
git add table/00_manifest/selected_death_modes_manifest.csv
git add table/00_manifest/source_trace_README.md
git status --short
git commit -m "Release regulated cell death atlas analysis code"
```

Do not run `git add .` until the large-data exclusion policy has been reviewed.

## 4A. Create and Push With GitHub CLI

Replace `<OWNER>` and `<REPO>` with your GitHub account/organization and repository name.

```bash
gh repo create <OWNER>/<REPO> --public --source=. --remote=origin --push
```

For a private repository:

```bash
gh repo create <OWNER>/<REPO> --private --source=. --remote=origin --push
```

## 4B. Push to an Existing Empty GitHub Repository

```bash
git remote add origin https://github.com/<OWNER>/<REPO>.git
git push -u origin main
```

## 5. If You Must Version Large Files

Ordinary GitHub rejects files larger than 100 MB. Prefer a DOI archive. If large files must be versioned in the repository, use Git LFS before adding them:

```bash
git lfs install
git lfs track "*.csv"
git lfs track "*.rds"
git lfs track "*.tif"
git lfs track "*.tiff"
git lfs track "*.ai"
git lfs track "*.zip"
git add .gitattributes
git commit -m "Track large data files with Git LFS"
```

Then add the selected large files intentionally, not with `git add .`.
