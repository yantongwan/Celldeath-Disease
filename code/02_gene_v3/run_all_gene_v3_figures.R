#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("CELLDEATH_ATLAS_ROOT", unset = getwd())

cmd_args <- commandArgs(FALSE)
file_idx <- grep("^--file=", cmd_args)
script_file <- if (length(file_idx) > 0) sub("^--file=", "", cmd_args[file_idx[1]]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file, mustWork = FALSE)) else ""
if (!nzchar(script_dir) || !dir.exists(script_dir)) {
  script_dir <- file.path(project_root, "reviewer_revision", "16_gene_protein_death_disease_layer_v3", "13_R_scripts")
}

source(file.path(script_dir, "00_setup_gene_v3.R"))
source(file.path(script_dir, "01_gene_summary_tables_v3.R"))
source(file.path(script_dir, "02_gene_figures_v3.R"))

paths <- gene_v3_paths(project_root)
ensure_gene_v3_dirs(paths)
summarise_gene_v3_tables(paths)
result <- generate_gene_v3_figures(paths)

cat("# Gene v3 figure run summary\n")
cat("generated_figures:", paste(result$generated, collapse = ","), "\n")
cat("skipped_figures:", paste(result$skipped, collapse = ","), "\n")
