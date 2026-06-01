#!/usr/bin/env Rscript

gene_v3_paths <- function(project_root) {
  module_root <- file.path(project_root, "reviewer_revision", "16_gene_protein_death_disease_layer_v3")
  list(
    root = module_root,
    fig_pdf = file.path(module_root, "10_publication_grade_figures", "pdf"),
    fig_png = file.path(module_root, "10_publication_grade_figures", "png"),
    fig_tiff = file.path(module_root, "10_publication_grade_figures", "tiff"),
    source_data = file.path(module_root, "10_publication_grade_figures", "source_data"),
    legends = file.path(module_root, "10_publication_grade_figures", "legends"),
    qc = file.path(module_root, "15_QC_reports")
  )
}

ensure_gene_v3_dirs <- function(paths) {
  for (p in unlist(paths[c("fig_pdf", "fig_png", "fig_tiff", "source_data", "legends")])) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

read_gene_v3_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

gene_v3_palette <- c(
  "#1B4965", "#5FA8D3", "#62B6CB", "#BEE9E8", "#F6AE2D",
  "#F26419", "#7A306C", "#2A9D8F", "#6A994E", "#BC4749"
)

