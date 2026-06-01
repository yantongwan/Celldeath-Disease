#!/usr/bin/env Rscript

summarise_gene_v3_tables <- function(paths) {
  match_summary <- read_gene_v3_csv(file.path(paths$root, "06_gene_comention_matching", "gene_comention_match_summary.csv"))
  vocab_summary <- read_gene_v3_csv(file.path(paths$root, "03_matching_vocabulary", "matching_vocabulary_summary.csv"))
  universe_summary <- read_gene_v3_csv(file.path(paths$root, "02_shared_gene_universe", "shared_gene_universe_summary.csv"))
  gate <- read_gene_v3_csv(file.path(paths$qc, "v3_gene_figure_gate_report.csv"))

  summary_rows <- rbind(
    data.frame(section = "universe", metric = universe_summary$metric, value = universe_summary$value),
    data.frame(section = "vocabulary", metric = vocab_summary$metric, value = vocab_summary$value),
    data.frame(section = "matching", metric = match_summary$metric, value = match_summary$value)
  )
  write.csv(summary_rows, file.path(paths$source_data, "S20_summary_metrics_source_data.csv"), row.names = FALSE)
  write.csv(gate, file.path(paths$source_data, "S20_S23_figure_gate_source_data.csv"), row.names = FALSE)
  invisible(summary_rows)
}

