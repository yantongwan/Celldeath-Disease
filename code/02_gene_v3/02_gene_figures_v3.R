#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

gate_passed <- function(gate, figure_id) {
  rows <- gate[gate$figure == figure_id, , drop = FALSE]
  nrow(rows) > 0 && rows$gate_result[1] == "PASS"
}

save_gene_v3_plot <- function(plot, paths, name, width = 13, height = 9) {
  pdf_path <- file.path(paths$fig_pdf, paste0(name, ".pdf"))
  png_path <- file.path(paths$fig_png, paste0(name, ".png"))
  tiff_path <- file.path(paths$fig_tiff, paste0(name, ".tiff"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 320)
  ggsave(tiff_path, plot, width = width, height = height, units = "in", dpi = 320, compression = "lzw")
  c(pdf = pdf_path, png = png_path, tiff = tiff_path)
}

theme_gene_v3 <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "#4A5568"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

plot_s20 <- function(paths, gate) {
  if (!gate_passed(gate, "S20")) return(NULL)
  match_summary <- read_gene_v3_csv(file.path(paths$root, "06_gene_comention_matching", "gene_comention_match_summary.csv"))
  vocab <- read_gene_v3_csv(file.path(paths$root, "03_matching_vocabulary", "shared_gene_matching_terms.csv"))
  rejected_terms <- read_gene_v3_csv(file.path(paths$root, "03_matching_vocabulary", "rejected_matching_terms.csv"))
  primary <- read_gene_v3_csv(file.path(paths$root, "06_gene_comention_matching", "gene_comention_matches_primary.csv"))

  vocab_panel <- data.frame(
    category = c("Retained matching terms", "Rejected matching terms"),
    n = c(nrow(vocab), nrow(rejected_terms))
  )
  ambiguity_panel <- vocab %>%
    count(ambiguity_status, sort = TRUE) %>%
    mutate(ambiguity_status = reorder(ambiguity_status, n))
  match_panel <- match_summary %>%
    filter(metric %in% c("primary_matches", "sensitivity_matches", "rejected_matches")) %>%
    mutate(metric = factor(metric, levels = c("primary_matches", "sensitivity_matches", "rejected_matches")))
  field_panel <- primary %>%
    count(match_field, sort = TRUE) %>%
    mutate(match_field = factor(match_field, levels = c("title", "abstract", "keywords")))

  write.csv(vocab_panel, file.path(paths$source_data, "S20_vocabulary_retention_source_data.csv"), row.names = FALSE)
  write.csv(ambiguity_panel, file.path(paths$source_data, "S20_ambiguity_status_source_data.csv"), row.names = FALSE)
  write.csv(match_panel, file.path(paths$source_data, "S20_match_layer_source_data.csv"), row.names = FALSE)
  write.csv(field_panel, file.path(paths$source_data, "S20_match_field_source_data.csv"), row.names = FALSE)

  p1 <- ggplot(vocab_panel, aes(category, n, fill = category)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    scale_fill_manual(values = gene_v3_palette[c(2, 6)]) +
    coord_flip() +
    labs(title = "A. Matching vocabulary audit", x = NULL, y = "Terms") +
    theme_gene_v3()

  p2 <- ggplot(ambiguity_panel, aes(ambiguity_status, n, fill = ambiguity_status)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = rep(gene_v3_palette, length.out = nrow(ambiguity_panel))) +
    coord_flip() +
    labs(title = "B. Retained term ambiguity status", x = NULL, y = "Terms") +
    theme_gene_v3()

  p3 <- ggplot(match_panel, aes(metric, as.numeric(value), fill = metric)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    scale_fill_manual(values = gene_v3_palette[c(3, 5, 8)]) +
    coord_flip() +
    labs(title = "C. Match layers", x = NULL, y = "Co-mention rows") +
    theme_gene_v3()

  p4 <- ggplot(field_panel, aes(match_field, n, fill = match_field)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    scale_fill_manual(values = gene_v3_palette[c(1, 5, 7)]) +
    labs(title = "D. Title/abstract/keyword contribution", x = NULL, y = "Primary matches") +
    theme_gene_v3()

  plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Supplementary Figure S20. Shared human-mouse gene co-mention layer")
  save_gene_v3_plot(plot, paths, "Supplementary_Figure_S20_shared_gene_comention_layer")
  writeLines(c(
    "# Supplementary Figure S20 Legend",
    "",
    "Shared human-mouse gene co-mention workflow and QC. Panels show retained/rejected matching terms, ambiguity status, primary/sensitivity/rejected match layers, and title/abstract/keyword contributions. These are literature-level co-mention signals only."
  ), file.path(paths$legends, "Supplementary_Figure_S20_legend.md"))
  "S20"
}

plot_s21 <- function(paths, gate) {
  if (!gate_passed(gate, "S21")) return(NULL)
  broad <- read_gene_v3_csv(file.path(paths$root, "07_gene_death_disease_comention", "broad_multi_death_gene_signals_primary.csv"))
  death <- read_gene_v3_csv(file.path(paths$root, "07_gene_death_disease_comention", "gene_death_summary_primary.csv"))
  selective <- read_gene_v3_csv(file.path(paths$root, "07_gene_death_disease_comention", "death_mode_selective_gene_signals_primary.csv"))

  top_broad <- broad %>% arrange(desc(death_mode_breadth), desc(pmid_count)) %>% slice_head(n = 20)
  top_genes <- top_broad$shared_display_symbol[1:min(25, nrow(top_broad))]
  heat <- death %>%
    filter(shared_display_symbol %in% top_genes) %>%
    mutate(shared_display_symbol = factor(shared_display_symbol, levels = rev(top_genes)))
  top_selective <- selective %>% arrange(desc(pmid_count)) %>% slice_head(n = 20)
  system_breadth <- broad %>% arrange(desc(disease_system_breadth), desc(pmid_count)) %>% slice_head(n = 20)

  write.csv(top_broad, file.path(paths$source_data, "S21_top_broad_genes_source_data.csv"), row.names = FALSE)
  write.csv(heat, file.path(paths$source_data, "S21_gene_death_heatmap_source_data.csv"), row.names = FALSE)
  write.csv(top_selective, file.path(paths$source_data, "S21_death_mode_selective_source_data.csv"), row.names = FALSE)
  write.csv(system_breadth, file.path(paths$source_data, "S21_disease_system_breadth_source_data.csv"), row.names = FALSE)

  p1 <- ggplot(top_broad, aes(reorder(shared_display_symbol, pmid_count), pmid_count, fill = death_mode_breadth)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#BEE9E8", high = "#1B4965") +
    coord_flip() +
    labs(title = "A. Broad multi-death shared genes", x = NULL, y = "PMID count", fill = "Death modes") +
    theme_gene_v3()

  p2 <- ggplot(top_selective, aes(reorder(shared_display_symbol, pmid_count), pmid_count, fill = disease_system_breadth)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#FDE68A", high = "#F26419") +
    coord_flip() +
    labs(title = "B. Death-mode selective shared genes", x = NULL, y = "PMID count", fill = "Disease systems") +
    theme_gene_v3()

  p3 <- ggplot(heat, aes(death_mode, shared_display_symbol, fill = log10(pmid_count + 1))) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "#F7FBFF", high = "#7A306C") +
    labs(title = "C. Gene by death-mode co-mention heatmap", x = NULL, y = NULL, fill = "log10(PMID+1)") +
    theme_gene_v3(9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  p4 <- ggplot(system_breadth, aes(reorder(shared_display_symbol, disease_system_breadth), disease_system_breadth, fill = pmid_count)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#C7F9CC", high = "#386641") +
    coord_flip() +
    labs(title = "D. Disease-system breadth", x = NULL, y = "Disease systems", fill = "PMIDs") +
    theme_gene_v3()

  plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Supplementary Figure S21. Broad and death-mode-selective shared gene co-mentions")
  save_gene_v3_plot(plot, paths, "Supplementary_Figure_S21_broad_and_selective_gene_comention")
  writeLines(c(
    "# Supplementary Figure S21 Legend",
    "",
    "Broad and selective shared gene co-mention signals. Ranking uses primary title/abstract/keyword matches aggregated by shared human-mouse gene ID. Breadth is a literature co-mention measure, not mechanistic evidence."
  ), file.path(paths$legends, "Supplementary_Figure_S21_legend.md"))
  "S21"
}

plot_s22 <- function(paths, gate) {
  if (!gate_passed(gate, "S22")) return(NULL)
  triads <- read_gene_v3_csv(file.path(paths$root, "07_gene_death_disease_comention", "gene_death_disease_triads_primary.csv"))
  strong <- triads %>%
    filter(evidence_grade %in% c("strong_comention_signal", "moderate_comention_signal")) %>%
    arrange(desc(pmid_count)) %>%
    slice_head(n = 25) %>%
    mutate(label = paste(shared_display_symbol, death_mode, disease_term, sep = " | "))
  system_dist <- triads %>%
    filter(evidence_grade %in% c("strong_comention_signal", "moderate_comention_signal")) %>%
    count(disease_system, sort = TRUE) %>%
    slice_head(n = 20)
  onc_dist <- triads %>%
    count(oncology_fraction = ifelse(as.numeric(oncology_fraction) >= 0.5, "oncology-enriched", "non-oncology/mixed"), evidence_grade)
  caution <- triads %>%
    mutate(ambiguity_flag_fraction = as.numeric(ambiguity_flag_fraction),
           low_count_warning_fraction = as.numeric(low_count_warning_fraction)) %>%
    sample_n(min(5000, nrow(.)))

  write.csv(strong, file.path(paths$source_data, "S22_top_triads_source_data.csv"), row.names = FALSE)
  write.csv(system_dist, file.path(paths$source_data, "S22_disease_system_distribution_source_data.csv"), row.names = FALSE)
  write.csv(onc_dist, file.path(paths$source_data, "S22_oncology_distribution_source_data.csv"), row.names = FALSE)
  write.csv(caution, file.path(paths$source_data, "S22_caution_scatter_source_data.csv"), row.names = FALSE)

  p1 <- ggplot(strong, aes(reorder(label, pmid_count), pmid_count, fill = evidence_grade)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c(strong_comention_signal = "#1B4965", moderate_comention_signal = "#F6AE2D")) +
    coord_flip() +
    labs(title = "A. Top gene-death-disease triads", x = NULL, y = "PMID count", fill = "Grade") +
    theme_gene_v3(8)

  p2 <- ggplot(system_dist, aes(reorder(disease_system, n), n, fill = disease_system)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = rep(gene_v3_palette, length.out = nrow(system_dist))) +
    coord_flip() +
    labs(title = "B. Disease-system distribution", x = NULL, y = "Strong/moderate triads") +
    theme_gene_v3()

  p3 <- ggplot(onc_dist, aes(oncology_fraction, n, fill = evidence_grade)) +
    geom_col(position = "stack", width = 0.72) +
    scale_fill_manual(values = c(strong_comention_signal = "#1B4965", moderate_comention_signal = "#F6AE2D", exploratory_comention_signal = "#BEE9E8", weak_or_ambiguous = "#BC4749")) +
    labs(title = "C. Oncology versus non-oncology triads", x = NULL, y = "Triads", fill = "Grade") +
    theme_gene_v3()

  p4 <- ggplot(caution, aes(low_count_warning_fraction, ambiguity_flag_fraction, color = evidence_grade)) +
    geom_point(alpha = 0.35, size = 1) +
    scale_color_manual(values = c(strong_comention_signal = "#1B4965", moderate_comention_signal = "#F6AE2D", exploratory_comention_signal = "#2A9D8F", weak_or_ambiguous = "#BC4749")) +
    labs(title = "D. Low-count and ambiguity caution", x = "Low-count warning fraction", y = "Ambiguity fraction", color = "Grade") +
    theme_gene_v3()

  plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Supplementary Figure S22. Gene-death-disease literature co-mention triads")
  save_gene_v3_plot(plot, paths, "Supplementary_Figure_S22_gene_death_disease_triads")
  writeLines(c(
    "# Supplementary Figure S22 Legend",
    "",
    "Gene-death-disease co-mention triads from primary title/abstract/keyword matches. Evidence grades reflect PMID support, title/abstract support, ambiguity fraction, and low-count warnings. They do not indicate causal gene regulation."
  ), file.path(paths$legends, "Supplementary_Figure_S22_legend.md"))
  "S22"
}

plot_s23 <- function(paths, gate) {
  if (!gate_passed(gate, "S23")) return(NULL)
  onc <- read_gene_v3_csv(file.path(paths$root, "07_gene_death_disease_comention", "oncology_nononcology_gene_summary_primary.csv")) %>%
    mutate(oncology_fraction = as.numeric(oncology_fraction),
           non_oncology_fraction = as.numeric(non_oncology_fraction),
           pmid_count = as.numeric(pmid_count))
  oncology_top <- onc %>% filter(oncology_fraction >= 0.75) %>% arrange(desc(pmid_count)) %>% slice_head(n = 20)
  non_top <- onc %>% filter(non_oncology_fraction >= 0.75) %>% arrange(desc(pmid_count)) %>% slice_head(n = 20)
  class_counts <- onc %>% count(oncology_enrichment_class)
  both <- onc %>% mutate(layer = case_when(
    oncology_fraction > 0 & non_oncology_fraction > 0 ~ "both",
    oncology_fraction > 0 ~ "oncology_only",
    TRUE ~ "non_oncology_only"
  )) %>% count(layer)

  write.csv(oncology_top, file.path(paths$source_data, "S23_oncology_enriched_source_data.csv"), row.names = FALSE)
  write.csv(non_top, file.path(paths$source_data, "S23_non_oncology_enriched_source_data.csv"), row.names = FALSE)
  write.csv(class_counts, file.path(paths$source_data, "S23_enrichment_class_counts_source_data.csv"), row.names = FALSE)
  write.csv(both, file.path(paths$source_data, "S23_overlap_layer_counts_source_data.csv"), row.names = FALSE)

  p1 <- ggplot(oncology_top, aes(reorder(shared_display_symbol, pmid_count), pmid_count, fill = oncology_fraction)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#FDE68A", high = "#BC4749") +
    coord_flip() +
    labs(title = "A. Oncology-enriched shared genes", x = NULL, y = "PMID count", fill = "Oncology fraction") +
    theme_gene_v3()

  p2 <- ggplot(non_top, aes(reorder(shared_display_symbol, pmid_count), pmid_count, fill = non_oncology_fraction)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#BEE9E8", high = "#1B4965") +
    coord_flip() +
    labs(title = "B. Non-oncology-enriched shared genes", x = NULL, y = "PMID count", fill = "Non-oncology fraction") +
    theme_gene_v3()

  p3 <- ggplot(class_counts, aes(oncology_enrichment_class, n, fill = oncology_enrichment_class)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = gene_v3_palette[c(6, 3, 1)]) +
    labs(title = "C. Enrichment class distribution", x = NULL, y = "Genes") +
    theme_gene_v3()

  p4 <- ggplot(onc, aes(oncology_fraction, pmid_count, color = oncology_enrichment_class)) +
    geom_point(alpha = 0.55, size = 1.6) +
    scale_y_log10() +
    scale_color_manual(values = c(oncology_enriched = "#BC4749", non_oncology_enriched = "#1B4965", mixed = "#2A9D8F")) +
    labs(title = "D. Oncology fraction versus PMID support", x = "Oncology fraction", y = "PMID count (log10)", color = "Class") +
    theme_gene_v3()

  plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = "Supplementary Figure S23. Oncology and non-oncology shared gene co-mention layers")
  save_gene_v3_plot(plot, paths, "Supplementary_Figure_S23_oncology_nononcology_gene_layer")
  writeLines(c(
    "# Supplementary Figure S23 Legend",
    "",
    "Oncology and non-oncology shared gene co-mention layers. Fractions are calculated from primary co-mention rows across death-mode-disease contexts. These are literature-level co-occurrence strata."
  ), file.path(paths$legends, "Supplementary_Figure_S23_legend.md"))
  "S23"
}

generate_gene_v3_figures <- function(paths) {
  gate <- read_gene_v3_csv(file.path(paths$qc, "v3_gene_figure_gate_report.csv"))
  generated <- c(
    plot_s20(paths, gate),
    plot_s21(paths, gate),
    plot_s22(paths, gate),
    plot_s23(paths, gate)
  )
  generated <- generated[!vapply(generated, is.null, logical(1))]
  skipped <- gate$figure[gate$gate_result != "PASS"]
  write.csv(data.frame(generated = generated), file.path(paths$source_data, "generated_figures_v3.csv"), row.names = FALSE)
  write.csv(data.frame(skipped = skipped), file.path(paths$source_data, "skipped_figures_v3.csv"), row.names = FALSE)
  list(generated = generated, skipped = skipped)
}

