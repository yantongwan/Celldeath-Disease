#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("CELLDEATH_ATLAS_ROOT", unset = getwd())
module_root <- file.path(project_root, "reviewer_revision", "16_gene_protein_death_disease_layer_v3")
out_root <- file.path(module_root, "17_vosviewer_enrichment_extensions", "00_v2style_figures")

for (d in c("pdf", "png", "tiff", "source_data", "legends", "15_QC_reports")) {
  dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

save_pub <- function(plot, name, width = 12, height = 8.5) {
  ggsave(file.path(out_root, "pdf", paste0(name, ".pdf")), plot, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(file.path(out_root, "png", paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 400)
  ggsave(file.path(out_root, "tiff", paste0(name, ".tiff")), plot, width = width, height = height, units = "in", dpi = 400, compression = "lzw")
}

theme_v2style <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = "#1F2933"),
      axis.ticks = element_line(linewidth = 0.22, colour = "#1F2933"),
      plot.title = element_text(face = "bold", size = base_size + 1.4),
      plot.subtitle = element_text(size = base_size, colour = "#4B5563"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#EEF2F6", colour = NA),
      strip.text = element_text(face = "bold")
    )
}

shorten <- function(x, n = 60) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

add_score <- function(df) {
  df %>%
    mutate(
      pmid_count = as.numeric(pmid_count),
      death_mode_breadth = as.numeric(death_mode_breadth),
      disease_term_breadth = as.numeric(disease_term_breadth),
      disease_system_breadth = as.numeric(disease_system_breadth),
      multi_death_signal_score = death_mode_breadth * log10(pmid_count + 1)
    )
}

palette_breadth <- function(values) {
  levels <- sort(unique(as.numeric(values)))
  setNames(grDevices::colorRampPalette(c("#D8ECF4", "#5DA5C8", "#173B57"))(length(levels)), as.character(levels))
}

module_palette <- c(
  "#173B57", "#5DA5C8", "#D8ECF4", "#F6AE2D", "#F26419", "#BC4749",
  "#2A9D8F", "#7A306C", "#386641", "#8D99AE", "#3D405B", "#E07A5F"
)

death_order_drug_figure <- c(
  "Ferroptosis",
  "Pyroptosis",
  "NETosis",
  "Necroptosis",
  "Immunogenic cell death",
  "Cuproptosis",
  "PANoptosis",
  "Disulfidptosis"
)

death_mode_palette_figure2 <- c(
  "Ferroptosis" = "#B9553C",
  "Pyroptosis" = "#D28A2E",
  "NETosis" = "#5C8F7B",
  "Necroptosis" = "#6B6BAE",
  "Immunogenic cell death" = "#A05B8F",
  "Cuproptosis" = "#3E7CB1",
  "PANoptosis" = "#1B9E77",
  "Disulfidptosis" = "#9C6A3D"
)

broad <- add_score(read_required(file.path(module_root, "07_gene_death_disease_comention", "broad_multi_death_gene_signals_primary.csv")))
death <- read_required(file.path(module_root, "07_gene_death_disease_comention", "gene_death_summary_primary.csv")) %>%
  mutate(pmid_count = as.numeric(pmid_count))
triads_raw <- read_required(file.path(module_root, "07_gene_death_disease_comention", "gene_death_disease_triads_primary.csv"))
onc <- add_score(read_required(file.path(module_root, "07_gene_death_disease_comention", "oncology_nononcology_gene_summary_primary.csv"))) %>%
  mutate(
    oncology_fraction = as.numeric(oncology_fraction),
    non_oncology_fraction = as.numeric(non_oncology_fraction)
  )
primary <- read_required(file.path(module_root, "06_gene_comention_matching", "gene_comention_matches_primary.csv"))
match_summary <- read_required(file.path(module_root, "06_gene_comention_matching", "gene_comention_match_summary.csv"))
vocab <- read_required(file.path(module_root, "03_matching_vocabulary", "shared_gene_matching_terms.csv"))
rejected_terms <- read_required(file.path(module_root, "03_matching_vocabulary", "rejected_matching_terms.csv"))

broad_ranked <- broad %>%
  arrange(desc(multi_death_signal_score), desc(death_mode_breadth), desc(pmid_count), shared_display_symbol)

broad_gene_score_lookup <- broad %>%
  transmute(
    shared_gene_id,
    gene_total_pmid_count = pmid_count,
    gene_death_mode_breadth = death_mode_breadth,
    gene_multi_death_signal_score = multi_death_signal_score
  )

triads <- triads_raw %>%
  mutate(
    triad_row_death_mode_breadth = as.numeric(death_mode_breadth),
    triad_pmid_count = as.numeric(pmid_count)
  ) %>%
  left_join(broad_gene_score_lookup, by = "shared_gene_id") %>%
  mutate(
    pmid_count = triad_pmid_count,
    death_mode_breadth = ifelse(is.na(gene_death_mode_breadth), triad_row_death_mode_breadth, gene_death_mode_breadth),
    disease_term_breadth = as.numeric(disease_term_breadth),
    disease_system_breadth = as.numeric(disease_system_breadth),
    multi_death_signal_score = death_mode_breadth * log10(pmid_count + 1)
  )

# -------------------------------------------------------------------------
# S20 v2style: workflow/QC plus score-ranked audit panel.
# -------------------------------------------------------------------------
s20_vocab <- data.frame(
  category = c("Retained matching terms", "Rejected matching terms"),
  n = c(nrow(vocab), nrow(rejected_terms))
)
s20_match_layer <- match_summary %>%
  filter(metric %in% c("primary_matches", "sensitivity_matches", "rejected_matches")) %>%
  mutate(
    value = as.numeric(value),
    metric = factor(metric, levels = c("primary_matches", "sensitivity_matches", "rejected_matches"))
  )
s20_field <- primary %>%
  count(match_field, sort = TRUE) %>%
  mutate(match_field = factor(match_field, levels = c("title", "abstract", "keywords")))
s20_top_score <- broad_ranked %>%
  slice_head(n = 20) %>%
  mutate(
    gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)),
    death_mode_breadth_factor = factor(death_mode_breadth, levels = sort(unique(death_mode_breadth)))
  )
s20_mode_palette <- palette_breadth(s20_top_score$death_mode_breadth)

write.csv(s20_vocab, file.path(out_root, "source_data", "S20_v2style_vocabulary_retention.csv"), row.names = FALSE)
write.csv(s20_match_layer, file.path(out_root, "source_data", "S20_v2style_match_layers.csv"), row.names = FALSE)
write.csv(s20_field, file.path(out_root, "source_data", "S20_v2style_match_fields.csv"), row.names = FALSE)
write.csv(s20_top_score, file.path(out_root, "source_data", "S20_v2style_top_score_audit.csv"), row.names = FALSE)

p20a <- ggplot(s20_vocab, aes(reorder(category, n), n, fill = category)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = module_palette[c(2, 6)]) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(title = "A Vocabulary audit", x = NULL, y = "Terms") +
  theme_v2style()

p20b <- ggplot(s20_match_layer, aes(metric, value, fill = metric)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = module_palette[c(1, 4, 6)]) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(title = "B Match-layer audit", x = NULL, y = "Co-mention rows") +
  theme_v2style()

p20c <- ggplot(s20_field, aes(match_field, n, fill = match_field)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = module_palette[c(2, 5, 8)]) +
  scale_y_continuous(labels = comma) +
  labs(title = "C Title/abstract/keyword contribution", x = NULL, y = "Primary matches") +
  theme_v2style()

p20d <- ggplot(s20_top_score, aes(gene_label, pmid_count, fill = death_mode_breadth_factor)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = sprintf("%.1f", multi_death_signal_score)), hjust = -0.05, size = 2.1, colour = "#263238") +
  coord_flip() +
  scale_fill_manual(values = s20_mode_palette, labels = paste0(names(s20_mode_palette), " modes")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  labs(title = "D Top score audit", subtitle = "Ranked by death-mode breadth x log10(PMID support)", x = NULL, y = "PMIDs", fill = "Breadth") +
  theme_v2style()

s20 <- (p20a | p20b) / (p20c | p20d) +
  plot_annotation(title = "Supplementary Figure S20 v2style. Shared gene co-mention workflow and score-ranked audit")
save_pub(s20, "Supplementary_Figure_S20_v3_v2style_shared_gene_comention_workflow", width = 10.8, height = 8.6)

# -------------------------------------------------------------------------
# S21 v2style: broad/selective gene landscape ranked by score.
# -------------------------------------------------------------------------
s21_top40 <- broad_ranked %>%
  slice_head(n = 40) %>%
  mutate(
    gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)),
    death_mode_breadth_factor = factor(death_mode_breadth, levels = sort(unique(death_mode_breadth)))
  )
s21_mode_palette <- palette_breadth(s21_top40$death_mode_breadth)

s21_heat_genes <- s21_top40 %>% slice_head(n = 30) %>% pull(shared_display_symbol)
s21_heat <- death %>%
  filter(shared_display_symbol %in% s21_heat_genes) %>%
  mutate(
    shared_display_symbol = factor(shared_display_symbol, levels = rev(s21_heat_genes)),
    death_mode = factor(death_mode, levels = death_order_drug_figure)
  )

s21_selective <- broad %>%
  filter(death_mode_breadth == 1) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), shared_display_symbol) %>%
  slice_head(n = 25) %>%
  mutate(
    gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)),
    death_mode_selective = factor(death_modes, levels = death_order_drug_figure),
    fill_colour = death_mode_palette_figure2[as.character(death_mode_selective)]
  )

s21_system <- broad_ranked %>%
  slice_head(n = 40) %>%
  arrange(desc(multi_death_signal_score)) %>%
  slice_head(n = 25) %>%
  mutate(gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)))

write.csv(s21_top40, file.path(out_root, "source_data", "S21_v2style_top40_broad_multi_death_signals.csv"), row.names = FALSE)
write.csv(s21_heat, file.path(out_root, "source_data", "S21_v2style_gene_death_heatmap.csv"), row.names = FALSE)
write.csv(s21_selective, file.path(out_root, "source_data", "S21_v2style_death_mode_selective_top25.csv"), row.names = FALSE)
write.csv(s21_system, file.path(out_root, "source_data", "S21_v2style_disease_system_breadth_top25.csv"), row.names = FALSE)
write.csv(data.frame(
  death_mode = death_order_drug_figure,
  left_to_right_order = seq_along(death_order_drug_figure),
  reference_figure = "Figure_same_drug_multi_death_induction_sensitization_unfiltered_reference_mapped"
), file.path(out_root, "source_data", "S21_v2style_panel_C_death_mode_order.csv"), row.names = FALSE)

p21a <- ggplot(s21_top40, aes(gene_label, pmid_count, fill = death_mode_breadth_factor)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = paste0(death_mode_breadth, " modes")), hjust = -0.05, size = 1.9, colour = "#263238") +
  coord_flip() +
  scale_fill_manual(values = s21_mode_palette, labels = paste0(names(s21_mode_palette), " modes")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.20))) +
  labs(title = "A Top 40 broad multi-death signals", subtitle = "Ranked by breadth x log10(PMID support)", x = NULL, y = "PMIDs", fill = "Breadth") +
  theme_v2style(base_size = 7)

p21b <- ggplot(s21_selective, aes(gene_label, pmid_count, fill = death_mode_selective)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_manual(values = death_mode_palette_figure2, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B Top death-mode-selective signals", subtitle = "Breadth = 1; ranked by the same score. Color follows Figure 2A death-mode palette.", x = NULL, y = "PMIDs", fill = "Death mode") +
  theme_v2style(base_size = 7)

p21c <- ggplot(s21_heat, aes(death_mode, shared_display_symbol, fill = log10(pmid_count + 1))) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient(low = "#F7FBFF", high = "#7A306C") +
  labs(title = "C Gene x death-mode heatmap", subtitle = "Top 30 genes from score ranking", x = NULL, y = NULL, fill = "log10(PMID+1)") +
  theme_v2style(base_size = 7) +
  theme(axis.text.x = element_text(angle = 38, hjust = 1))

p21d <- ggplot(s21_system, aes(gene_label, disease_system_breadth, fill = multi_death_signal_score)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_gradient(low = "#C7F9CC", high = "#386641") +
  labs(title = "D Disease-system breadth among score-ranked genes", x = NULL, y = "Disease systems", fill = "Score") +
  theme_v2style(base_size = 7)

s21 <- (p21a | p21b) / (p21c | p21d) +
  plot_annotation(title = "Supplementary Figure S21 v2style. Broad and selective shared gene co-mention signals")
save_pub(s21, "Supplementary_Figure_S21_v3_v2style_broad_and_selective_gene_comention", width = 11.2, height = 10)

# -------------------------------------------------------------------------
# S22 v2style: gene-death-disease triads ranked by score.
# -------------------------------------------------------------------------
s22_triads_ranked <- triads %>%
  mutate(
    title_match_count = as.numeric(title_match_count),
    abstract_match_count = as.numeric(abstract_match_count),
    keyword_match_count = as.numeric(keyword_match_count),
    ambiguity_flag_fraction = as.numeric(ambiguity_flag_fraction),
    low_count_warning_fraction = as.numeric(low_count_warning_fraction),
    label = shorten(paste(shared_display_symbol, death_mode, disease_term, sep = " | "), 72)
  ) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), desc(title_match_count + abstract_match_count), shared_display_symbol)

s22_top <- s22_triads_ranked %>%
  filter(evidence_grade %in% c("strong_comention_signal", "moderate_comention_signal", "exploratory_comention_signal")) %>%
  slice_head(n = 30) %>%
  mutate(label = factor(label, levels = rev(label)))

s22_system <- s22_triads_ranked %>%
  slice_head(n = 1000) %>%
  group_by(disease_system) %>%
  summarise(
    triads = n(),
    median_score = median(multi_death_signal_score, na.rm = TRUE),
    top_score = max(multi_death_signal_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(top_score), desc(triads)) %>%
  slice_head(n = 20) %>%
  mutate(disease_system = factor(disease_system, levels = rev(disease_system)))

s22_onc <- s22_triads_ranked %>%
  slice_head(n = 1500) %>%
  mutate(layer = ifelse(as.numeric(oncology_fraction) >= 0.5, "oncology-enriched", "non-oncology/mixed")) %>%
  group_by(layer, evidence_grade) %>%
  summarise(triads = n(), median_score = median(multi_death_signal_score, na.rm = TRUE), .groups = "drop")

s22_caution <- s22_triads_ranked %>%
  slice_head(n = min(5000, nrow(s22_triads_ranked)))

write.csv(s22_top, file.path(out_root, "source_data", "S22_v2style_top_score_triads.csv"), row.names = FALSE)
write.csv(s22_system, file.path(out_root, "source_data", "S22_v2style_top_score_disease_system_distribution.csv"), row.names = FALSE)
write.csv(s22_onc, file.path(out_root, "source_data", "S22_v2style_top_score_oncology_distribution.csv"), row.names = FALSE)
write.csv(s22_caution, file.path(out_root, "source_data", "S22_v2style_top_score_caution_scatter.csv"), row.names = FALSE)

p22a <- ggplot(s22_top, aes(label, pmid_count, fill = evidence_grade)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_manual(values = c(strong_comention_signal = "#173B57", moderate_comention_signal = "#F6AE2D", exploratory_comention_signal = "#5DA5C8", weak_or_ambiguous = "#BC4749")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A Top score-ranked gene-death-disease triads", subtitle = "Ranked by breadth x log10(triad PMID support)", x = NULL, y = "Triad PMIDs", fill = "Grade") +
  theme_v2style(base_size = 6.7)

p22b <- ggplot(s22_system, aes(disease_system, triads, fill = top_score)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_gradient(low = "#D8ECF4", high = "#173B57") +
  scale_y_continuous(labels = comma) +
  labs(title = "B Disease systems in top score-ranked triads", x = NULL, y = "Top-ranked triads", fill = "Top score") +
  theme_v2style(base_size = 7)

p22c <- ggplot(s22_onc, aes(layer, triads, fill = evidence_grade)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = c(strong_comention_signal = "#173B57", moderate_comention_signal = "#F6AE2D", exploratory_comention_signal = "#5DA5C8", weak_or_ambiguous = "#BC4749")) +
  scale_y_continuous(labels = comma) +
  labs(title = "C Oncology layer among top score-ranked triads", x = NULL, y = "Triads", fill = "Grade") +
  theme_v2style(base_size = 7)

p22d <- ggplot(s22_caution, aes(low_count_warning_fraction, ambiguity_flag_fraction, color = multi_death_signal_score)) +
  geom_point(alpha = 0.35, size = 1) +
  scale_color_gradient(low = "#5DA5C8", high = "#BC4749") +
  labs(title = "D Top-score triad caution map", x = "Low-count warning fraction", y = "Ambiguity fraction", color = "Score") +
  theme_v2style(base_size = 7)

s22 <- (p22a | p22b) / (p22c | p22d) +
  plot_annotation(title = "Supplementary Figure S22 v2style. Score-ranked gene-death-disease co-mention triads")
save_pub(s22, "Supplementary_Figure_S22_v3_v2style_gene_death_disease_triads", width = 11.2, height = 10)

# -------------------------------------------------------------------------
# S23 v2style: oncology/non-oncology layers ranked by score.
# -------------------------------------------------------------------------
s23_onc_top <- onc %>%
  filter(oncology_fraction >= 0.75) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), shared_display_symbol) %>%
  slice_head(n = 30) %>%
  mutate(gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)))

s23_non_top <- onc %>%
  filter(non_oncology_fraction >= 0.75) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), shared_display_symbol) %>%
  slice_head(n = 30) %>%
  mutate(gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)))

s23_mixed_top <- onc %>%
  filter(oncology_fraction > 0.25, oncology_fraction < 0.75) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), shared_display_symbol) %>%
  slice_head(n = 25) %>%
  mutate(gene_label = factor(shared_display_symbol, levels = rev(shared_display_symbol)))

s23_class <- onc %>%
  group_by(oncology_enrichment_class) %>%
  summarise(genes = n(), median_score = median(multi_death_signal_score, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_score))

write.csv(s23_onc_top, file.path(out_root, "source_data", "S23_v2style_oncology_enriched_top_score_genes.csv"), row.names = FALSE)
write.csv(s23_non_top, file.path(out_root, "source_data", "S23_v2style_non_oncology_enriched_top_score_genes.csv"), row.names = FALSE)
write.csv(s23_mixed_top, file.path(out_root, "source_data", "S23_v2style_mixed_top_score_genes.csv"), row.names = FALSE)
write.csv(s23_class, file.path(out_root, "source_data", "S23_v2style_enrichment_class_score_summary.csv"), row.names = FALSE)

p23a <- ggplot(s23_onc_top, aes(gene_label, pmid_count, fill = multi_death_signal_score)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_gradient(low = "#FDE68A", high = "#BC4749") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A Oncology-enriched genes", subtitle = "Ranked by multi-death signal score", x = NULL, y = "PMIDs", fill = "Score") +
  theme_v2style(base_size = 7)

p23b <- ggplot(s23_non_top, aes(gene_label, pmid_count, fill = multi_death_signal_score)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_gradient(low = "#D8ECF4", high = "#173B57") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B Non-oncology-enriched genes", subtitle = "Ranked by multi-death signal score", x = NULL, y = "PMIDs", fill = "Score") +
  theme_v2style(base_size = 7)

p23c <- ggplot(s23_mixed_top, aes(gene_label, pmid_count, fill = multi_death_signal_score)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_gradient(low = "#C7F9CC", high = "#386641") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "C Mixed oncology/non-oncology genes", subtitle = "Ranked by multi-death signal score", x = NULL, y = "PMIDs", fill = "Score") +
  theme_v2style(base_size = 7)

p23d <- ggplot(onc, aes(oncology_fraction, multi_death_signal_score, color = oncology_enrichment_class)) +
  geom_point(alpha = 0.55, size = 1.5) +
  scale_color_manual(values = c(oncology_enriched = "#BC4749", non_oncology_enriched = "#173B57", mixed = "#2A9D8F")) +
  labs(title = "D Oncology fraction versus multi-death score", x = "Oncology fraction", y = "Multi-death signal score", color = "Class") +
  theme_v2style(base_size = 7)

s23 <- (p23a | p23b) / (p23c | p23d) +
  plot_annotation(title = "Supplementary Figure S23 v2style. Score-ranked oncology and non-oncology shared gene layers")
save_pub(s23, "Supplementary_Figure_S23_v3_v2style_oncology_nononcology_gene_layer", width = 11.2, height = 10)

writeLines(c(
  "# v3 v2style S20-S23 figure legend",
  "",
  "All ranked gene and triad panels use multi_death_signal_score = death_mode_breadth x log10(PMID count + 1).",
  "S20 retains workflow/QC panels and adds a score-ranked top co-mention audit panel.",
  "S21 ranks broad, selective, heatmap, and disease-system-breadth displays using the same score where gene ranking is required. Panel C uses the same left-to-right death-mode order as Figure_same_drug_multi_death_induction_sensitization_unfiltered_reference_mapped.",
  "S22 ranks gene-death-disease triads by the same score, using the triad PMID count and gene death-mode breadth.",
  "S23 ranks oncology-enriched, non-oncology-enriched, and mixed genes by the same score.",
  "",
  "Interpretation: all panels show literature-level shared gene co-mention signals only. They do not establish causal gene function or pathway regulation."
), file.path(out_root, "legends", "Supplementary_Figures_S20_S23_v2style_legend.md"))

manifest <- data.frame(
  figure = c("S20_v2style", "S21_v2style", "S22_v2style", "S23_v2style"),
  pdf = file.path(out_root, "pdf", c(
    "Supplementary_Figure_S20_v3_v2style_shared_gene_comention_workflow.pdf",
    "Supplementary_Figure_S21_v3_v2style_broad_and_selective_gene_comention.pdf",
    "Supplementary_Figure_S22_v3_v2style_gene_death_disease_triads.pdf",
    "Supplementary_Figure_S23_v3_v2style_oncology_nononcology_gene_layer.pdf"
  )),
  sorting_rule = "multi_death_signal_score = death_mode_breadth * log10(PMID count + 1)",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_root, "15_QC_reports", "S20_S23_v2style_generated_manifest.csv"), row.names = FALSE)

cat("Generated v2style S20-S23 figures in: ", out_root, "\n", sep = "")
cat("Ranking rule: multi_death_signal_score = death_mode_breadth * log10(PMID count + 1)\n")
