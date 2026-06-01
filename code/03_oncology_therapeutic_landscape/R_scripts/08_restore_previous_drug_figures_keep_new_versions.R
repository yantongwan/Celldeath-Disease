script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
dir_create(fig_dir)

class_labels <- c(
  "chemotherapy" = "Chemotherapy",
  "targeted_therapy" = "Targeted therapy",
  "immunotherapy_checkpoint" = "Checkpoint immunotherapy",
  "radiotherapy" = "Radiotherapy",
  "ferroptosis_inducer_tool_compound" = "Ferroptosis tool compound",
  "kinase_inhibitor" = "Kinase inhibitor",
  "proteasome_inhibitor" = "Proteasome inhibitor",
  "metabolic_modulator" = "Metabolic modulator",
  "natural_product" = "Natural product",
  "experimental_tool_compound" = "Experimental tool compound",
  "unclear_or_mixed" = "Unclear/mixed"
)

copy_export_set <- function(from_base, to_base) {
  for (ext in c("pdf", "png", "tiff")) {
    src <- paste0(from_base, ".", ext)
    dst <- paste0(to_base, ".", ext)
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

short_death <- function(x) {
  recode(as.character(x), "Immunogenic cell death" = "ICD", .default = as.character(x))
}

death_mode_levels_short <- unique(short_death(names(death_mode_palette)))

make_previous_global_triads <- function() {
  top_path <- file.path(ROOT, "06_drug_death_mode_triads/top_oncology_drug_death_triads.csv")
  triad_path <- file.path(ROOT, "06_drug_death_mode_triads/oncology_drug_death_triads.csv")
  triads <- if (file.exists(top_path)) read_csv_required(top_path) else read_csv_required(triad_path)

  triads <- triads %>%
    mutate(
      induction_sensitization_count = coalesce(induction_count, 0) + coalesce(sensitization_count, 0),
      triad_label = paste(tumor_family, drug_display, short_death(death_mode), sep = " | "),
      triad_label = str_trunc(triad_label, 58),
      death_mode_plot = factor(short_death(death_mode), levels = death_mode_levels_short),
      therapy_class = if_else(is.na(therapy_class) | therapy_class == "", "unclear_or_mixed", therapy_class)
    )

  top_triads <- triads %>%
    arrange(desc(pmid_count), desc(induction_sensitization_count), desc(high_confidence_relation_count)) %>%
    slice_head(n = 30) %>%
    mutate(triad_label = factor(triad_label, levels = rev(unique(triad_label))))

  p1 <- ggplot(top_triads, aes(x = pmid_count, y = triad_label, fill = therapy_class)) +
    geom_col(width = 0.72, alpha = 0.88) +
    geom_point(aes(x = induction_sensitization_count), shape = 21, fill = "white",
               colour = "grey15", stroke = 0.25, size = 1.9) +
    scale_fill_manual(values = therapy_palette, labels = class_labels[names(therapy_palette)],
                      guide = "none", na.value = "grey70") +
    labs(
      title = "A  Highest-support tumor-drug-death triads",
      x = "PMID count; dot = induction/sensitization PMIDs",
      y = NULL
    ) +
    theme_oncology_atlas(base_size = 6.3)

  p2 <- top_triads %>%
    count(death_mode_plot, therapy_class, wt = pmid_count, name = "pmid_count") %>%
    ggplot(aes(x = death_mode_plot, y = class_labels[therapy_class], fill = log1p(pmid_count))) +
    geom_tile(colour = "white", linewidth = 0.28) +
    scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p PMID count") +
    labs(title = "B  Class composition", x = NULL, y = NULL) +
    theme_oncology_atlas(base_size = 6.5) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  p3 <- top_triads %>%
    count(death_mode_plot, wt = pmid_count, name = "pmid_count") %>%
    ggplot(aes(x = reorder(death_mode_plot, pmid_count), y = pmid_count, fill = death_mode_plot)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = death_mode_palette, na.value = "grey70") +
    labs(title = "C  Death-mode contribution", x = NULL, y = "PMID count") +
    theme_oncology_atlas(base_size = 6.5)

  plot <- (p1 | (p2 / p3)) +
    plot_layout(widths = c(1.35, 1)) +
    plot_annotation(
      title = "Oncology tumor-drug-death triads remain literature-level hypotheses",
      subtitle = "Restored global triad view. Signals are PMID-supported co-mention and sentence-rule candidates, not validated drug-induced mechanisms.",
      theme = theme(plot.title = element_text(face = "bold", size = 10),
                    plot.subtitle = element_text(size = 8, colour = "grey30"))
    )

  save_pub_r(plot, file.path(fig_dir, "Figure_oncology_drug_death_triads"),
             width_mm = 183, height_mm = 165)
  save_pub_r(plot, file.path(fig_dir, "Figure_oncology_drug_death_triads_previous_global"),
             width_mm = 183, height_mm = 165)
}

make_previous_polydeath <- function() {
  poly <- read_csv_required(file.path(ROOT, "07_same_drug_multi_death_signals/drug_polydeath_signal_table.csv"))
  triads <- read_csv_required(file.path(ROOT, "06_drug_death_mode_triads/oncology_drug_death_triads.csv")) %>%
    mutate(
      induction_sensitization_count = coalesce(induction_count, 0) + coalesce(sensitization_count, 0),
      death_mode_plot = factor(short_death(death_mode), levels = death_mode_levels_short)
    )

  display_drugs <- poly %>%
    arrange(desc(polydeath_index), desc(death_mode_breadth), desc(high_confidence_relation_count)) %>%
    slice_head(n = 20) %>%
    mutate(
      therapy_class = if_else(is.na(therapy_class) | therapy_class == "", "unclear_or_mixed", therapy_class),
      drug_display = factor(drug_display, levels = rev(drug_display))
    )

  p1 <- ggplot(display_drugs, aes(x = polydeath_index, y = drug_display, fill = therapy_class)) +
    geom_col(width = 0.7, alpha = 0.88) +
    geom_point(aes(x = death_mode_breadth), shape = 21, fill = "white",
               colour = "grey15", stroke = 0.25, size = 2.0) +
    scale_fill_manual(values = therapy_palette, labels = class_labels[names(therapy_palette)],
                      name = "Therapy class", na.value = "grey70") +
    labs(
      title = "A  Same-drug multi-death literature ranking",
      x = "Polydeath literature index; white dot = death-mode breadth",
      y = NULL
    ) +
    theme_oncology_atlas(base_size = 6.4) +
    theme(legend.position = "bottom")

  mode_matrix <- triads %>%
    filter(drug_normalized %in% display_drugs$drug_normalized) %>%
    group_by(drug_normalized, drug_display, death_mode_plot) %>%
    summarise(
      pmid_count = sum(pmid_count, na.rm = TRUE),
      high_confidence_relation_count = sum(high_confidence_relation_count, na.rm = TRUE),
      induction_sensitization_count = sum(induction_sensitization_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(drug_display = factor(drug_display, levels = rev(as.character(display_drugs$drug_display))))

  p2 <- ggplot(mode_matrix, aes(x = death_mode_plot, y = drug_display, fill = log1p(pmid_count))) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_point(aes(size = induction_sensitization_count), shape = 21, fill = NA,
               colour = "grey25", stroke = 0.2) +
    scale_fill_gradient(low = "grey95", high = "#4C78A8", name = "log1p PMID count") +
    scale_size_continuous(range = c(0.6, 4), name = "Induction/sensitization candidates") +
    labs(title = "B  Drug x death-mode support profile", x = NULL, y = NULL) +
    theme_oncology_atlas(base_size = 6.4) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  plot <- p1 / p2 +
    plot_layout(heights = c(0.9, 1.15)) +
    plot_annotation(
      title = "Same drugs recur across multiple RCD-associated oncology literatures",
      subtitle = "Restored polydeath-index view. This supports hypothesis generation only and does not prove multi-pathway induction.",
      theme = theme(plot.title = element_text(face = "bold", size = 10),
                    plot.subtitle = element_text(size = 8, colour = "grey30"))
    )

  save_pub_r(plot, file.path(fig_dir, "Figure_same_drug_multi_death_signals"),
             width_mm = 183, height_mm = 165)
  save_pub_r(plot, file.path(fig_dir, "Figure_same_drug_multi_death_signals_previous_polydeath"),
             width_mm = 183, height_mm = 165)
}

copy_export_set(
  file.path(fig_dir, "Figure_oncology_drug_death_triads"),
  file.path(fig_dir, "Figure_oncology_drug_death_triads_by_death_mode")
)
copy_export_set(
  file.path(fig_dir, "Figure_same_drug_multi_death_signals"),
  file.path(fig_dir, "Figure_same_drug_multi_death_induction_sensitization")
)

make_previous_global_triads()
make_previous_polydeath()

write_text_safe(c(
  "# Restored Previous Drug Figures",
  "",
  "- The currently stratified-by-death-mode triad figure was preserved as `Figure_oncology_drug_death_triads_by_death_mode.*`.",
  "- The current induction/sensitization-focused same-drug figure was preserved as `Figure_same_drug_multi_death_induction_sensitization.*`.",
  "- The canonical `Figure_oncology_drug_death_triads.*` export was restored to a global top-triad view.",
  "- The canonical `Figure_same_drug_multi_death_signals.*` export was restored to a polydeath-index view.",
  "- Restored figures remain literature-level diagnostics and should not be used to claim validated drug-induced mechanisms."
), file.path(ROOT, "15_QC_reports/restored_previous_drug_figures_note.md"))

message("Restored previous drug figures and preserved revised versions with explicit suffixes.")
