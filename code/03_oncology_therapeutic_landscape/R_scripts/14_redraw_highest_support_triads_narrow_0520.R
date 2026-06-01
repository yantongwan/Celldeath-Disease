script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

source_dir <- file.path(ROOT, "10_publication_grade_figures", "source_data")
out_dir <- file.path(ROOT, "10_publication_grade_figures", "new_figure_0520")
dir_create(out_dir)

death_order <- SELECTED_DEATH_MODES
death_short <- function(x) recode(as.character(x), "Immunogenic cell death" = "ICD", .default = as.character(x))

source_path <- file.path(source_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined_source.csv")
triad_summary <- read_csv_required(source_path)

top_triads <- triad_summary %>%
  arrange(desc(pmid_count), desc(induction_sensitization_pmids), desc(high_confidence_pmids), desc(context_count)) %>%
  slice_head(n = 32) %>%
  mutate(
    death_mode = factor(as.character(death_mode), levels = death_order),
    triad_label = paste(tumor_family, canonical_drug_display, death_short(death_mode), sep = " | "),
    triad_label = str_trunc(triad_label, 56),
    triad_label = factor(triad_label, levels = rev(unique(triad_label)))
  )

p_triads <- ggplot(top_triads, aes(x = pmid_count, y = triad_label, fill = death_mode)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_point(
    aes(x = induction_sensitization_pmids),
    shape = 21,
    fill = "white",
    colour = "grey15",
    stroke = 0.22,
    size = 1.55
  ) +
  scale_fill_manual(values = death_mode_palette, name = "Death mode") +
  labs(
    title = "Highest-support tumor-drug-death triads",
    subtitle = "Refined approved clinical drugs only; dot = induction/sensitization candidate PMIDs.",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 5.6) +
  theme(
    legend.position = "bottom",
    legend.key.size = grid::unit(3.0, "mm"),
    plot.title = element_text(size = 7.2, face = "bold"),
    plot.subtitle = element_text(size = 5.8, colour = "grey30"),
    axis.text.y = element_text(size = 4.9),
    axis.text.x = element_text(size = 5.2),
    plot.margin = margin(4, 5, 4, 4)
  )

base <- file.path(out_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined")
save_pub_r(p_triads, base, width_mm = 135, height_mm = 150)

write_csv_safe(
  top_triads,
  file.path(out_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined_narrow_source.csv")
)

write_text_safe(
  c(
    "# Narrow Redraw Note",
    "",
    "- Source: `10_publication_grade_figures/source_data/Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined_source.csv`.",
    "- Output folder: `10_publication_grade_figures/new_figure_0520`.",
    "- Width changed from the previous 165 mm export to 135 mm.",
    "- Labels were truncated to 56 characters and text was slightly reduced to preserve readability.",
    "- No original figure files were overwritten."
  ),
  file.path(out_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined_narrow_note.md")
)

message("Generated narrow highest-support tumor-drug-death triads figure.")
