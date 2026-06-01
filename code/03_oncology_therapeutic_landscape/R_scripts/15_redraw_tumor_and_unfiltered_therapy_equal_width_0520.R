script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
out_dir <- file.path(fig_dir, "new_figure_0520")
dir_create(out_dir)

death_order <- SELECTED_DEATH_MODES

target_width_mm <- 135

td_path <- file.path(source_dir, "Figure_tumor_family_death_mode_selectivity_latest_source.csv")
td_plot <- read_csv_required(td_path) %>%
  mutate(
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
    death_mode = factor(death_mode, levels = death_order)
  )

tumor_levels <- td_plot %>%
  group_by(tumor_family) %>%
  summarise(total_pair_count = sum(pair_count, na.rm = TRUE), .groups = "drop") %>%
  arrange(total_pair_count) %>%
  pull(tumor_family)

td_plot <- td_plot %>%
  mutate(tumor_family = factor(tumor_family, levels = tumor_levels))

p_tumor <- ggplot(td_plot, aes(x = death_mode, y = tumor_family, fill = pmax(pmin(log2_tumor_death_selectivity, 3), -3))) +
  geom_tile(colour = "white", linewidth = 0.24) +
  geom_point(aes(size = pair_count), shape = 21, stroke = 0.22, colour = "grey20", fill = NA) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#B44D3A", midpoint = 0, limits = c(-3, 3), name = "log2 selectivity") +
  scale_size_continuous(range = c(0.5, 3.6), name = "Pair count") +
  labs(
    title = "Tumor-family x death-mode selectivity",
    subtitle = "Normalized to oncology background; literature signal only.",
    x = NULL,
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 5.6) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 4.9),
    axis.text.y = element_text(size = 4.9),
    legend.title = element_text(size = 5.1),
    legend.text = element_text(size = 4.8),
    plot.title = element_text(size = 7.0, face = "bold"),
    plot.subtitle = element_text(size = 5.5, colour = "grey30"),
    legend.key.size = grid::unit(2.7, "mm"),
    plot.margin = margin(4, 5, 4, 4)
  )

save_pub_r(
  p_tumor,
  file.path(out_dir, "Figure_tumor_family_death_mode_selectivity_latest_equal_width_narrow"),
  width_mm = target_width_mm,
  height_mm = 105
)

class_path <- file.path(source_dir, "Figure_therapy_class_death_mode_landscape_unfiltered_reference_mapped_source.csv")
class_mode <- read_csv_required(class_path) %>%
  mutate(
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
    death_mode = factor(death_mode, levels = death_order)
  )

class_order <- c(
  "Approved clinical drug",
  "FDA product/brand",
  "Natural product candidate",
  "Local curated specific drug",
  "Generic therapy term",
  "Unmapped PubMed/PubTator entity"
)

class_palette <- c(
  "Approved clinical drug" = "#B9553C",
  "FDA product/brand" = "#D28A2E",
  "Natural product candidate" = "#5C8F7B",
  "Local curated specific drug" = "#6B6BAE",
  "Generic therapy term" = "#A05B8F",
  "Unmapped PubMed/PubTator entity" = "#8E8E8E"
)

class_mode <- class_mode %>%
  mutate(
    reference_class_label = factor(reference_class_label, levels = rev(class_order)),
    clipped_selectivity = pmax(pmin(clipped_selectivity, 3), -3)
  )

p_class_a <- ggplot(class_mode, aes(x = death_mode, y = reference_class_label, fill = clipped_selectivity)) +
  geom_tile(colour = "white", linewidth = 0.26) +
  geom_point(
    data = class_mode %>% filter(pmid_count > 0),
    aes(size = log1p(pmid_count)),
    shape = 21,
    colour = "grey20",
    fill = NA,
    stroke = 0.22
  ) +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#B44D3A",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey92",
    name = "log2 selectivity"
  ) +
  scale_size_continuous(range = c(0.6, 3.4), name = "log1p PMID") +
  labs(title = "A  Reference class x death-mode candidate landscape", x = NULL, y = NULL) +
  theme_oncology_atlas(base_size = 5.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 4.8),
    axis.text.y = element_text(size = 4.8),
    legend.title = element_text(size = 5.0),
    legend.text = element_text(size = 4.7),
    plot.title = element_text(size = 6.2, face = "bold"),
    legend.key.size = grid::unit(2.7, "mm")
  )

p_class_b <- class_mode %>%
  group_by(reference_class_label) %>%
  summarise(captured_contexts = sum(candidate_context_count, na.rm = TRUE), .groups = "drop") %>%
  mutate(reference_class_label = factor(reference_class_label, levels = rev(class_order))) %>%
  ggplot(aes(x = captured_contexts, y = reference_class_label, fill = reference_class_label)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = class_palette, drop = FALSE) +
  labs(title = "B  Captured contexts", x = "Contexts (log10)", y = NULL) +
  theme_oncology_atlas(base_size = 5.5) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 4.8),
    plot.title = element_text(size = 6.2, face = "bold")
  )

p_class <- p_class_a | p_class_b +
  plot_layout(widths = c(1.2, 0.55), guides = "collect") +
  plot_annotation(
    title = "Reference-mapped PubMed/PubTator candidate landscape",
    subtitle = "Unfiltered captured candidate view; not restricted to approved drugs.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 7.0),
      plot.subtitle = element_text(size = 5.5, colour = "grey30")
    )
  ) &
  theme(legend.position = "bottom")

save_pub_r(
  p_class,
  file.path(out_dir, "Figure_therapy_class_death_mode_landscape_unfiltered_reference_mapped_equal_width_narrow"),
  width_mm = target_width_mm,
  height_mm = 95
)

p_class_a_standalone <- p_class_a +
  labs(
    title = "Reference class x death-mode candidate landscape",
    subtitle = "Unfiltered PubMed/PubTator captured candidate view; not restricted to approved drugs."
  ) +
  theme(
    plot.title = element_text(size = 7.0, face = "bold"),
    plot.subtitle = element_text(size = 5.5, colour = "grey30"),
    legend.position = "right"
  )

save_pub_r(
  p_class_a_standalone,
  file.path(out_dir, "Figure_reference_class_death_mode_candidate_landscape_panelA_same_aspect_ratio_narrow"),
  width_mm = target_width_mm,
  height_mm = 105
)

p_class_a_tile_matched <- p_class_a +
  labs(
    title = "Reference class x death-mode candidate landscape",
    subtitle = NULL
  ) +
  theme(
    plot.title = element_text(size = 6.6, face = "bold"),
    legend.position = "right",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 4.4),
    axis.text.y = element_text(size = 4.4),
    legend.title = element_text(size = 4.7),
    legend.text = element_text(size = 4.4),
    legend.key.size = grid::unit(2.3, "mm"),
    plot.margin = margin(3, 4, 3, 3)
  )

save_pub_r(
  p_class_a_tile_matched,
  file.path(out_dir, "Figure_reference_class_death_mode_candidate_landscape_panelA_tile_matched_narrow"),
  width_mm = target_width_mm,
  height_mm = 55
)

write_text_safe(
  c(
    "# Equal-Width Narrow Figure Redraw Note",
    "",
    paste0("- Target width: ", target_width_mm, " mm for both figures."),
    "- Source data were reused from existing latest/unfiltered figure source CSVs.",
    "- Original figures were not overwritten.",
    "- Outputs are in `10_publication_grade_figures/new_figure_0520`.",
    "- A standalone Panel A reference-class landscape was also exported with the same 135 x 105 mm aspect ratio as the tumor-family selectivity figure.",
    "- A tile-matched Panel A was exported at 135 x 55 mm so its heatmap cells are closer in size to the 18-row tumor-family matrix when the two panels are assembled together."
  ),
  file.path(out_dir, "equal_width_narrow_tumor_therapy_redraw_note.md")
)

message("Generated equal-width narrow tumor and unfiltered therapy-class figures.")
