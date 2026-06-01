#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("CELLDEATH_ATLAS_ROOT", unset = getwd())
module_root <- file.path(project_root, "reviewer_revision", "16_gene_protein_death_disease_layer_v3")
extension_root <- file.path(module_root, "17_vosviewer_enrichment_extensions")
input_path <- file.path(extension_root, "06_system_broad_top20", "v3_system_disease_broad_multi_death_gene_signals_top20.csv")
out_root <- file.path(extension_root, "07_system_broad_top20_figures")

dirs <- c("pdf", "png", "tiff", "source_data", "per_system_pdf", "per_system_png", "15_QC_reports")
for (d in dirs) dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Missing required input: ", input_path)
}

system_top20 <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("disease_system", "shared_display_symbol", "pmid_count", "death_mode_breadth", "broad_score", "system_rank")
missing_cols <- setdiff(required_cols, colnames(system_top20))
if (length(missing_cols) > 0) {
  stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
}

theme_system <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#1F2933"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#1F2933"),
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size, colour = "#4B5563"),
      strip.background = element_rect(fill = "#F4F6F7", colour = "#D5DBDB", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = base_size - 0.2, colour = "#1F2933"),
      legend.position = "none",
      panel.spacing = unit(1.1, "lines")
    )
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "system", x)
}

lighten_colour <- function(colour, amount) {
  rgb_val <- grDevices::col2rgb(colour)
  mixed <- rgb_val + (255 - rgb_val) * amount
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ], maxColorValue = 255)
}

system_levels <- system_top20 %>%
  group_by(disease_system) %>%
  summarise(max_score = max(as.numeric(broad_score), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_score), disease_system) %>%
  pull(disease_system)

base_cols <- setNames(
  grDevices::hcl.colors(length(system_levels), palette = "Dark 3"),
  system_levels
)

plot_df <- system_top20 %>%
  mutate(
    pmid_count = as.numeric(pmid_count),
    death_mode_breadth = as.numeric(death_mode_breadth),
    broad_score = as.numeric(broad_score),
    system_rank = as.integer(system_rank),
    disease_system = factor(disease_system, levels = system_levels)
  ) %>%
  arrange(disease_system, system_rank) %>%
  group_by(disease_system) %>%
  mutate(
    system_base_colour = base_cols[as.character(disease_system)],
    breadth_scaled = {
      breadth_range <- range(death_mode_breadth, na.rm = TRUE)
      if (breadth_range[[1]] == breadth_range[[2]]) {
        rep(1, dplyr::n())
      } else {
        scales::rescale(death_mode_breadth, to = c(0, 1), from = breadth_range)
      }
    },
    fill_colour = mapply(
      function(colour, scaled) lighten_colour(colour, amount = 0.58 - 0.42 * scaled),
      system_base_colour,
      breadth_scaled,
      USE.NAMES = FALSE
    ),
    panel_gene_label = paste(shared_display_symbol, disease_system, sep = "___")
  ) %>%
  ungroup()

panel_levels <- plot_df %>%
  group_by(disease_system) %>%
  arrange(system_rank, .by_group = TRUE) %>%
  summarise(levels = list(rev(panel_gene_label)), .groups = "drop") %>%
  pull(levels) %>%
  unlist(use.names = FALSE)

plot_df <- plot_df %>%
  mutate(
    panel_gene_label = factor(panel_gene_label, levels = panel_levels),
    source_note = "Primary v3 shared human-mouse title/abstract/keyword co-mention"
  )

write.csv(
  plot_df,
  file.path(out_root, "source_data", "v3_system_broad_top20_figure_source_data.csv"),
  row.names = FALSE
)

label_lookup <- setNames(plot_df$shared_display_symbol, as.character(plot_df$panel_gene_label))
label_gene <- function(x) unname(label_lookup[as.character(x)])

combined_plot <- ggplot(plot_df, aes(x = panel_gene_label, y = pmid_count, fill = fill_colour)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.08) +
  geom_text(aes(label = paste0(death_mode_breadth, " modes")), hjust = -0.05, size = 1.45, colour = "#263238") +
  coord_flip() +
  facet_wrap(~ disease_system, scales = "free_y", ncol = 3) +
  scale_x_discrete(labels = label_gene) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.28))) +
  scale_fill_identity() +
  labs(
    title = "System-specific Top 20 broad multi-death gene signals",
    subtitle = "Ortholog-aware display; ranked within each disease system by breadth x log10(PMID support). Hue encodes disease system; darker shade indicates broader death-mode breadth.",
    x = NULL,
    y = "PMIDs"
  ) +
  theme_system(base_size = 6.2)

ggsave(
  file.path(out_root, "pdf", "Supplementary_Figure_System_Broad_Top20_v3_faceted.pdf"),
  combined_plot, width = 16, height = 24, units = "in", useDingbats = FALSE
)
ggsave(
  file.path(out_root, "png", "Supplementary_Figure_System_Broad_Top20_v3_faceted.png"),
  combined_plot, width = 16, height = 24, units = "in", dpi = 400
)
ggsave(
  file.path(out_root, "tiff", "Supplementary_Figure_System_Broad_Top20_v3_faceted.tiff"),
  combined_plot, width = 16, height = 24, units = "in", dpi = 400, compression = "lzw"
)

make_system_plot <- function(df, system_name) {
  single <- df %>%
    filter(as.character(disease_system) == system_name) %>%
    arrange(system_rank) %>%
    mutate(gene_factor = factor(shared_display_symbol, levels = rev(shared_display_symbol)))

  ggplot(single, aes(x = gene_factor, y = pmid_count, fill = fill_colour)) +
    geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
    geom_text(aes(label = paste0(death_mode_breadth, " modes")), hjust = -0.05, size = 2.0, colour = "#263238") +
    coord_flip() +
    scale_fill_identity() +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = paste0("Top 20 broad multi-death signals: ", system_name),
      subtitle = "Ortholog-aware display; ranked by breadth x log10(PMID support)",
      x = NULL,
      y = "PMIDs"
    ) +
    theme_system(base_size = 7)
}

for (system_name in system_levels) {
  p <- make_system_plot(plot_df, system_name)
  file_stub <- paste0("System_Broad_Top20_v3_", sanitize_filename(system_name))
  ggsave(file.path(out_root, "per_system_pdf", paste0(file_stub, ".pdf")), p, width = 7.2, height = 6.2, units = "in", useDingbats = FALSE)
  ggsave(file.path(out_root, "per_system_png", paste0(file_stub, ".png")), p, width = 7.2, height = 6.2, units = "in", dpi = 400)
}

qc <- data.frame(
  metric = c(
    "disease_systems",
    "rows",
    "rows_per_system_min",
    "rows_per_system_max",
    "combined_faceted_pdf",
    "per_system_figures"
  ),
  value = c(
    length(system_levels),
    nrow(plot_df),
    min(table(plot_df$disease_system)),
    max(table(plot_df$disease_system)),
    file.path(out_root, "pdf", "Supplementary_Figure_System_Broad_Top20_v3_faceted.pdf"),
    length(system_levels)
  )
)
write.csv(qc, file.path(out_root, "15_QC_reports", "v3_system_broad_top20_figure_QC.csv"), row.names = FALSE)

summary_lines <- c(
  "# v3 System-Specific Broad Multi-Death Top20 Figure Summary",
  "",
  paste0("- Disease systems plotted: ", length(system_levels)),
  paste0("- Total rows plotted: ", nrow(plot_df)),
  paste0("- Combined faceted figure: ", file.path(out_root, "pdf", "Supplementary_Figure_System_Broad_Top20_v3_faceted.pdf")),
  paste0("- Per-system figures: ", length(system_levels), " PDF/PNG pairs"),
  "",
  "Visual encoding: x-axis is PMID count; rows are ranked within each disease system by death-mode breadth x log10(PMID support). Hue encodes disease system, and darker shades within a system indicate broader death-mode breadth.",
  "",
  "Interpretation boundary: these are v3 shared human-mouse title/abstract/keyword co-mention signals, not causal gene-function evidence."
)
writeLines(summary_lines, file.path(out_root, "15_QC_reports", "v3_system_broad_top20_figure_summary.md"))

cat(paste(summary_lines, collapse = "\n"), "\n")
