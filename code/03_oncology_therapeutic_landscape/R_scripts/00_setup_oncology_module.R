suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

resolve_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script_path) == 0) {
    return(normalizePath("reviewer_revision/15_oncology_therapeutic_death_landscape", mustWork = FALSE))
  }
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}

ROOT <- resolve_root()
PROJECT_ROOT <- normalizePath(
  Sys.getenv("CELLDEATH_ATLAS_ROOT", unset = file.path(ROOT, "..", "..")),
  mustWork = FALSE
)
V5_ROOT <- Sys.getenv(
  "CELLDEATH_ATLAS_PACKAGE_DIR",
  unset = file.path(
    PROJECT_ROOT,
    Sys.getenv("CELLDEATH_ATLAS_PACKAGE", unset = "selected_death_modes_2000_2025")
  )
)
SELECTED_DEATH_MODES <- c(
  "Ferroptosis",
  "Pyroptosis",
  "NETosis",
  "Necroptosis",
  "Immunogenic cell death",
  "Cuproptosis",
  "PANoptosis",
  "Disulfidptosis"
)

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_csv_safe <- function(x, path) {
  dir_create(dirname(path))
  readr::write_csv(x, path, na = "")
}

write_text_safe <- function(lines, path) {
  dir_create(dirname(path))
  writeLines(lines, path, useBytes = TRUE)
}

zscore <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(0, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

safe_div <- function(num, den, eps = 0) {
  ifelse(is.na(den) | den == 0, NA_real_, (num + eps) / (den + eps))
}

theme_oncology_atlas <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "grey15"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey20"),
      axis.text = element_text(colour = "grey15"),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = grid::unit(3.5, "mm"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "grey30"),
      plot.caption = element_text(size = base_size - 1, colour = "grey35", hjust = 0),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "grey90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

death_mode_palette <- c(
  "Ferroptosis" = "#B9553C",
  "Pyroptosis" = "#D28A2E",
  "NETosis" = "#5C8F7B",
  "Necroptosis" = "#6B6BAE",
  "Immunogenic cell death" = "#A05B8F",
  "ICD" = "#A05B8F",
  "Cuproptosis" = "#3E7CB1",
  "PANoptosis" = "#1B9E77",
  "Disulfidptosis" = "#9C6A3D"
)

therapy_palette <- c(
  "chemotherapy" = "#4C78A8",
  "targeted_therapy" = "#F58518",
  "immunotherapy_checkpoint" = "#54A24B",
  "radiotherapy" = "#B279A2",
  "ferroptosis_inducer_tool_compound" = "#E45756",
  "kinase_inhibitor" = "#72B7B2",
  "proteasome_inhibitor" = "#9D755D",
  "metabolic_modulator" = "#BAB0AC",
  "natural_product" = "#8CD17D",
  "experimental_tool_compound" = "#FF9DA6",
  "unclear_or_mixed" = "#B9B9B9"
)

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 125, dpi = 600) {
  dir_create(dirname(filename))
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  ggsave(paste0(filename, ".pdf"), plot = plot, width = w, height = h, device = cairo_pdf, family = "Arial")
  ggsave(paste0(filename, ".png"), plot = plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(filename, ".tiff"), plot = plot, width = w, height = h, dpi = dpi, bg = "white", compression = "lzw")
}

read_csv_required <- function(path) {
  if (!file.exists(path)) stop("Required input missing: ", path)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

clean_bool <- function(x) {
  lx <- tolower(as.character(x))
  lx %in% c("true", "yes", "1", "y")
}

norm_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  str_squish(x)
}

regex_any <- function(patterns) {
  paste0("(", paste(patterns, collapse = "|"), ")")
}
