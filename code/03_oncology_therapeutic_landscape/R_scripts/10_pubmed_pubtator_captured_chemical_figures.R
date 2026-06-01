script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
dir_create(fig_dir)
dir_create(source_dir)

captured_path <- file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions.csv")
if (!file.exists(captured_path)) {
  stop("Missing PubMed/PubTator merged drug/chemical table: ", captured_path)
}

death_order <- SELECTED_DEATH_MODES

captured <- read_csv_required(captured_path) %>%
  mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    mention_source = coalesce(mention_source, ""),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
    death_mode = factor(death_mode, levels = death_order),
    mapping_status = case_when(
      is_specific_drug_yes_no == "yes" ~ "specific mapped drug",
      is_generic_therapy_term_yes_no == "yes" | is_specific_drug_yes_no == "no" ~ "generic therapy term",
      TRUE ~ "needs-review chemical/entity"
    ),
    captured_source_group = case_when(
      str_detect(mention_source, "PubMed_MeSH_ChemicalList") & str_detect(mention_source, "PubTator3_chemical") ~ "PubMed + PubTator",
      str_detect(mention_source, "PubMed_MeSH_ChemicalList") ~ "PubMed MeSH only",
      str_detect(mention_source, "PubTator3_chemical") ~ "PubTator only",
      TRUE ~ "not PubMed/PubTator"
    )
  ) %>%
  filter(!is.na(death_mode), captured_source_group != "not PubMed/PubTator") %>%
  distinct(pmid, tumor_family, disease_term, death_mode, drug_normalized, .keep_all = TRUE)

if (nrow(captured) == 0) {
  stop("No PubMed/PubTator captured chemical/drug candidates were available after filtering.")
}

status_levels <- c("specific mapped drug", "generic therapy term", "needs-review chemical/entity")
status_palette <- c(
  "specific mapped drug" = "#B9553C",
  "generic therapy term" = "#D28A2E",
  "needs-review chemical/entity" = "#8E8E8E"
)

panel_a <- captured %>%
  count(death_mode, mapping_status, name = "captured_contexts") %>%
  complete(
    death_mode = factor(death_order, levels = death_order),
    mapping_status = status_levels,
    fill = list(captured_contexts = 0)
  ) %>%
  mutate(
    mapping_status = factor(mapping_status, levels = status_levels),
    log10_contexts = log10(captured_contexts + 1)
  )

entity_summary <- captured %>%
  group_by(drug_normalized, drug_display) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    death_mode_breadth = n_distinct(death_mode),
    disease_count = n_distinct(disease_term),
    tumor_family_breadth = n_distinct(tumor_family),
    mapping_status = names(sort(table(mapping_status), decreasing = TRUE))[1],
    captured_sources = paste(sort(unique(captured_source_group)), collapse = "; "),
    death_modes = paste(sort(unique(as.character(death_mode))), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    ranking_score = death_mode_breadth * log10(pmid_count + 1),
    entity_label = ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display)
  ) %>%
  arrange(desc(ranking_score), desc(pmid_count), entity_label)

panel_b <- entity_summary %>%
  slice_head(n = 40) %>%
  mutate(entity_label = factor(entity_label, levels = rev(entity_label)))

panel_c <- captured %>%
  filter(death_mode == "Disulfidptosis") %>%
  group_by(drug_normalized, drug_display) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    disease_count = n_distinct(disease_term),
    tumor_family_breadth = n_distinct(tumor_family),
    mapping_status = names(sort(table(mapping_status), decreasing = TRUE))[1],
    captured_sources = paste(sort(unique(captured_source_group)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(entity_label = ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display)) %>%
  arrange(desc(pmid_count), entity_label) %>%
  slice_head(n = 35) %>%
  mutate(entity_label = factor(entity_label, levels = rev(entity_label)))

example_terms <- c("auranofin", "gaudichaudione h")
panel_d <- captured %>%
  filter(drug_normalized %in% example_terms | str_to_lower(drug_display) %in% example_terms) %>%
  mutate(
    entity_label = case_when(
      drug_normalized == "auranofin" ~ "Auranofin",
      drug_normalized == "gaudichaudione h" ~ "Gaudichaudione H",
      TRUE ~ ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display)
    )
  ) %>%
  group_by(drug_normalized, entity_label, death_mode, mapping_status) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    disease_count = n_distinct(disease_term),
    tumor_family_breadth = n_distinct(tumor_family),
    captured_sources = paste(sort(unique(captured_source_group)), collapse = "; "),
    .groups = "drop"
  )

write_csv_safe(captured, file.path(source_dir, "PubMed_PubTator_captured_chemical_candidate_contexts_unfiltered.csv"))
write_csv_safe(panel_a, file.path(source_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_panel_A_mapping_status.csv"))
write_csv_safe(entity_summary, file.path(source_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_all_entity_summary.csv"))
write_csv_safe(panel_b, file.path(source_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_panel_B_top40_broad.csv"))
write_csv_safe(panel_c, file.path(source_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_panel_C_disulfidptosis_top35.csv"))
write_csv_safe(panel_d, file.path(source_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_panel_D_user_checked_examples.csv"))

p_a <- ggplot(panel_a, aes(x = log10_contexts, y = death_mode, fill = mapping_status)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.22, colour = "white", linewidth = 0.15) +
  scale_x_continuous(breaks = 0:5) +
  scale_fill_manual(values = status_palette, drop = FALSE, name = NULL) +
  labs(
    title = "A  PubMed/PubTator captured candidates by death mode",
    x = "log10(captured PMID-chemical contexts + 1)",
    y = NULL
  ) +
  theme_oncology_atlas()

p_b <- ggplot(panel_b, aes(x = pmid_count, y = entity_label, fill = death_mode_breadth)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_x_continuous(labels = comma) +
  scale_fill_gradient(low = "#CFE3E0", high = "#1B6A6F", name = "Death-mode breadth") +
  labs(
    title = "B  Top 40 broad captured chemical/entity signals",
    subtitle = "Ranked by death-mode breadth x log10(PMID support)",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.4)

p_c <- ggplot(panel_c, aes(x = pmid_count, y = entity_label, fill = mapping_status)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = status_palette, drop = FALSE, name = NULL) +
  labs(
    title = "C  Disulfidptosis captured candidates",
    subtitle = "Unfiltered PubMed/PubTator chemical/entity layer",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.4)

if (nrow(panel_d) > 0) {
  p_d <- ggplot(panel_d, aes(x = pmid_count, y = entity_label, fill = death_mode)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.15) +
    facet_wrap(~ death_mode, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = death_mode_palette, guide = "none") +
    scale_x_continuous(labels = comma) +
    labs(
      title = "D  User-checked examples are captured before curation",
      subtitle = "Auranofin and Gaudichaudione H in PubMed/PubTator outputs",
      x = "Unique PMIDs",
      y = NULL
    ) +
    theme_oncology_atlas(base_size = 6.4) +
    theme(strip.text = element_text(size = 5.8))
} else {
  p_d <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Auranofin / Gaudichaudione H not found in PubMed/PubTator captured table", size = 3) +
    labs(title = "D  User-checked examples") +
    theme_void(base_family = "Arial")
}

main_plot <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(0.95, 1.35), heights = c(1, 1.15)) +
  plot_annotation(
    title = "Broad PubMed/PubTator captured chemical/entity layer in oncology death-mode literature",
    subtitle = "This is an intentionally unfiltered candidate view. Needs-review chemical/entity terms are not interpreted as specific therapies.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, colour = "grey30")
    )
  )

base <- file.path(fig_dir, "Figure_pubmed_pubtator_captured_chemical_landscape_unfiltered")
save_pub_r(main_plot, base, width_mm = 210, height_mm = 180)

focused_plot <- p_c / p_d +
  plot_layout(heights = c(1.25, 0.75)) +
  plot_annotation(
    title = "Disulfidptosis-associated PubMed/PubTator captured candidates",
    subtitle = "No curated-drug filter applied; candidate status should be resolved during downstream manual screening.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, colour = "grey30")
    )
  )

focused_base <- file.path(fig_dir, "Figure_disulfidptosis_pubmed_pubtator_captured_candidates_unfiltered")
save_pub_r(focused_plot, focused_base, width_mm = 183, height_mm = 150)

summary_lines <- c(
  "# PubMed/PubTator Captured Chemical Candidate Figure Summary",
  "",
  paste0("- PubMed/PubTator captured contexts: ", comma(nrow(captured))),
  paste0("- Unique captured PMIDs: ", comma(n_distinct(captured$pmid))),
  paste0("- Unique captured chemical/entity labels: ", comma(n_distinct(captured$drug_normalized))),
  paste0("- Specific mapped drug contexts: ", comma(sum(captured$is_specific_drug_yes_no == "yes"))),
  paste0("- Generic therapy-term contexts: ", comma(sum(captured$is_generic_therapy_term_yes_no == "yes" | captured$is_specific_drug_yes_no == "no"))),
  paste0("- Needs-review chemical/entity contexts: ", comma(sum(captured$mapping_status == "needs-review chemical/entity"))),
  paste0("- Disulfidptosis captured contexts: ", comma(sum(captured$death_mode == "Disulfidptosis"))),
  paste0("- User-checked example rows: ", comma(nrow(panel_d))),
  "",
  "Interpretation boundary: these figures show PubMed/PubTator-captured chemical/entity candidates before manual therapy filtering. They should not be interpreted as validated drug-death mechanisms."
)
write_text_safe(summary_lines, file.path(ROOT, "15_QC_reports/pubmed_pubtator_captured_chemical_candidate_figure_summary.md"))

message("Generated: ", base, ".pdf")
message("Generated: ", focused_base, ".pdf")
