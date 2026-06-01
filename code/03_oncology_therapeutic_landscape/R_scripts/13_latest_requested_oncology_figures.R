script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
qc_dir <- file.path(ROOT, "15_QC_reports")
ref_dir <- file.path(ROOT, "05_drug_therapy_extraction", "drug_reference")
dir_create(fig_dir)
dir_create(source_dir)
dir_create(qc_dir)

death_order <- SELECTED_DEATH_MODES
death_short <- function(x) recode(as.character(x), "Immunogenic cell death" = "ICD", .default = as.character(x))

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
  "experimental_tool_compound" = "Experimental/tool compound",
  "unmapped_chemical" = "Approved drug not locally classed",
  "unclear_or_mixed" = "Unclear/mixed"
)

class_palette_latest <- c(
  therapy_palette,
  "unmapped_chemical" = "#8A8A8A"
)

approved_path <- file.path(ref_dir, "pubmed_pubtator_approved_clinical_drug_refined_contexts.csv")
if (!file.exists(approved_path)) {
  stop(
    "Missing refined approved clinical drug table: ", approved_path,
    "\nRun 13_R_scripts/12_approved_clinical_drug_refined_figures.R first."
  )
}

approved <- read_csv_required(approved_path) %>%
  mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
    death_mode = factor(death_mode, levels = death_order),
    canonical_drug_display = coalesce(canonical_drug_display, drug_display, drug_normalized),
    therapy_class = ifelse(is.na(therapy_class) | therapy_class == "", "unmapped_chemical", therapy_class),
    therapy_class = ifelse(therapy_class %in% names(class_labels), therapy_class, "unmapped_chemical")
  ) %>%
  filter(!is.na(death_mode), drug_normalized != "")

relations_path <- file.path(ROOT, "09_sentence_level_relation_audit", "drug_death_relation_sentences.csv")
relations <- if (file.exists(relations_path)) {
  read_csv_required(relations_path) %>%
    mutate(
      pmid = as.character(pmid),
      death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
      relation_key = paste(pmid, death_mode, drug_normalized, sep = "||"),
      induction_sensitization_flag = relation_category %in% c("induction_or_activation", "sensitization") &
        confidence_high_medium_low %in% c("high", "medium")
    )
} else {
  tibble(relation_key = character(), induction_sensitization_flag = logical(), confidence_high_medium_low = character())
}

approved_relation_keys <- approved %>%
  mutate(relation_key = paste(pmid, as.character(death_mode), drug_normalized, sep = "||")) %>%
  select(relation_key, pmid, tumor_family, death_mode, drug_normalized, canonical_drug_display, therapy_class) %>%
  distinct()

relations_by_key <- relations %>%
  group_by(relation_key) %>%
  summarise(
    induction_sensitization_flag = any(induction_sensitization_flag, na.rm = TRUE),
    confidence_high_medium_low = case_when(
      any(confidence_high_medium_low == "high", na.rm = TRUE) ~ "high",
      any(confidence_high_medium_low == "medium", na.rm = TRUE) ~ "medium",
      any(confidence_high_medium_low == "low", na.rm = TRUE) ~ "low",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

approved_relations <- approved_relation_keys %>%
  inner_join(relations_by_key, by = "relation_key")

approved_relation_summary <- approved_relations %>%
  group_by(tumor_family, death_mode, drug_normalized, canonical_drug_display) %>%
  summarise(
    relation_candidate_pmids = n_distinct(pmid),
    induction_sensitization_pmids = n_distinct(pmid[induction_sensitization_flag]),
    high_confidence_pmids = n_distinct(pmid[confidence_high_medium_low == "high"]),
    .groups = "drop"
  )

approved_death_summary <- approved %>%
  group_by(death_mode) %>%
  summarise(
    approved_drug_pmid_count = n_distinct(pmid),
    approved_drug_context_count = n(),
    approved_drug_count = n_distinct(drug_normalized),
    .groups = "drop"
  )

hot <- read_csv_required(file.path(ROOT, "03_oncology_death_mode_hotness", "oncology_death_mode_hotness.csv")) %>%
  mutate(death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode)) %>%
  left_join(approved_death_summary, by = "death_mode") %>%
  mutate(
    across(c(approved_drug_pmid_count, approved_drug_context_count, approved_drug_count), ~coalesce(.x, 0)),
    approved_drug_pmid_fraction = ifelse(oncology_unique_pmids > 0, approved_drug_pmid_count / oncology_unique_pmids, 0),
    heat_index_latest_approved_drug_refined =
      zscore(log1p(oncology_unique_pmids)) +
      zscore(tumor_family_breadth) +
      zscore(recent_growth_3yr) +
      zscore(approved_drug_pmid_fraction),
    heat_rank_latest_approved_drug_refined = min_rank(desc(heat_index_latest_approved_drug_refined)),
    death_mode_factor = factor(death_mode, levels = death_mode[order(heat_index_latest_approved_drug_refined)])
  )

write_csv_safe(
  hot,
  file.path(source_dir, "Figure_oncology_literature_heat_index_latest_approved_drug_refined_source.csv")
)

p_heat <- ggplot(hot, aes(x = death_mode_factor, y = heat_index_latest_approved_drug_refined, fill = death_mode)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = death_mode_palette, guide = "none") +
  labs(
    title = "Oncology literature heat index",
    subtitle = "Approved-drug-calibrated component-wise literature index.",
    x = NULL,
    y = "Component-wise heat index"
  ) +
  theme_oncology_atlas(base_size = 6.8)

save_pub_r(
  p_heat,
  file.path(fig_dir, "Figure_oncology_literature_heat_index_latest_approved_drug_refined"),
  width_mm = 145,
  height_mm = 88
)

hot_label <- hot %>%
  mutate(
    label_short = death_short(death_mode),
    death_mode = factor(death_mode, levels = death_order)
  )

write_csv_safe(
  hot_label,
  file.path(source_dir, "Figure_oncology_fraction_drug_vocabulary_latest_approved_drug_refined_source.csv")
)

p_fraction <- ggplot(
  hot_label,
  aes(
    x = fraction_of_death_mode_literature_in_oncology,
    y = approved_drug_pmid_fraction,
    colour = death_mode
  )
) +
  geom_point(aes(size = approved_drug_count), alpha = 0.9) +
  geom_text(aes(label = label_short), nudge_x = 0.018, hjust = 0, size = 2.25, show.legend = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(0.9, max(hot_label$approved_drug_pmid_fraction, na.rm = TRUE) * 1.12))) +
  scale_colour_manual(values = death_mode_palette, guide = "none") +
  scale_size_continuous(range = c(2.2, 6), name = "Approved drugs") +
  labs(
    title = "Oncology fraction and approved-drug vocabulary",
    subtitle = "Drug vocabulary is limited to refined approved clinical drugs; common ions/metabolites are audit-only.",
    x = "Fraction of death-mode records in oncology",
    y = "Approved-drug PMID fraction"
  ) +
  coord_cartesian(clip = "off") +
  theme_oncology_atlas(base_size = 7.2) +
  theme(plot.margin = margin(5.5, 25, 5.5, 5.5), legend.position = "bottom")

save_pub_r(
  p_fraction,
  file.path(fig_dir, "Figure_oncology_fraction_drug_vocabulary_latest_approved_drug_refined"),
  width_mm = 120,
  height_mm = 88
)

triad_summary <- approved %>%
  group_by(tumor_family, death_mode, drug_normalized, canonical_drug_display, therapy_class) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    context_count = n(),
    disease_count = n_distinct(disease_term),
    first_year = suppressWarnings(min(year, na.rm = TRUE)),
    last_year = suppressWarnings(max(year, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    approved_relation_summary,
    by = c("tumor_family", "death_mode", "drug_normalized", "canonical_drug_display")
  ) %>%
  mutate(
    across(c(relation_candidate_pmids, induction_sensitization_pmids, high_confidence_pmids), ~coalesce(.x, 0L)),
    death_mode_plot = factor(death_short(death_mode), levels = death_short(death_order)),
    therapy_class_label = class_labels[therapy_class],
    therapy_class_label = ifelse(is.na(therapy_class_label), "Approved drug not locally classed", therapy_class_label),
    triad_label = paste(tumor_family, canonical_drug_display, death_short(death_mode), sep = " | "),
    triad_label = str_trunc(triad_label, 66)
  ) %>%
  arrange(desc(pmid_count), desc(induction_sensitization_pmids), desc(high_confidence_pmids), desc(context_count))

write_csv_safe(
  triad_summary,
  file.path(source_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined_source.csv")
)

top_triads <- triad_summary %>%
  slice_head(n = 32) %>%
  mutate(
    triad_label = factor(triad_label, levels = rev(unique(triad_label))),
    death_mode = factor(as.character(death_mode), levels = death_order)
  )

p_triads <- ggplot(top_triads, aes(x = pmid_count, y = triad_label, fill = death_mode)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_point(aes(x = induction_sensitization_pmids), shape = 21, fill = "white", colour = "grey15", stroke = 0.25, size = 1.9) +
  scale_fill_manual(values = death_mode_palette, name = "Death mode") +
  labs(
    title = "Highest-support tumor-drug-death triads",
    subtitle = "Refined approved clinical drugs only; dot = induction/sensitization sentence-candidate PMIDs.",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.4) +
  theme(legend.position = "bottom")

save_pub_r(
  p_triads,
  file.path(fig_dir, "Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined"),
  width_mm = 165,
  height_mm = 150
)

td <- read_csv_required(file.path(ROOT, "04_tumor_specific_death_profiles", "tumor_death_selectivity.csv")) %>%
  mutate(death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode))

top_tumors <- td %>%
  group_by(tumor_family) %>%
  summarise(total = sum(pair_count), breadth = sum(pair_count >= 3), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 18) %>%
  pull(tumor_family)

td_plot <- td %>%
  filter(tumor_family %in% top_tumors) %>%
  mutate(
    tumor_family = factor(tumor_family, levels = rev(top_tumors)),
    death_mode = factor(death_mode, levels = death_order)
  )

write_csv_safe(
  td_plot,
  file.path(source_dir, "Figure_tumor_family_death_mode_selectivity_latest_source.csv")
)

p_tumor <- ggplot(td_plot, aes(x = death_mode, y = tumor_family, fill = pmax(pmin(log2_tumor_death_selectivity, 3), -3))) +
  geom_tile(colour = "white", linewidth = 0.28) +
  geom_point(aes(size = pair_count), shape = 21, stroke = 0.25, colour = "grey20", fill = NA) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#B44D3A", midpoint = 0, limits = c(-3, 3), name = "log2 selectivity") +
  scale_size_continuous(range = c(0.8, 4.8), name = "Pair count") +
  labs(
    title = "Tumor-family x death-mode selectivity",
    subtitle = "Selectivity is normalized to the oncology background; it is not proof of tumor biology.",
    x = NULL,
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

save_pub_r(
  p_tumor,
  file.path(fig_dir, "Figure_tumor_family_death_mode_selectivity_latest"),
  width_mm = 165,
  height_mm = 115
)

drug_mode_class <- approved %>%
  group_by(therapy_class, death_mode) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    context_count = n(),
    drug_count = n_distinct(drug_normalized),
    tumor_family_breadth = n_distinct(tumor_family),
    .groups = "drop"
  )

relation_class <- approved_relations %>%
  group_by(therapy_class, death_mode) %>%
  summarise(
    relation_candidate_pmids = n_distinct(pmid),
    induction_sensitization_pmids = n_distinct(pmid[induction_sensitization_flag]),
    high_confidence_pmids = n_distinct(pmid[confidence_high_medium_low == "high"]),
    .groups = "drop"
  )

class_mode <- drug_mode_class %>%
  left_join(relation_class, by = c("therapy_class", "death_mode")) %>%
  mutate(across(c(relation_candidate_pmids, induction_sensitization_pmids, high_confidence_pmids), ~coalesce(.x, 0L)))

class_order <- class_mode %>%
  group_by(therapy_class) %>%
  summarise(total = sum(pmid_count), .groups = "drop") %>%
  arrange(total) %>%
  pull(therapy_class)

class_totals <- class_mode %>% group_by(therapy_class) %>% summarise(class_total = sum(pmid_count), .groups = "drop")
death_totals <- class_mode %>% group_by(death_mode) %>% summarise(death_total = sum(pmid_count), .groups = "drop")
all_total <- sum(class_mode$pmid_count)

class_mode <- class_mode %>%
  left_join(class_totals, by = "therapy_class") %>%
  left_join(death_totals, by = "death_mode") %>%
  mutate(
    therapy_class_label = class_labels[therapy_class],
    therapy_class_label = ifelse(is.na(therapy_class_label), "Approved drug not locally classed", therapy_class_label),
    therapy_class_label = factor(therapy_class_label, levels = class_labels[class_order]),
    death_mode = factor(death_mode, levels = death_order),
    therapy_death_selectivity = log2((pmid_count / class_total + 1e-6) / (death_total / all_total + 1e-6)),
    clipped_selectivity = pmax(pmin(therapy_death_selectivity, 3), -3)
  ) %>%
  complete(
    therapy_class_label = factor(class_labels[class_order], levels = class_labels[class_order]),
    death_mode = factor(death_order, levels = death_order),
    fill = list(
      pmid_count = 0,
      context_count = 0,
      drug_count = 0,
      tumor_family_breadth = 0,
      relation_candidate_pmids = 0,
      induction_sensitization_pmids = 0,
      high_confidence_pmids = 0,
      clipped_selectivity = NA_real_
    )
  )

write_csv_safe(
  class_mode,
  file.path(source_dir, "Figure_therapy_class_death_mode_relation_landscape_latest_approved_drug_refined_source.csv")
)

p_class <- ggplot(class_mode, aes(x = death_mode, y = therapy_class_label, fill = clipped_selectivity)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data = class_mode %>% filter(induction_sensitization_pmids > 0),
    aes(size = induction_sensitization_pmids),
    shape = 21, colour = "grey20", fill = NA, stroke = 0.25
  ) +
  scale_fill_gradient2(
    low = "#4C78A8", mid = "white", high = "#B44D3A",
    midpoint = 0, limits = c(-3, 3), na.value = "grey92",
    name = "log2 selectivity"
  ) +
  scale_size_continuous(range = c(0.7, 4), name = "Induction/sensitization PMIDs") +
  labs(
    title = "Therapy class x death-mode relation landscape",
    subtitle = "Refined approved clinical drugs only; classes use the local oncology therapy class map where available.",
    x = NULL,
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

save_pub_r(
  p_class,
  file.path(fig_dir, "Figure_therapy_class_death_mode_relation_landscape_latest_approved_drug_refined"),
  width_mm = 165,
  height_mm = 105
)

summary_lines <- c(
  "# Latest Requested Oncology Figure Set",
  "",
  "## Scope",
  "These figures were regenerated as standalone panel-style exports. Drug-linked panels use the refined approved clinical drug layer and exclude common ion/nutrient/metabolite-like terms from plotted drug signals.",
  "",
  "## Generated figures",
  "- Figure_oncology_literature_heat_index_latest_approved_drug_refined",
  "- Figure_oncology_fraction_drug_vocabulary_latest_approved_drug_refined",
  "- Figure_highest_support_tumor_drug_death_triads_latest_approved_drug_refined",
  "- Figure_tumor_family_death_mode_selectivity_latest",
  "- Figure_therapy_class_death_mode_relation_landscape_latest_approved_drug_refined",
  "",
  "## Key counts",
  paste0("- Refined approved clinical drug contexts: ", comma(nrow(approved))),
  paste0("- Refined approved clinical drug PMIDs: ", comma(n_distinct(approved$pmid))),
  paste0("- Refined approved clinical drug entities: ", comma(n_distinct(approved$drug_normalized))),
  paste0("- Top tumor-drug-death triads exported: ", comma(nrow(top_triads))),
  paste0("- Tumor families in selectivity panel: ", comma(length(top_tumors))),
  paste0("- Therapy classes in relation landscape: ", comma(n_distinct(class_mode$therapy_class_label))),
  "",
  "## Interpretation boundary",
  "All panels are literature-signal summaries. They do not establish clinical efficacy, drug-induced cell-death mechanisms, or tumor biology."
)

write_text_safe(summary_lines, file.path(qc_dir, "latest_requested_oncology_figure_set_summary.md"))

message("Generated latest requested oncology figure set.")
