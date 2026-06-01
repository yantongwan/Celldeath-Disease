script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
dir_create(fig_dir)
dir_create(source_dir)

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

tier_labels <- c(
  "Tier1_PubMed_MeSH_ChemicalList" = "Tier 1 PubMed MeSH ChemicalList",
  "Tier2_PubTator3_chemical" = "Tier 2 PubTator3 chemical",
  "Tier3_local_lexicon_fallback" = "Tier 3 title/abstract lexical fallback"
)

tier_palette <- c(
  "Tier 1 PubMed MeSH ChemicalList" = "#4C78A8",
  "Tier 2 PubTator3 chemical" = "#72B7B2",
  "Tier 3 title/abstract lexical fallback" = "#F58518"
)

status_palette <- c(
  "Specific mapped drug" = "#B9553C",
  "Generic therapy term" = "#D28A2E",
  "Needs-review chemical/entity" = "#8E8E8E"
)

death_order <- SELECTED_DEATH_MODES

mentions <- read_csv_required(file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions.csv")) %>%
  mutate(
    pmid = as.character(pmid),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
    death_mode = factor(death_mode, levels = death_order),
    evidence_tier = coalesce(cooccurrence_tier, ""),
    mapping_status = case_when(
      is_specific_drug_yes_no == "yes" ~ "Specific mapped drug",
      is_generic_therapy_term_yes_no == "yes" | is_specific_drug_yes_no == "no" ~ "Generic therapy term",
      TRUE ~ "Needs-review chemical/entity"
    ),
    context_key = paste(pmid, death_mode, disease_term, drug_normalized, sep = "||")
  ) %>%
  filter(!is.na(death_mode))

tier_source <- mentions %>%
  select(context_key, pmid, death_mode, drug_normalized, drug_display, mapping_status, evidence_tier) %>%
  tidyr::separate_rows(evidence_tier, sep = ";") %>%
  mutate(evidence_tier = str_squish(evidence_tier)) %>%
  filter(evidence_tier %in% names(tier_labels)) %>%
  mutate(evidence_tier_label = tier_labels[evidence_tier])

panel_a_data <- tier_source %>%
  group_by(death_mode, evidence_tier_label) %>%
  summarise(
    evidence_contexts = n_distinct(context_key),
    unique_pmids = n_distinct(pmid),
    unique_entities = n_distinct(drug_normalized),
    .groups = "drop"
  ) %>%
  mutate(
    death_mode = factor(death_mode, levels = death_order),
    evidence_tier_label = factor(evidence_tier_label, levels = tier_labels)
  )

panel_b_data <- mentions %>%
  group_by(death_mode, mapping_status) %>%
  summarise(
    contexts = n_distinct(context_key),
    unique_pmids = n_distinct(pmid),
    unique_entities = n_distinct(drug_normalized),
    .groups = "drop"
  ) %>%
  group_by(death_mode) %>%
  mutate(fraction = contexts / sum(contexts)) %>%
  ungroup() %>%
  mutate(
    death_mode = factor(death_mode, levels = death_order),
    mapping_status = factor(mapping_status, levels = names(status_palette))
  )

panel_c_data <- mentions %>%
  group_by(drug_normalized, drug_display, mapping_status) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    evidence_contexts = n_distinct(context_key),
    death_mode_breadth = n_distinct(death_mode),
    tier_count = n_distinct(unlist(str_split(paste(evidence_tier, collapse = ";"), ";"))),
    top_death_modes = paste(head(names(sort(table(as.character(death_mode)), decreasing = TRUE)), 5), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(pmid_count), desc(death_mode_breadth), desc(evidence_contexts)) %>%
  slice_head(n = 25) %>%
  mutate(
    entity_label = ifelse(nchar(drug_display) > 42, paste0(str_sub(drug_display, 1, 39), "..."), drug_display),
    entity_label = factor(entity_label, levels = rev(entity_label)),
    mapping_status = factor(mapping_status, levels = names(status_palette))
  )

mat <- read_csv_required(file.path(ROOT, "08_therapy_class_landscape/therapy_class_death_mode_matrix.csv")) %>%
  mutate(
    therapy_class_raw = therapy_class,
    therapy_class = class_labels[therapy_class_raw],
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
    death_mode = factor(death_mode, levels = death_order)
  ) %>%
  filter(!is.na(therapy_class), !is.na(death_mode))

class_order <- mat %>%
  group_by(therapy_class) %>%
  summarise(total = sum(pmid_count), .groups = "drop") %>%
  arrange(total) %>%
  pull(therapy_class)

panel_d_data <- mat %>%
  mutate(
    therapy_class = factor(therapy_class, levels = class_order),
    clipped_selectivity = pmax(pmin(therapy_death_selectivity, 3), -3)
  )

write_csv_safe(panel_a_data, file.path(source_dir, "Figure_therapy_tier123_broad_evidence_panel_A_tier_by_death_mode.csv"))
write_csv_safe(panel_b_data, file.path(source_dir, "Figure_therapy_tier123_broad_evidence_panel_B_mapping_status.csv"))
write_csv_safe(panel_c_data, file.path(source_dir, "Figure_therapy_tier123_broad_evidence_panel_C_top_entities.csv"))
write_csv_safe(panel_d_data, file.path(source_dir, "Figure_therapy_tier123_broad_evidence_panel_D_specific_drug_matrix.csv"))

p1 <- ggplot(panel_a_data, aes(x = evidence_contexts, y = death_mode, fill = evidence_tier_label)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.18) +
  scale_x_continuous(breaks = breaks_extended(n = 4), labels = function(x) paste0(round(x / 1000), "k")) +
  scale_fill_manual(values = tier_palette, name = NULL) +
  labs(
    title = "A  Tier 1/2/3 evidence support by death mode",
    x = "Evidence-supported PMID-chemical contexts",
    y = NULL
  ) +
  theme_oncology_atlas()

p2 <- ggplot(panel_b_data, aes(x = fraction, y = death_mode, fill = mapping_status)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.18) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = status_palette, name = NULL) +
  labs(
    title = "B  Mapping status of all chemical evidence",
    x = "Fraction of evidence contexts",
    y = NULL
  ) +
  theme_oncology_atlas()

p3 <- ggplot(panel_c_data, aes(x = pmid_count, y = entity_label, fill = mapping_status)) +
  geom_col(width = 0.72) +
  geom_point(aes(size = death_mode_breadth), shape = 21, colour = "grey20", fill = "white", stroke = 0.25) +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = status_palette, name = NULL) +
  scale_size_continuous(range = c(1.2, 3.6), name = "Death-mode breadth") +
  labs(
    title = "C  Top broad Tier-supported chemical/entity terms",
    x = "Unique oncology PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.4)

p4 <- ggplot(panel_d_data, aes(x = death_mode, y = therapy_class, fill = clipped_selectivity)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(aes(size = high_confidence_relation_count), shape = 21, colour = "grey20", fill = NA, stroke = 0.25) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#B9553C", midpoint = 0, limits = c(-3, 3), name = "log2 selectivity") +
  scale_size_continuous(range = c(0.6, 3.6), name = "High-confidence relations") +
  labs(
    title = "D  Conservative specific-drug therapy-class layer",
    x = NULL,
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.4) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Tiered chemical/drug evidence before conservative therapy-class restriction",
    subtitle = "Tier 1 = PubMed MeSH ChemicalList; Tier 2 = PubTator3 chemical; Tier 3 = title/abstract lexical fallback.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, colour = "grey30")
    )
  )

base <- file.path(fig_dir, "Figure_therapy_class_death_mode_tier123_broad_evidence")
save_pub_r(plot, base, width_mm = 183, height_mm = 160)

summary_lines <- c(
  "# Tier 1/2/3 Therapy Evidence Figure Summary",
  "",
  paste0("- Total evidence contexts: ", scales::comma(n_distinct(mentions$context_key))),
  paste0("- Unique PMIDs: ", scales::comma(n_distinct(mentions$pmid))),
  paste0("- Unique normalized chemical/entity terms: ", scales::comma(n_distinct(mentions$drug_normalized))),
  paste0("- Specific mapped drug contexts: ", scales::comma(sum(mentions$mapping_status == "Specific mapped drug"))),
  paste0("- Generic therapy-term contexts: ", scales::comma(sum(mentions$mapping_status == "Generic therapy term"))),
  paste0("- Needs-review chemical/entity contexts: ", scales::comma(sum(mentions$mapping_status == "Needs-review chemical/entity"))),
  "",
  "Interpretation boundary: Tier-supported chemical/entity annotations are literature-level evidence. Only the conservative specific-drug layer should be read as mapped therapy-class evidence, and even that remains a sentence-rule candidate rather than drug efficacy or mechanism proof."
)
write_text_safe(summary_lines, file.path(ROOT, "15_QC_reports/therapy_tier123_broad_evidence_figure_summary.md"))
