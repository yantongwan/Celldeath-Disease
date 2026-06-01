script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
ref_dir <- file.path(ROOT, "05_drug_therapy_extraction", "drug_reference")
qc_dir <- file.path(ROOT, "15_QC_reports")
dir_create(fig_dir)
dir_create(source_dir)
dir_create(ref_dir)
dir_create(qc_dir)

normalize_drug_name <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- str_replace_all(x, "<[^>]+>", " ")
  x <- str_to_lower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_squish(x)
}

canonical_display <- function(x) {
  x <- str_squish(as.character(x))
  x <- ifelse(is.na(x) | x == "", NA_character_, x)
  x
}

death_order <- SELECTED_DEATH_MODES

mapped_path <- file.path(ref_dir, "pubmed_pubtator_candidates_reference_mapped.csv")
if (!file.exists(mapped_path)) {
  stop(
    "Missing reference-mapped PubMed/PubTator candidate table: ", mapped_path,
    "\nRun 13_R_scripts/11_drug_reference_mapped_unfiltered_figures.R first."
  )
}

mapped <- read_csv_required(mapped_path) %>%
  mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
    death_mode = factor(death_mode, levels = death_order),
    drug_normalized = normalize_drug_name(coalesce(drug_normalized, candidate_normalized, drug_display)),
    drug_display = canonical_display(drug_display),
    reference_display = canonical_display(reference_display),
    reference_class_label = coalesce(reference_class_label, "")
  ) %>%
  filter(!is.na(death_mode))

approved_raw <- mapped %>%
  filter(reference_class_label == "Approved clinical drug", drug_normalized != "")

# These terms are FDA/Drugs@FDA active ingredients or products, but in the oncology
# death-mode landscape they mostly behave as generic ions, nutrients, metabolites,
# solvents, excipients, or broad physiology terms rather than specific clinical drugs.
common_substance_terms <- c(
  "iron", "copper", "glucose", "oxygen", "calcium", "magnesium", "potassium", "sodium",
  "zinc", "manganese", "chloride", "phosphate", "sulfate", "sulphate", "water", "saline",
  "dextrose", "sucrose", "lactose", "glycerol", "ethanol", "edta", "atp", "adenosine",
  "glutathione", "cysteine", "cystine", "glycine", "alanine", "arginine", "lysine",
  "tryptophan", "tyrosine", "histidine", "methionine", "phenylalanine",
  "glutamine", "glutamic acid", "serine", "amino acids", "adenine",
  "choline", "citrulline", "asparagine", "aspartic acid",
  "hydrogen peroxide", "nitric oxide", "hyaluronic acid", "lactic acid",
  "silicon dioxide", "ferric ammonium citrate", "prussian blue", "folic acid",
  "oleic acid", "tannic acid", "glutathione disulfide", "alcohol", "dopamine",
  "nitrogen", "ammonia", "ammonium lactate", "acetic acid", "acetylcholine",
  "aluminum chloride", "aluminium chloride", "aluminum acetate", "aluminum hydroxide",
  "aluminium hydroxide", "alpha tocopherol acetate", "alpha-tocopherol acetate",
  "reactive oxygen species", "ros", "n acetylcysteine", "acetylcysteine",
  "vitamin c", "ascorbic acid", "ascorbate", "vitamin e", "tocopherol"
)

canonical_lookup <- approved_raw %>%
  mutate(display_candidate = coalesce(reference_display, drug_display, drug_normalized)) %>%
  count(drug_normalized, display_candidate, sort = TRUE, name = "display_context_count") %>%
  group_by(drug_normalized) %>%
  arrange(desc(display_context_count), nchar(display_candidate), display_candidate, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    drug_normalized,
    canonical_drug_display = display_candidate
  )

approved_annotated <- approved_raw %>%
  left_join(canonical_lookup, by = "drug_normalized") %>%
  mutate(
    canonical_drug_display = coalesce(canonical_drug_display, drug_display, reference_display, drug_normalized),
    common_substance_excluded_yes_no = ifelse(drug_normalized %in% common_substance_terms, "yes", "no"),
    exclusion_reason = ifelse(
      common_substance_excluded_yes_no == "yes",
      "Common ion/nutrient/metabolite/excipient-like term captured as approved active ingredient; retained in audit only.",
      ""
    )
  )

approved_exclusion_audit <- approved_annotated %>%
  filter(common_substance_excluded_yes_no == "yes") %>%
  group_by(drug_normalized, canonical_drug_display, exclusion_reason) %>%
  summarise(
    context_count = n(),
    pmid_count = n_distinct(pmid),
    death_modes = paste(sort(unique(as.character(death_mode))), collapse = "; "),
    top_mentions = paste(head(names(sort(table(drug_display), decreasing = TRUE)), 10), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(context_count), drug_normalized)

approved_refined <- approved_annotated %>%
  filter(common_substance_excluded_yes_no == "no") %>%
  distinct(pmid, tumor_family, disease_term, death_mode, drug_normalized, .keep_all = TRUE)

write_csv_safe(
  approved_exclusion_audit,
  file.path(ref_dir, "approved_clinical_drug_common_substance_exclusion_audit.csv")
)

write_csv_safe(
  approved_refined,
  file.path(ref_dir, "pubmed_pubtator_approved_clinical_drug_refined_contexts.csv")
)

relations_path <- file.path(ROOT, "09_sentence_level_relation_audit", "drug_death_relation_sentences.csv")
relations <- if (file.exists(relations_path)) {
  read_csv_required(relations_path) %>%
    mutate(
      pmid = as.character(pmid),
      death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", as.character(death_mode)),
      drug_normalized = normalize_drug_name(drug_normalized),
      relation_key = paste(pmid, death_mode, drug_normalized, sep = "||")
    )
} else {
  tibble(
    pmid = character(), death_mode = character(), drug_normalized = character(),
    relation_category = character(), confidence_high_medium_low = character(),
    relation_key = character()
  )
}

approved_keys <- approved_refined %>%
  mutate(relation_key = paste(pmid, as.character(death_mode), drug_normalized, sep = "||")) %>%
  select(relation_key, canonical_drug_display) %>%
  distinct()

relation_summary <- relations %>%
  inner_join(approved_keys, by = "relation_key") %>%
  mutate(
    induction_sensitization_flag = relation_category %in% c("induction_or_activation", "sensitization") &
      confidence_high_medium_low %in% c("high", "medium")
  ) %>%
  group_by(death_mode, drug_normalized) %>%
  summarise(
    relation_candidate_pmids = n_distinct(pmid),
    induction_sensitization_pmids = n_distinct(pmid[induction_sensitization_flag]),
    high_confidence_pmids = n_distinct(pmid[confidence_high_medium_low == "high"]),
    .groups = "drop"
  )

approved_by_mode <- approved_refined %>%
  group_by(death_mode, drug_normalized, canonical_drug_display) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    captured_context_count = n(),
    tumor_family_breadth = n_distinct(tumor_family),
    disease_count = n_distinct(disease_term),
    .groups = "drop"
  ) %>%
  left_join(relation_summary, by = c("death_mode", "drug_normalized")) %>%
  mutate(
    across(c(relation_candidate_pmids, induction_sensitization_pmids, high_confidence_pmids), ~coalesce(.x, 0L)),
    death_mode = factor(death_mode, levels = death_order)
  )

write_csv_safe(
  approved_by_mode,
  file.path(source_dir, "Figure_oncology_drug_death_triads_by_death_mode_approved_drugs_only_refined_source.csv")
)

top_by_mode <- approved_by_mode %>%
  group_by(death_mode) %>%
  arrange(desc(pmid_count), desc(induction_sensitization_pmids), desc(high_confidence_pmids), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  mutate(
    rank_within_death_mode = row_number(),
    entity_label_unique = paste(death_mode, sprintf("%02d", rank_within_death_mode), canonical_drug_display, sep = "___"),
    death_mode_plot = recode(as.character(death_mode), "Immunogenic cell death" = "ICD", .default = as.character(death_mode)),
    death_mode_plot = factor(death_mode_plot, levels = recode(death_order, "Immunogenic cell death" = "ICD", .default = death_order)),
    entity_label_unique = factor(entity_label_unique, levels = rev(entity_label_unique))
  ) %>%
  ungroup()

p_triads <- ggplot(top_by_mode, aes(x = pmid_count, y = entity_label_unique)) +
  geom_col(aes(fill = death_mode), width = 0.72, alpha = 0.9, show.legend = FALSE) +
  geom_point(aes(x = induction_sensitization_pmids), shape = 21, size = 1.7, stroke = 0.25, colour = "grey15", fill = "white") +
  facet_wrap(~ death_mode_plot, scales = "free", ncol = 4) +
  scale_y_discrete(labels = function(x) str_replace(x, "^.*___[0-9]+___", "")) +
  scale_fill_manual(values = death_mode_palette, drop = FALSE) +
  labs(
    title = "Approved clinical drugs captured by PubMed/PubTator across oncology death modes",
    subtitle = "Case/synonym-collapsed approved-drug layer; white dot = induction/sensitization sentence-candidate PMIDs.",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 5.7) +
  theme(strip.text = element_text(size = 5.8))

save_pub_r(
  p_triads,
  file.path(fig_dir, "Figure_oncology_drug_death_triads_by_death_mode_approved_drugs_only_refined"),
  width_mm = 210,
  height_mm = 235
)

multi <- approved_by_mode %>%
  group_by(drug_normalized, canonical_drug_display) %>%
  summarise(
    death_mode_breadth_induction_sensitization = sum(induction_sensitization_pmids > 0),
    death_mode_breadth_relation = sum(high_confidence_pmids > 0),
    death_mode_breadth_pmid = n_distinct(death_mode[pmid_count > 0]),
    total_induction_sensitization_count = sum(induction_sensitization_pmids, na.rm = TRUE),
    total_high_confidence_relation_count = sum(high_confidence_pmids, na.rm = TRUE),
    total_pmid_count = sum(pmid_count, na.rm = TRUE),
    top_death_modes = paste(as.character(death_mode[order(-pmid_count, -induction_sensitization_pmids)][1:min(5, n())]), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    multi_death_signal_score = zscore(death_mode_breadth_pmid) +
      zscore(log1p(total_pmid_count)) +
      zscore(death_mode_breadth_induction_sensitization) +
      zscore(log1p(total_induction_sensitization_count))
  ) %>%
  arrange(desc(death_mode_breadth_pmid), desc(total_pmid_count), desc(total_induction_sensitization_count))

write_csv_safe(
  multi,
  file.path(source_dir, "Figure_same_drug_multi_death_induction_sensitization_approved_drugs_only_refined_source.csv")
)

display_entities <- multi %>%
  filter(total_pmid_count > 0) %>%
  slice_head(n = 25) %>%
  mutate(canonical_drug_display = factor(canonical_drug_display, levels = rev(canonical_drug_display)))

p_same_a <- ggplot(display_entities, aes(x = total_pmid_count, y = canonical_drug_display)) +
  geom_col(aes(fill = death_mode_breadth_pmid), width = 0.7) +
  geom_point(
    aes(
      x = death_mode_breadth_induction_sensitization * max(total_pmid_count, na.rm = TRUE) /
        max(death_mode_breadth_induction_sensitization, na.rm = TRUE)
    ),
    shape = 21, fill = "white", colour = "grey20", stroke = 0.25, size = 2.1
  ) +
  scale_fill_viridis_c(option = "C", begin = 0.18, end = 0.88, breaks = 1:length(death_order), name = "Death-mode breadth") +
  labs(
    title = "A  Same approved drug across multiple death modes",
    x = "Total unique PMIDs; white dot = scaled induction/sensitization breadth",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.2) +
  theme(legend.position = "bottom")

p_same_b <- approved_by_mode %>%
  filter(drug_normalized %in% display_entities$drug_normalized) %>%
  mutate(
    canonical_drug_display = factor(canonical_drug_display, levels = rev(as.character(display_entities$canonical_drug_display))),
    death_mode = factor(death_mode, levels = death_order)
  ) %>%
  ggplot(aes(x = death_mode, y = canonical_drug_display, fill = log1p(pmid_count))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(size = induction_sensitization_pmids), shape = 21, fill = NA, colour = "grey25", stroke = 0.2) +
  scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p PMID") +
  scale_size_continuous(range = c(0.7, 4), name = "Induction/sensitization PMIDs") +
  labs(title = "B  Approved drug x death-mode support profile", x = NULL, y = NULL) +
  theme_oncology_atlas(base_size = 6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_same <- p_same_a / p_same_b +
  plot_layout(heights = c(0.95, 1.15)) +
  plot_annotation(
    title = "Approved clinical drugs with multi-death oncology literature signals",
    subtitle = "Approved-drug layer only; common ion/nutrient/metabolite-like entities are excluded from plotted drug signals.",
    theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
  )

save_pub_r(
  p_same,
  file.path(fig_dir, "Figure_same_drug_multi_death_induction_sensitization_approved_drugs_only_refined"),
  width_mm = 210,
  height_mm = 180
)

top_matrix_drugs <- multi %>%
  arrange(desc(death_mode_breadth_pmid), desc(total_pmid_count), desc(total_induction_sensitization_count)) %>%
  slice_head(n = 30) %>%
  mutate(canonical_drug_display = factor(canonical_drug_display, levels = rev(canonical_drug_display))) %>%
  select(drug_normalized, canonical_drug_display)

drug_mode_matrix <- top_matrix_drugs %>%
  left_join(approved_by_mode, by = c("drug_normalized", "canonical_drug_display")) %>%
  mutate(
    canonical_drug_display = factor(as.character(canonical_drug_display), levels = levels(top_matrix_drugs$canonical_drug_display)),
    death_mode = factor(death_mode, levels = death_order)
  )

death_mode_totals <- approved_by_mode %>%
  group_by(death_mode) %>%
  summarise(
    pmid_count = n_distinct(drug_normalized[pmin(pmid_count, 1) > 0]),
    total_pmid_support = sum(pmid_count, na.rm = TRUE),
    induction_sensitization_pmids = sum(induction_sensitization_pmids, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(total_pmid_support) %>%
  mutate(death_mode_ranked = factor(as.character(death_mode), levels = as.character(death_mode)))

write_csv_safe(
  drug_mode_matrix,
  file.path(source_dir, "Figure_therapy_class_death_mode_landscape_approved_drugs_only_refined_source.csv")
)

p_class_a <- ggplot(drug_mode_matrix, aes(x = death_mode, y = canonical_drug_display, fill = log1p(pmid_count))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(size = induction_sensitization_pmids), shape = 21, fill = NA, colour = "grey25", stroke = 0.2) +
  scale_fill_gradient(low = "grey96", high = "#B44D3A", name = "log1p PMID") +
  scale_size_continuous(range = c(0.4, 3.5), name = "Induction/sensitization PMIDs") +
  labs(title = "A  Top approved clinical drugs x death modes", x = NULL, y = NULL) +
  theme_oncology_atlas(base_size = 6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_class_b <- ggplot(death_mode_totals, aes(x = total_pmid_support, y = death_mode_ranked, fill = death_mode)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = comma(total_pmid_support)), hjust = -0.08, size = 2.1, family = "Arial") +
  scale_fill_manual(values = death_mode_palette, drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(title = "B  Approved-drug PMID support by death mode", x = "Summed drug-mode PMID support", y = NULL) +
  theme_oncology_atlas(base_size = 6.2)

p_class <- p_class_a | p_class_b +
  plot_layout(widths = c(1.35, 0.85)) +
  plot_annotation(
    title = "Approved clinical drug death-mode landscape",
    subtitle = "Canonicalized approved-drug names are shown after excluding common substance-like active ingredients from main figures.",
    theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
  )

save_pub_r(
  p_class,
  file.path(fig_dir, "Figure_therapy_class_death_mode_landscape_approved_drugs_only_refined"),
  width_mm = 210,
  height_mm = 140
)

approved_distribution <- approved_refined %>%
  count(canonical_drug_display, drug_normalized, sort = TRUE, name = "context_count") %>%
  summarise(
    top_examples = paste(head(paste0(canonical_drug_display, " (", comma(context_count), ")"), 15), collapse = "; "),
    .groups = "drop"
  ) %>%
  pull(top_examples)

retained_examples <- c("sorafenib", "auranofin", "cisplatin", "oxaliplatin", "doxorubicin", "paclitaxel", "temozolomide")
retained_check <- approved_refined %>%
  filter(drug_normalized %in% retained_examples) %>%
  count(drug_normalized, canonical_drug_display, sort = TRUE, name = "contexts") %>%
  mutate(line = paste0("- ", canonical_drug_display, " [", drug_normalized, "]: ", comma(contexts), " contexts")) %>%
  pull(line)

if (length(retained_check) == 0) retained_check <- "- No retained example drugs found in refined table."

summary_lines <- c(
  "# Approved Clinical Drug Refined Mapping Summary",
  "",
  "## Scope",
  "This refined layer keeps only PubMed/PubTator-captured entities mapped to the approved clinical drug reference class, then collapses case/synonym variants by normalized drug key.",
  "",
  "Common ion, nutrient, metabolite, solvent, excipient, or broad physiology terms are excluded from plotted drug signals and retained in an audit table. This prevents terms such as iron, copper, and glucose from dominating oncology therapy figures despite being FDA/Drugs@FDA active ingredients or products in some contexts.",
  "",
  "## Counts",
  paste0("- Approved clinical drug raw contexts before common-substance exclusion: ", comma(nrow(approved_raw))),
  paste0("- Approved clinical drug normalized entities before common-substance exclusion: ", comma(n_distinct(approved_raw$drug_normalized))),
  paste0("- Excluded common-substance contexts: ", comma(nrow(approved_annotated %>% filter(common_substance_excluded_yes_no == "yes")))),
  paste0("- Excluded common-substance normalized entities: ", comma(n_distinct((approved_annotated %>% filter(common_substance_excluded_yes_no == "yes"))$drug_normalized))),
  paste0("- Refined approved clinical drug contexts used in figures: ", comma(nrow(approved_refined))),
  paste0("- Refined approved clinical drug normalized entities used in figures: ", comma(n_distinct(approved_refined$drug_normalized))),
  paste0("- Refined approved clinical drug PMIDs: ", comma(n_distinct(approved_refined$pmid))),
  "",
  "## Top retained approved-drug examples",
  approved_distribution,
  "",
  "## Sanity-check retained examples",
  retained_check,
  "",
  "## Main common-substance exclusions",
  approved_exclusion_audit %>%
    slice_head(n = 20) %>%
    transmute(line = paste0("- ", canonical_drug_display, " [", drug_normalized, "]: ", comma(context_count), " contexts; ", comma(pmid_count), " PMIDs")) %>%
    pull(line),
  "",
  "## Interpretation boundary",
  "These outputs are approved clinical drug literature-signal figures, not evidence that a drug induces or treats a death mode. Sentence-candidate dots are relation candidates and require manual review for mechanism-level claims."
)

write_text_safe(summary_lines, file.path(qc_dir, "approved_clinical_drug_refined_mapping_summary.md"))

message("Generated refined approved clinical drug-only figures.")
