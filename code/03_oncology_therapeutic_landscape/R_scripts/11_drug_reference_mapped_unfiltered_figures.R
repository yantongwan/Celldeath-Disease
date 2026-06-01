script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

fig_dir <- file.path(ROOT, "10_publication_grade_figures")
source_dir <- file.path(fig_dir, "source_data")
ref_dir <- file.path(ROOT, "05_drug_therapy_extraction/drug_reference")
dir_create(fig_dir)
dir_create(source_dir)

normalize_drug_name <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- str_replace_all(x, "<[^>]+>", " ")
  x <- str_to_lower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_squish(x)
}

strip_salt_suffix <- function(x) {
  salt_terms <- c(
    "hydrochloride", "dihydrochloride", "hydrobromide", "bromide", "chloride",
    "sodium", "potassium", "calcium", "magnesium", "disodium", "monosodium",
    "phosphate", "diphosphate", "sulfate", "sulphate", "mesylate", "tosylate",
    "besylate", "acetate", "maleate", "tartrate", "citrate", "fumarate",
    "succinate", "lactate", "nitrate", "oxide", "hydrate", "monohydrate"
  )
  pat <- paste0("\\s+(", paste(salt_terms, collapse = "|"), ")(\\s+(", paste(salt_terms, collapse = "|"), "))*$")
  stripped <- str_squish(str_replace(x, pat, ""))
  ifelse(stripped == "", x, stripped)
}

make_ref_row <- function(term, display, reference_class, source, priority) {
  tibble(
    reference_term = term,
    reference_display = display,
    term_normalized = normalize_drug_name(term),
    term_normalized_salt_stripped = strip_salt_suffix(normalize_drug_name(term)),
    reference_class = reference_class,
    reference_source = source,
    priority = priority
  )
}

read_fda_approved <- function(path) {
  if (!file.exists(path)) return(tibble())
  read_csv(
    path,
    col_names = c("drugcentral_id", "drug_name"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(!is.na(drug_name), drug_name != "") %>%
    transmute(
      reference_term = drug_name,
      reference_display = drug_name,
      reference_class = "approved_clinical_drug",
      reference_source = "DrugCentral_FDA_Approved",
      priority = 1
    )
}

split_ingredients <- function(x) {
  x <- ifelse(is.na(x), "", x)
  parts <- unlist(str_split(x, "\\s*;\\s*|\\s+AND\\s+|\\s*/\\s*"))
  parts <- str_squish(parts)
  parts[parts != ""]
}

read_drugs_at_fda <- function(product_path) {
  if (!file.exists(product_path)) return(tibble())
  products <- read_tsv(
    product_path,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  active_terms <- products %>%
    select(DrugName, ActiveIngredient) %>%
    mutate(row_id = row_number()) %>%
    tidyr::separate_longer_delim(ActiveIngredient, delim = ";") %>%
    mutate(ActiveIngredient = str_squish(ActiveIngredient)) %>%
    filter(!is.na(ActiveIngredient), ActiveIngredient != "") %>%
    transmute(
      reference_term = ActiveIngredient,
      reference_display = str_to_title(ActiveIngredient),
      reference_class = "approved_clinical_drug",
      reference_source = "DrugsAtFDA_ActiveIngredient",
      priority = 2
    )
  product_terms <- products %>%
    filter(!is.na(DrugName), DrugName != "") %>%
    transmute(
      reference_term = DrugName,
      reference_display = str_to_title(DrugName),
      reference_class = "approved_product_or_brand",
      reference_source = "DrugsAtFDA_DrugName",
      priority = 3
    )
  bind_rows(active_terms, product_terms)
}

read_coconut <- function(zip_path) {
  if (!file.exists(zip_path)) return(tibble())
  csv_name <- "coconut_csv_lite-05-2026.csv"
  read_csv(
    unz(zip_path, csv_name),
    col_select = any_of(c("identifier", "name")),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(!is.na(name), name != "") %>%
    transmute(
      reference_term = name,
      reference_display = name,
      reference_class = "natural_product_candidate",
      reference_source = "COCONUT_2026_05_name",
      priority = 5
    )
}

read_local_lexicon <- function(path) {
  if (!file.exists(path)) return(tibble())
  read_csv_required(path) %>%
    filter(!is.na(drug_normalized), drug_normalized != "") %>%
    mutate(
      reference_class = case_when(
        is_specific_drug_yes_no == "yes" ~ "local_curated_specific_drug",
        is_generic_therapy_term_yes_no == "yes" | is_specific_drug_yes_no == "no" ~ "generic_therapy_term",
        TRUE ~ "local_curated_other"
      ),
      priority = case_when(
        reference_class == "local_curated_specific_drug" ~ 4,
        reference_class == "generic_therapy_term" ~ 6,
        TRUE ~ 7
      )
    ) %>%
    transmute(
      reference_term = drug_display,
      reference_display = drug_display,
      reference_class,
      reference_source = "local_curated_oncology_therapy_lexicon",
      priority
    )
}

reference_universe <- bind_rows(
  read_fda_approved(file.path(ref_dir, "FDA_Approved.csv")),
  read_drugs_at_fda(file.path(ref_dir, "datdaf20260519", "Products.txt")),
  read_coconut(file.path(ref_dir, "coconut_csv_lite-05-2026.zip")),
  read_local_lexicon(file.path(ROOT, "05_drug_therapy_extraction/drug_lexicon_normalized.csv"))
) %>%
  mutate(
    reference_term = str_squish(reference_term),
    reference_display = str_squish(reference_display),
    term_normalized = normalize_drug_name(reference_term),
    term_normalized_salt_stripped = strip_salt_suffix(term_normalized)
  ) %>%
  filter(term_normalized != "", nchar(term_normalized) > 1) %>%
  distinct(term_normalized, reference_class, reference_source, .keep_all = TRUE)

reference_lookup <- reference_universe %>%
  select(term_normalized, reference_display, reference_class, reference_source, priority) %>%
  bind_rows(
    reference_universe %>%
      transmute(
        term_normalized = term_normalized_salt_stripped,
        reference_display,
        reference_class,
        reference_source = paste0(reference_source, "_salt_stripped"),
        priority = priority + 0.1
      )
  ) %>%
  filter(term_normalized != "", nchar(term_normalized) > 1) %>%
  arrange(priority) %>%
  distinct(term_normalized, .keep_all = TRUE)

write_csv_safe(reference_universe, file.path(ref_dir, "drug_reference_universe_from_downloaded_sources.csv"))
write_csv_safe(reference_lookup, file.path(ref_dir, "drug_reference_lookup_terms_from_downloaded_sources.csv"))

death_order <- SELECTED_DEATH_MODES

captured <- read_csv_required(file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions.csv")) %>%
  mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
    death_mode = factor(death_mode, levels = death_order),
    mention_source = coalesce(mention_source, ""),
    candidate_normalized = normalize_drug_name(coalesce(drug_normalized, drug_display)),
    local_mapping_status = case_when(
      is_specific_drug_yes_no == "yes" ~ "local specific",
      is_generic_therapy_term_yes_no == "yes" | is_specific_drug_yes_no == "no" ~ "generic therapy term",
      TRUE ~ "needs review"
    )
  ) %>%
  filter(!is.na(death_mode), str_detect(mention_source, "PubMed_MeSH_ChemicalList|PubTator3_chemical")) %>%
  left_join(reference_lookup, by = c("candidate_normalized" = "term_normalized")) %>%
  mutate(
    reference_class = coalesce(reference_class, ifelse(local_mapping_status == "generic therapy term", "generic_therapy_term", "unmapped_pubmed_pubtator_entity")),
    reference_source = coalesce(reference_source, ifelse(local_mapping_status == "generic therapy term", "local_generic_therapy_term", "no_downloaded_reference_match")),
    reference_display = coalesce(reference_display, drug_display, drug_normalized),
    reference_class_label = recode(reference_class,
      approved_clinical_drug = "Approved clinical drug",
      approved_product_or_brand = "FDA product/brand",
      natural_product_candidate = "Natural product candidate",
      local_curated_specific_drug = "Local curated specific drug",
      generic_therapy_term = "Generic therapy term",
      unmapped_pubmed_pubtator_entity = "Unmapped PubMed/PubTator entity",
      .default = reference_class
    )
  ) %>%
  distinct(pmid, tumor_family, disease_term, death_mode, drug_normalized, .keep_all = TRUE)

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

write_csv_safe(captured, file.path(ref_dir, "pubmed_pubtator_candidates_reference_mapped.csv"))

relations <- read_csv_required(file.path(ROOT, "09_sentence_level_relation_audit/drug_death_relation_sentences.csv")) %>%
  mutate(
    pmid = as.character(pmid),
    death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
    relation_key = paste(pmid, death_mode, drug_normalized, sep = "||")
  )
captured_keys <- captured %>%
  mutate(relation_key = paste(pmid, death_mode, drug_normalized, sep = "||")) %>%
  select(relation_key, reference_class_label, reference_class, reference_source, reference_display) %>%
  distinct()
relations_mapped <- relations %>%
  inner_join(captured_keys, by = "relation_key")

relation_summary <- relations_mapped %>%
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

candidate_by_mode <- captured %>%
  group_by(death_mode, drug_normalized, drug_display, reference_class_label, reference_source) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    candidate_context_count = n(),
    tumor_family_breadth = n_distinct(tumor_family),
    disease_count = n_distinct(disease_term),
    .groups = "drop"
  ) %>%
  left_join(relation_summary, by = c("death_mode", "drug_normalized")) %>%
  mutate(across(c(relation_candidate_pmids, induction_sensitization_pmids, high_confidence_pmids), ~coalesce(.x, 0L)))

write_csv_safe(candidate_by_mode, file.path(source_dir, "Figure_oncology_drug_death_triads_by_death_mode_unfiltered_reference_mapped_source.csv"))

top_by_mode <- candidate_by_mode %>%
  group_by(death_mode) %>%
  arrange(desc(pmid_count), desc(induction_sensitization_pmids), desc(high_confidence_pmids), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  mutate(
    rank_within_death_mode = row_number(),
    entity_label = ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display),
    entity_label_unique = paste(death_mode, sprintf("%02d", rank_within_death_mode), entity_label, sep = "___"),
    death_mode_plot = recode(as.character(death_mode), "Immunogenic cell death" = "ICD", .default = as.character(death_mode)),
    death_mode_plot = factor(death_mode_plot, levels = recode(death_order, "Immunogenic cell death" = "ICD", .default = death_order)),
    reference_class_label = factor(reference_class_label, levels = class_order),
    entity_label_unique = factor(entity_label_unique, levels = rev(entity_label_unique))
  ) %>%
  ungroup()

p_triads <- ggplot(top_by_mode, aes(x = pmid_count, y = entity_label_unique, fill = reference_class_label)) +
  geom_col(width = 0.72, alpha = 0.88) +
  geom_point(aes(x = induction_sensitization_pmids), shape = 21, size = 1.7, stroke = 0.25, colour = "grey15", fill = "white") +
  facet_wrap(~ death_mode_plot, scales = "free", ncol = 4) +
  scale_y_discrete(labels = function(x) str_replace(x, "^.*___[0-9]+___", "")) +
  scale_fill_manual(values = class_palette, drop = FALSE, name = "Reference class") +
  labs(
    title = "Top PubMed/PubTator-captured candidates by death mode",
    subtitle = "Unfiltered captured layer; white dot = induction/sensitization sentence-candidate PMIDs.",
    x = "Unique PMIDs",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 5.7) +
  theme(legend.position = "bottom", strip.text = element_text(size = 5.8))

save_pub_r(
  p_triads,
  file.path(fig_dir, "Figure_oncology_drug_death_triads_by_death_mode_unfiltered_reference_mapped"),
  width_mm = 210,
  height_mm = 235
)

drug_mode <- candidate_by_mode
multi <- drug_mode %>%
  group_by(drug_normalized, drug_display, reference_class_label) %>%
  summarise(
    death_mode_breadth_induction_sensitization = sum(induction_sensitization_pmids > 0),
    death_mode_breadth_relation = sum(high_confidence_pmids > 0),
    death_mode_breadth_pmid = n_distinct(death_mode[pmid_count > 0]),
    total_induction_sensitization_count = sum(induction_sensitization_pmids, na.rm = TRUE),
    total_high_confidence_relation_count = sum(high_confidence_pmids, na.rm = TRUE),
    total_pmid_count = sum(pmid_count, na.rm = TRUE),
    top_death_modes = paste(as.character(death_mode[order(-induction_sensitization_pmids, -pmid_count)][1:min(5, n())]), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    entity_label = ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display),
    multi_death_induction_score = zscore(death_mode_breadth_induction_sensitization) +
      zscore(log1p(total_induction_sensitization_count)) +
      zscore(death_mode_breadth_relation) +
      zscore(log1p(total_pmid_count))
  ) %>%
  arrange(desc(death_mode_breadth_induction_sensitization), desc(total_induction_sensitization_count), desc(total_pmid_count))

write_csv_safe(multi, file.path(source_dir, "Figure_same_drug_multi_death_induction_sensitization_unfiltered_reference_mapped_source.csv"))

display_entities <- multi %>%
  filter(death_mode_breadth_induction_sensitization >= 1, total_induction_sensitization_count > 0) %>%
  slice_head(n = 25) %>%
  mutate(
    entity_label = factor(entity_label, levels = rev(entity_label)),
    reference_class_label = factor(reference_class_label, levels = class_order)
  )

p_same_a <- ggplot(display_entities, aes(x = total_induction_sensitization_count, y = entity_label, fill = reference_class_label)) +
  geom_col(width = 0.7) +
  geom_point(
    aes(x = death_mode_breadth_induction_sensitization * max(total_induction_sensitization_count, na.rm = TRUE) /
          max(death_mode_breadth_induction_sensitization, na.rm = TRUE)),
    shape = 21, fill = "white", colour = "grey20", stroke = 0.25, size = 2.1
  ) +
  scale_fill_manual(values = class_palette, drop = FALSE, name = "Reference class") +
  labs(
    title = "A  Multi-death induction/sensitization candidate breadth",
    x = "Induction/sensitization candidate PMIDs; white dot = scaled death-mode breadth",
    y = NULL
  ) +
  theme_oncology_atlas(base_size = 6.2) +
  theme(legend.position = "bottom")

p_same_b <- drug_mode %>%
  filter(drug_normalized %in% display_entities$drug_normalized) %>%
  mutate(
    entity_label = ifelse(is.na(drug_display) | drug_display == "", drug_normalized, drug_display),
    entity_label = factor(entity_label, levels = rev(as.character(display_entities$entity_label))),
    death_mode = factor(death_mode, levels = death_order)
  ) %>%
  ggplot(aes(x = death_mode, y = entity_label, fill = log1p(induction_sensitization_pmids))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(size = pmid_count), shape = 21, fill = NA, colour = "grey25", stroke = 0.2) +
  scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p induction/sensitization") +
  scale_size_continuous(range = c(0.7, 4), name = "PMID count") +
  labs(title = "B  Candidate x death-mode support profile", x = NULL, y = NULL) +
  theme_oncology_atlas(base_size = 6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_same <- p_same_a / p_same_b +
  plot_layout(heights = c(0.95, 1.15)) +
  plot_annotation(
    title = "Unfiltered same-candidate multi-death induction/sensitization signals",
    subtitle = "This is a PubMed/PubTator candidate view. Non-drug chemicals and protein entities remain visible for downstream screening.",
    theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
  )

save_pub_r(
  p_same,
  file.path(fig_dir, "Figure_same_drug_multi_death_induction_sensitization_unfiltered_reference_mapped"),
  width_mm = 210,
  height_mm = 180
)

class_mode <- captured %>%
  group_by(reference_class_label, death_mode) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    candidate_context_count = n(),
    entity_count = n_distinct(drug_normalized),
    tumor_family_breadth = n_distinct(tumor_family),
    .groups = "drop"
  )

class_totals <- class_mode %>% group_by(reference_class_label) %>% summarise(class_total = sum(pmid_count), .groups = "drop")
death_totals <- class_mode %>% group_by(death_mode) %>% summarise(death_total = sum(pmid_count), .groups = "drop")
all_total <- sum(class_mode$pmid_count)

class_mode <- class_mode %>%
  left_join(class_totals, by = "reference_class_label") %>%
  left_join(death_totals, by = "death_mode") %>%
  mutate(
    reference_class_label = factor(reference_class_label, levels = class_order),
    death_mode = factor(death_mode, levels = death_order),
    class_death_selectivity = log2((pmid_count / class_total + 1e-6) / (death_total / all_total + 1e-6)),
    clipped_selectivity = pmax(pmin(class_death_selectivity, 3), -3)
  ) %>%
  complete(
    reference_class_label = factor(class_order, levels = class_order),
    death_mode = factor(death_order, levels = death_order),
    fill = list(pmid_count = 0, candidate_context_count = 0, entity_count = 0, tumor_family_breadth = 0, clipped_selectivity = NA_real_)
  )

write_csv_safe(class_mode, file.path(source_dir, "Figure_therapy_class_death_mode_landscape_unfiltered_reference_mapped_source.csv"))

p_class_a <- ggplot(class_mode, aes(x = death_mode, y = reference_class_label, fill = clipped_selectivity)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data = class_mode %>% filter(pmid_count > 0),
    aes(size = log1p(pmid_count)),
    shape = 21, colour = "grey20", fill = NA, stroke = 0.25
  ) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#B44D3A", midpoint = 0, limits = c(-3, 3), na.value = "grey92", name = "log2 selectivity") +
  scale_size_continuous(range = c(0.7, 4), name = "log1p PMID count") +
  labs(title = "A  Reference class x death-mode candidate landscape", x = NULL, y = NULL) +
  theme_oncology_atlas() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_class_b <- captured %>%
  count(reference_class_label, name = "captured_contexts") %>%
  mutate(reference_class_label = factor(reference_class_label, levels = rev(class_order))) %>%
  ggplot(aes(x = captured_contexts, y = reference_class_label, fill = reference_class_label)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = class_palette, drop = FALSE) +
  labs(title = "B  Captured contexts by reference class", x = "Captured PMID-entity contexts (log10)", y = NULL) +
  theme_oncology_atlas()

p_class <- p_class_a | p_class_b +
  plot_layout(widths = c(1.15, 0.85)) +
  plot_annotation(
    title = "Reference-mapped PubMed/PubTator candidate landscape",
    subtitle = "Unfiltered candidates are shown before downstream manual or database-based exclusion.",
    theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
  )

save_pub_r(
  p_class,
  file.path(fig_dir, "Figure_therapy_class_death_mode_landscape_unfiltered_reference_mapped"),
  width_mm = 210,
  height_mm = 120
)

reference_class_distribution <- captured %>%
  count(reference_class_label, sort = TRUE) %>%
  transmute(line = paste0("- ", reference_class_label, ": ", comma(n))) %>%
  pull(line)

summary_lines <- c(
  "# Downloaded Drug Reference Mapping Summary",
  "",
  paste0("- Reference terms loaded: ", comma(nrow(reference_lookup))),
  paste0("- PubMed/PubTator captured candidate contexts: ", comma(nrow(captured))),
  paste0("- Unique captured PMIDs: ", comma(n_distinct(captured$pmid))),
  paste0("- Unique captured entities: ", comma(n_distinct(captured$drug_normalized))),
  "",
  "## Reference class distribution",
  reference_class_distribution,
  "",
  "Interpretation boundary: these figures are unfiltered PubMed/PubTator captured candidate views mapped to the downloaded FDA/Drugs@FDA/COCONUT references. They should not be interpreted as validated therapy-death mechanisms."
)
write_text_safe(summary_lines, file.path(ROOT, "15_QC_reports/downloaded_drug_reference_mapping_summary.md"))

message("Generated unfiltered reference-mapped drug candidate figures.")
