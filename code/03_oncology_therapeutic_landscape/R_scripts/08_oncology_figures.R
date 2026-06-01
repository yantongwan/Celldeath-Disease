script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

gate_path <- file.path(ROOT, "15_QC_reports/oncology_analysis_gate_decision_table.csv")
if (!file.exists(gate_path)) stop("Gate table missing. Run 01_oncology_analysis_pipeline.R first.")
gates <- read_csv_required(gate_path)
fig_dir <- file.path(ROOT, "10_publication_grade_figures")
dir_create(fig_dir)

gate_passes <- function(module) {
  val <- gates %>% filter(.data$module == !!module) %>% pull(passes_gate_yes_no)
  length(val) > 0 && val[1] == "yes"
}

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

figure_log <- tibble(
  figure = character(), status = character(), reason = character(),
  pdf = character(), png = character(), tiff = character()
)
add_log <- function(figure, status, reason, base = "") {
  figure_log <<- bind_rows(figure_log, tibble(
    figure = figure,
    status = status,
    reason = reason,
    pdf = ifelse(base == "", "", paste0(base, ".pdf")),
    png = ifelse(base == "", "", paste0(base, ".png")),
    tiff = ifelse(base == "", "", paste0(base, ".tiff"))
  ))
}

make_hotness <- function() {
  hot <- read_csv_required(file.path(ROOT, "03_oncology_death_mode_hotness/oncology_death_mode_hotness.csv"))
  hot <- hot %>%
    mutate(death_mode = factor(death_mode, levels = rev(death_mode[order(heat_index)])))

  p1 <- ggplot(hot, aes(x = death_mode, y = heat_index, fill = as.character(death_mode))) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = death_mode_palette, guide = "none") +
    labs(title = "A  Oncology literature heat index", x = NULL, y = "Component-wise heat index") +
    theme_oncology_atlas()

  p2 <- ggplot(hot, aes(x = oncology_unique_pmids, y = recent_growth_3yr, colour = as.character(death_mode), size = tumor_family_breadth)) +
    geom_point(alpha = 0.86) +
    scale_x_log10(labels = comma) +
    scale_colour_manual(values = death_mode_palette, guide = "none") +
    scale_size_continuous(range = c(2, 6), name = "Tumor-family breadth") +
    labs(title = "B  Volume and recent growth", x = "Oncology unique PMIDs (log10)", y = "2023-2025 / 2020-2022") +
    theme_oncology_atlas()

  component <- hot %>%
    select(death_mode, oncology_unique_pmids, tumor_family_breadth, recent_growth_3yr, drug_associated_pmid_fraction) %>%
    mutate(death_mode = as.character(death_mode),
           oncology_unique_pmids = scales::rescale(log1p(oncology_unique_pmids)),
           tumor_family_breadth = scales::rescale(tumor_family_breadth),
           recent_growth_3yr = scales::rescale(recent_growth_3yr),
           drug_associated_pmid_fraction = scales::rescale(drug_associated_pmid_fraction)) %>%
    pivot_longer(-death_mode, names_to = "component", values_to = "scaled_value") %>%
    mutate(component = recode(component,
      oncology_unique_pmids = "PMID support",
      tumor_family_breadth = "Tumor breadth",
      recent_growth_3yr = "Recent growth",
      drug_associated_pmid_fraction = "Drug vocabulary"
    ))

  p3 <- ggplot(component, aes(x = factor(death_mode, levels = levels(hot$death_mode)), y = scaled_value, fill = component)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(title = "C  Heat-index components", x = NULL, y = "Scaled component") +
    theme_oncology_atlas()

  hot_label <- hot %>%
    mutate(label_short = recode(as.character(death_mode), "Immunogenic cell death" = "ICD", .default = as.character(death_mode)))
  p4 <- ggplot(hot_label, aes(x = fraction_of_death_mode_literature_in_oncology, y = drug_associated_pmid_fraction, colour = as.character(death_mode))) +
    geom_point(size = 3.0) +
    geom_text(aes(label = label_short), nudge_x = 0.018, hjust = 0, size = 2.05, show.legend = FALSE) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.9)) +
    scale_colour_manual(values = death_mode_palette, guide = "none") +
    labs(title = "D  Oncology fraction and drug vocabulary", x = "Fraction of death-mode records in oncology", y = "Drug-associated PMID fraction") +
    coord_cartesian(clip = "off") +
    theme_oncology_atlas() +
    theme(plot.margin = margin(5.5, 24, 5.5, 5.5))

  plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Oncology RCD death-mode hotness is descriptive and component-wise",
      subtitle = "Signals summarize literature attention, tumor breadth, recent growth, and drug-associated vocabulary; they are not target-quality rankings.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_oncology_death_mode_hotness")
  save_pub_r(plot, base, width_mm = 183, height_mm = 135)
  add_log("Figure_oncology_death_mode_hotness", "generated", "death_mode_hotness gate passed", base)
}

make_tumor_landscape <- function() {
  td <- read_csv_required(file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_selectivity.csv"))
  top_tumors <- td %>%
    group_by(tumor_family) %>%
    summarise(total = sum(pair_count), breadth = sum(pair_count >= 3), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = 16) %>%
    pull(tumor_family)
  td_plot <- td %>%
    filter(tumor_family %in% top_tumors) %>%
    mutate(
      tumor_family = factor(tumor_family, levels = rev(top_tumors)),
      death_mode = factor(death_mode, levels = names(death_mode_palette)[names(death_mode_palette) %in% unique(death_mode)]),
      support_label = ifelse(pair_count >= 3, "", "low")
    )
  p1 <- ggplot(td_plot, aes(x = death_mode, y = tumor_family, fill = pmax(pmin(log2_tumor_death_selectivity, 3), -3))) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_point(aes(size = pair_count), shape = 21, stroke = 0.25, colour = "grey20", fill = NA) +
    scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#B44D3A", midpoint = 0, limits = c(-3, 3), name = "log2 selectivity") +
    scale_size_continuous(range = c(0.8, 4.5), name = "Pair count") +
    labs(title = "A  Tumor-family x death-mode selectivity", x = NULL, y = NULL) +
    theme_oncology_atlas() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  breadth <- td %>%
    filter(tumor_family %in% top_tumors) %>%
    group_by(tumor_family) %>%
    summarise(death_mode_breadth = sum(pair_count >= 3), total_pair_count = sum(pair_count), .groups = "drop") %>%
    mutate(tumor_family = factor(tumor_family, levels = rev(top_tumors)))

  p2 <- ggplot(breadth, aes(x = death_mode_breadth, y = tumor_family)) +
    geom_col(width = 0.7, fill = "#727272") +
    labs(title = "B  Death-mode breadth by tumor family", x = "Death modes with Pair_Count >= 3", y = NULL) +
    theme_oncology_atlas()

  examples <- td %>%
    filter(tumor_family %in% c("liver cancer", "breast cancer", "lung cancer", "glioblastoma", "colorectal cancer", "melanoma", "renal cancer", "pancreatic cancer")) %>%
    group_by(tumor_family) %>%
    slice_max(log2_tumor_death_selectivity, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = paste0(death_mode, " (", round(log2_tumor_death_selectivity, 1), ")"),
           tumor_family = factor(tumor_family, levels = rev(unique(tumor_family))))

  p3 <- ggplot(examples, aes(x = log2_tumor_death_selectivity, y = reorder(label, log2_tumor_death_selectivity), colour = death_mode)) +
    geom_point(size = 2.2) +
    facet_wrap(~tumor_family, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = death_mode_palette, guide = "none") +
    labs(title = "C  Representative selective tumor-death profiles", x = "log2 selectivity", y = NULL) +
    theme_oncology_atlas(base_size = 6.6)

  plot <- (p1 | p2) / p3 +
    plot_layout(heights = c(1.45, 1)) +
    plot_annotation(
      title = "Tumor families show heterogeneous RCD-associated literature profiles",
      subtitle = "Selectivity is normalized to the oncology background and should not be read as proven tumor biology.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_tumor_specific_death_landscape")
  save_pub_r(plot, base, width_mm = 183, height_mm = 180)
  add_log("Figure_tumor_specific_death_landscape", "generated", "tumor_specificity gate passed", base)
}

make_same_drug <- function() {
  triads <- read_csv_required(file.path(ROOT, "06_drug_death_mode_triads/oncology_drug_death_triads.csv")) %>%
    mutate(
      induction_sensitization_count = induction_count + sensitization_count,
      death_mode = factor(death_mode, levels = names(death_mode_palette)[names(death_mode_palette) %in% unique(death_mode)])
    )

  drug_mode <- triads %>%
    group_by(drug_normalized, drug_display, therapy_class, death_mode) %>%
    summarise(
      pmid_count = sum(pmid_count, na.rm = TRUE),
      high_confidence_relation_count = sum(high_confidence_relation_count, na.rm = TRUE),
      induction_sensitization_count = sum(induction_sensitization_count, na.rm = TRUE),
      tumor_family_breadth = n_distinct(tumor_family),
      top_tumor_families = paste(head(unique(tumor_family[order(-pmid_count)]), 4), collapse = "; "),
      .groups = "drop"
    )

  multi <- drug_mode %>%
    group_by(drug_normalized, drug_display, therapy_class) %>%
    summarise(
      death_mode_breadth_induction_sensitization = sum(induction_sensitization_count > 0),
      death_mode_breadth_relation = sum(high_confidence_relation_count > 0),
      death_mode_breadth_pmid = n_distinct(death_mode[pmid_count > 0]),
      total_induction_sensitization_count = sum(induction_sensitization_count, na.rm = TRUE),
      total_high_confidence_relation_count = sum(high_confidence_relation_count, na.rm = TRUE),
      total_pmid_count = sum(pmid_count, na.rm = TRUE),
      top_death_modes = paste(as.character(death_mode[order(-induction_sensitization_count, -high_confidence_relation_count, -pmid_count)][1:min(5, n())]), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      multi_death_induction_score = zscore(death_mode_breadth_induction_sensitization) +
        zscore(log1p(total_induction_sensitization_count)) +
        zscore(death_mode_breadth_relation) +
        zscore(log1p(total_pmid_count))
    ) %>%
    arrange(desc(death_mode_breadth_induction_sensitization), desc(total_induction_sensitization_count), desc(total_pmid_count))

  write_csv_safe(multi, file.path(ROOT, "07_same_drug_multi_death_signals/drug_multi_death_induction_sensitization_table.csv"))

  display_drugs <- multi %>%
    filter(death_mode_breadth_induction_sensitization >= 2, total_induction_sensitization_count > 0) %>%
    slice_head(n = 20)
  if (nrow(display_drugs) < 3) {
    add_log("Figure_same_drug_multi_death_signals", "rejected", "fewer than three drugs with induction/sensitization support across at least two death modes")
    return(invisible(NULL))
  }

  top_drugs <- display_drugs$drug_normalized

  p1 <- display_drugs %>%
    mutate(drug_display = factor(drug_display, levels = rev(drug_display))) %>%
    ggplot(aes(x = total_induction_sensitization_count, y = drug_display, fill = therapy_class)) +
    geom_col(width = 0.7) +
    geom_point(aes(x = death_mode_breadth_induction_sensitization * max(total_induction_sensitization_count, na.rm = TRUE) / max(death_mode_breadth_induction_sensitization, na.rm = TRUE)),
               shape = 21, fill = "white", colour = "grey20", stroke = 0.25, size = 2.2) +
    scale_fill_manual(values = therapy_palette, labels = class_labels[names(therapy_palette)], name = "Therapy class") +
    labs(title = "A  Same-drug induction/sensitization breadth", x = "Induction/sensitization candidate PMIDs; white dot = scaled death-mode breadth", y = NULL) +
    theme_oncology_atlas()

  p2 <- drug_mode %>%
    filter(drug_normalized %in% top_drugs) %>%
    mutate(
      drug_display = factor(drug_display, levels = rev(display_drugs$drug_display)),
      death_mode = factor(death_mode, levels = names(death_mode_palette)[names(death_mode_palette) %in% as.character(unique(death_mode))])
    ) %>%
    ggplot(aes(x = death_mode, y = drug_display, fill = log1p(induction_sensitization_count))) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_point(aes(size = pmid_count), shape = 21, fill = NA, colour = "grey25", stroke = 0.2) +
    scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p induction/sensitization candidates") +
    scale_size_continuous(range = c(0.8, 4), name = "PMID count") +
    labs(title = "B  Drug x death-mode candidate support", x = NULL, y = NULL) +
    theme_oncology_atlas() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  plot <- p1 / p2 +
    plot_annotation(
      title = "Some specific drugs recur in induction/sensitization candidates across multiple RCD concepts",
      subtitle = "This is a literature-level sentence-rule signal; it does not prove that a drug induces multiple death pathways.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_same_drug_multi_death_signals")
  save_pub_r(plot, base, width_mm = 183, height_mm = 165)
  add_log("Figure_same_drug_multi_death_signals", "generated", "same-drug figure rebuilt around induction/sensitization support across death modes", base)
}

make_triads <- function() {
  all_mentions <- read_csv_required(file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions.csv"))
  all_relations <- read_csv_required(file.path(ROOT, "09_sentence_level_relation_audit/drug_death_relation_sentences.csv"))
  triads <- read_csv_required(file.path(ROOT, "06_drug_death_mode_triads/oncology_drug_death_triads.csv")) %>%
    mutate(induction_sensitization_count = induction_count + sensitization_count)

  relation_summary_all <- all_relations %>%
    mutate(induction_sensitization_flag = relation_category %in% c("induction_or_activation", "sensitization") &
             confidence_high_medium_low %in% c("high", "medium")) %>%
    group_by(death_mode, drug_normalized) %>%
    summarise(
      relation_candidate_pmids = n_distinct(pmid),
      induction_sensitization_pmids = n_distinct(pmid[induction_sensitization_flag]),
      high_confidence_pmids = n_distinct(pmid[confidence_high_medium_low == "high"]),
      .groups = "drop"
    )

  all_candidate_summary <- all_mentions %>%
    group_by(death_mode, drug_normalized, drug_display, therapy_class, is_specific_drug_yes_no, mapping_confidence) %>%
    summarise(
      pmid_count = n_distinct(pmid),
      tumor_family_breadth = n_distinct(tumor_family),
      .groups = "drop"
    ) %>%
    left_join(relation_summary_all, by = c("death_mode", "drug_normalized")) %>%
    mutate(across(c(relation_candidate_pmids, induction_sensitization_pmids, high_confidence_pmids), ~coalesce(.x, 0L))) %>%
    arrange(death_mode, desc(induction_sensitization_pmids), desc(high_confidence_pmids), desc(pmid_count))
  write_csv_safe(all_candidate_summary, file.path(ROOT, "05_drug_therapy_extraction/all_pubmed_pubtator_chemical_candidate_summary_by_death_mode.csv"))

  specific_scope_summary <- all_candidate_summary %>%
    count(is_specific_drug_yes_no, name = "death_mode_drug_candidate_count") %>%
    arrange(desc(death_mode_drug_candidate_count))
  write_csv_safe(specific_scope_summary, file.path(ROOT, "05_drug_therapy_extraction/drug_scope_specificity_summary.csv"))
  write_text_safe(c(
    "# Drug Scope Audit",
    "",
    paste0("- PubMed/PubTator/local merged context rows: ", nrow(all_mentions), "."),
    paste0("- Specific mapped drug context rows used for triad figures: ", sum(all_mentions$is_specific_drug_yes_no == "yes", na.rm = TRUE), "."),
    paste0("- Needs-review chemical context rows not promoted to specific-drug triads: ", sum(all_mentions$is_specific_drug_yes_no == "needs_review", na.rm = TRUE), "."),
    "- Current triad figures use specific mapped drugs only, not every PubTator chemical annotation.",
    "- Reason: PubTator chemical annotations include ions, metabolites, lipids, ROS, ATP, copper, iron, and other chemical entities that are not necessarily drugs.",
    "- The full all-candidate audit table is saved as `05_drug_therapy_extraction/all_pubmed_pubtator_chemical_candidate_summary_by_death_mode.csv`."
  ), file.path(ROOT, "15_QC_reports/drug_scope_audit_all_candidates_vs_specific_mapped.md"))

  top_by_mode <- triads %>%
    group_by(death_mode, drug_normalized, drug_display, therapy_class) %>%
    summarise(
      pmid_count = sum(pmid_count, na.rm = TRUE),
      high_confidence_relation_count = sum(high_confidence_relation_count, na.rm = TRUE),
      induction_sensitization_count = sum(induction_sensitization_count, na.rm = TRUE),
      tumor_family_breadth = n_distinct(tumor_family),
      top_tumor_families = paste(head(unique(tumor_family[order(-pmid_count)]), 4), collapse = "; "),
      strongest_grade = ifelse(any(triad_evidence_grade == "strong_literature_signal"), "strong_literature_signal",
        ifelse(any(triad_evidence_grade == "moderate_literature_signal"), "moderate_literature_signal", "exploratory")),
      .groups = "drop"
    ) %>%
    group_by(death_mode) %>%
    arrange(desc(induction_sensitization_count), desc(high_confidence_relation_count), desc(pmid_count), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    mutate(
      rank_within_death_mode = row_number(),
      drug_label = paste0(rank_within_death_mode, ". ", drug_display),
      drug_label_unique = paste(death_mode, sprintf("%02d", rank_within_death_mode), drug_display, sep = "___")
    ) %>%
    ungroup()
  write_csv_safe(top_by_mode, file.path(ROOT, "06_drug_death_mode_triads/top10_specific_drugs_by_death_mode.csv"))

  if (nrow(top_by_mode) < 10) {
    add_log("Figure_oncology_drug_death_triads", "rejected", "fewer than ten death-mode grouped specific-drug candidates")
    return(invisible(NULL))
  }

  top_by_mode <- top_by_mode %>%
    mutate(
      death_mode = factor(death_mode, levels = names(death_mode_palette)[names(death_mode_palette) %in% unique(death_mode)]),
      death_mode_plot = recode(as.character(death_mode), "Immunogenic cell death" = "ICD", .default = as.character(death_mode)),
      death_mode_plot = factor(death_mode_plot, levels = recode(names(death_mode_palette)[names(death_mode_palette) %in% as.character(unique(death_mode))], "Immunogenic cell death" = "ICD", .default = names(death_mode_palette)[names(death_mode_palette) %in% as.character(unique(death_mode))])),
      drug_label_unique = factor(drug_label_unique, levels = rev(drug_label_unique))
    )

  p1 <- ggplot(top_by_mode, aes(x = pmid_count, y = drug_label_unique, fill = therapy_class)) +
    geom_col(width = 0.72, alpha = 0.86) +
    geom_point(aes(x = induction_sensitization_count), shape = 21, size = 1.8, stroke = 0.25, colour = "grey15", fill = "white") +
    facet_wrap(~ death_mode_plot, scales = "free", ncol = 4) +
    scale_y_discrete(labels = function(x) str_replace(x, "^.*___[0-9]+___", "")) +
    scale_fill_manual(values = therapy_palette, labels = class_labels[names(therapy_palette)], name = "Therapy class") +
    labs(title = "A  Top 10 specific drugs by death mode", x = "PMID count; white dot = induction/sensitization candidate PMIDs", y = NULL) +
    theme_oncology_atlas(base_size = 5.8) +
    theme(legend.position = "bottom", strip.text = element_text(size = 5.8))

  p2 <- top_by_mode %>%
    group_by(therapy_class, death_mode) %>%
    summarise(
      induction_sensitization_count = sum(induction_sensitization_count, na.rm = TRUE),
      pmid_count = sum(pmid_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(therapy_class = class_labels[therapy_class],
           death_mode_plot = recode(as.character(death_mode), "Immunogenic cell death" = "ICD", .default = as.character(death_mode)),
           death_mode_plot = factor(death_mode_plot, levels = recode(names(death_mode_palette)[names(death_mode_palette) %in% as.character(unique(death_mode))], "Immunogenic cell death" = "ICD", .default = names(death_mode_palette)[names(death_mode_palette) %in% as.character(unique(death_mode))]))) %>%
    ggplot(aes(x = death_mode_plot, y = therapy_class, fill = log1p(induction_sensitization_count))) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_point(aes(size = pmid_count), shape = 21, fill = NA, colour = "grey25", stroke = 0.2) +
    scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p induction/sensitization candidates") +
    scale_size_continuous(range = c(0.7, 3.2), name = "PMID count") +
    labs(title = "B  Top-drug support by therapy class and death mode", x = NULL, y = NULL) +
    theme_oncology_atlas(base_size = 6.6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  plot <- p1 / p2 +
    plot_layout(heights = c(1.85, 0.75)) +
    plot_annotation(
      title = "Drug signals stratified by RCD concept",
      subtitle = "Specific mapped drugs only; white dots mark rule-based induction/sensitization candidates requiring manual validation.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_oncology_drug_death_triads")
  save_pub_r(plot, base, width_mm = 183, height_mm = 230)
  add_log("Figure_oncology_drug_death_triads", "generated", "triad figure rebuilt as top 10 specific drugs by death mode", base)
}

make_sentence_audit <- function() {
  rel <- read_csv_required(file.path(ROOT, "09_sentence_level_relation_audit/drug_death_relation_sentences.csv"))
  if (nrow(rel) < 500) {
    add_log("Figure_sentence_level_relation_audit", "rejected", "fewer than 500 sentence-level relation candidates")
    return(invisible(NULL))
  }
  rel <- rel %>%
    mutate(
      relation_label = recode(relation_category,
        "induction_or_activation" = "Induction/activation",
        "sensitization" = "Sensitization",
        "inhibition_or_protection" = "Inhibition/protection",
        "resistance_escape" = "Resistance/escape",
        "combination_therapy" = "Combination",
        "unclear_co_mention" = "Unclear co-mention"
      ),
      confidence_high_medium_low = factor(confidence_high_medium_low, levels = c("high", "medium", "low"))
    )
  p1 <- rel %>%
    count(relation_label, confidence_high_medium_low, name = "records") %>%
    mutate(relation_label = factor(relation_label, levels = rev(c("Induction/activation", "Sensitization", "Inhibition/protection", "Combination", "Resistance/escape", "Unclear co-mention")))) %>%
    ggplot(aes(x = records, y = relation_label, fill = confidence_high_medium_low)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c(high = "#B44D3A", medium = "#D9A066", low = "#BDBDBD"), name = "Confidence") +
    labs(title = "A  Rule-based sentence relation candidates", x = "Sentence-level candidates", y = NULL) +
    theme_oncology_atlas()

  p2 <- rel %>%
    filter(confidence_high_medium_low %in% c("high", "medium")) %>%
    count(death_mode, relation_label, name = "records") %>%
    mutate(death_mode = factor(death_mode, levels = names(death_mode_palette)[names(death_mode_palette) %in% unique(death_mode)]),
           relation_label = factor(relation_label, levels = c("Induction/activation", "Sensitization", "Inhibition/protection", "Combination", "Resistance/escape", "Unclear co-mention"))) %>%
    ggplot(aes(x = relation_label, y = death_mode, fill = log1p(records))) +
    geom_tile(colour = "white", linewidth = 0.28) +
    scale_fill_gradient(low = "grey95", high = "#B44D3A", name = "log1p records") +
    labs(title = "B  Death mode x relation cue structure", x = NULL, y = NULL) +
    theme_oncology_atlas() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  p3 <- rel %>%
    summarise(
      sentence_candidates = n(),
      high_confidence_induction_or_sensitization = sum(confidence_high_medium_low == "high" & relation_category %in% c("induction_or_activation", "sensitization")),
      manual_validation_sample = nrow(read_csv_required(file.path(ROOT, "09_sentence_level_relation_audit/drug_relation_manual_validation_sample.csv")))
    ) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
      sentence_candidates = "Sentence candidates",
      high_confidence_induction_or_sensitization = "High-confidence induction/sensitization",
      manual_validation_sample = "Manual-validation sample"
    )) %>%
    ggplot(aes(x = value, y = factor(metric, levels = rev(metric)))) +
    geom_col(width = 0.65, fill = "#727272") +
    labs(title = "C  Audit thresholds", x = "Records", y = NULL) +
    theme_oncology_atlas()

  plot <- (p1 | p2) / p3 +
    plot_layout(heights = c(1, 0.58)) +
    plot_annotation(
      title = "Sentence-level relation audit supports cautious therapy-linked RCD vocabulary",
      subtitle = "Candidates are rule-based and require manual validation before drug-induction claims.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_sentence_level_relation_audit")
  save_pub_r(plot, base, width_mm = 183, height_mm = 150)
  add_log("Figure_sentence_level_relation_audit", "generated", "sentence_relation_audit gate passed for supplementary audit figure", base)
}

make_therapy_class <- function() {
  mat <- read_csv_required(file.path(ROOT, "08_therapy_class_landscape/therapy_class_death_mode_matrix.csv"))
  if (nrow(mat) == 0 || n_distinct(mat$therapy_class) < 3 || n_distinct(mat$death_mode) < 4) {
    add_log("Figure_therapy_class_death_mode_landscape", "rejected", "insufficient therapy-class/death-mode relation support")
    return(invisible(NULL))
  }
  death_mode_order <- names(death_mode_palette)
  death_mode_order <- intersect(SELECTED_DEATH_MODES, death_mode_order)
  class_order <- mat %>%
    group_by(therapy_class) %>%
    summarise(total = sum(pmid_count), .groups = "drop") %>%
    arrange(total) %>%
    pull(therapy_class)
  mat <- mat %>%
    mutate(therapy_class_raw = therapy_class) %>%
    complete(
      therapy_class_raw = class_order,
      death_mode = death_mode_order,
      fill = list(
        pmid_count = 0,
        high_confidence_relation_count = 0,
        tumor_family_breadth = 0,
        relation_category_distribution = "",
        therapy_death_selectivity = NA_real_,
        recent_growth = 0,
        low_count_warning_yes_no = "no"
      )
    ) %>%
    mutate(
           therapy_class = factor(class_labels[therapy_class_raw], levels = class_labels[class_order]),
           death_mode = factor(death_mode, levels = death_mode_order),
           clipped_selectivity = pmax(pmin(therapy_death_selectivity, 3), -3)
    )
  p1 <- ggplot(mat, aes(x = death_mode, y = therapy_class, fill = clipped_selectivity)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(
      data = mat %>% filter(high_confidence_relation_count > 0),
      aes(size = high_confidence_relation_count),
      shape = 21, colour = "grey20", fill = NA, stroke = 0.25
    ) +
    scale_fill_gradient2(
      low = "#4C78A8", mid = "white", high = "#B44D3A",
      midpoint = 0, limits = c(-3, 3), na.value = "grey92",
      name = "log2 selectivity"
    ) +
    scale_size_continuous(range = c(0.6, 4), name = "High-confidence relations") +
    labs(title = "A  Therapy class x death-mode relation landscape", x = NULL, y = NULL) +
    theme_oncology_atlas() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  p2 <- mat %>%
    group_by(therapy_class) %>%
    summarise(total_relations = sum(high_confidence_relation_count), tumor_breadth = max(tumor_family_breadth), .groups = "drop") %>%
    ggplot(aes(x = total_relations, y = therapy_class, fill = therapy_class)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = setNames(therapy_palette[names(class_labels)], class_labels[names(class_labels)]), na.value = "grey70") +
    labs(title = "B  Relation support by class", x = "High-confidence relation candidates", y = NULL) +
    theme_oncology_atlas()

  plot <- p1 | p2 +
    plot_annotation(
      title = "Therapy classes differ in death-mode relation-candidate contexts",
      subtitle = "Counts are title/abstract sentence-rule candidates, not evidence of clinical efficacy.",
      theme = theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey30"))
    )
  base <- file.path(fig_dir, "Figure_therapy_class_death_mode_landscape")
  save_pub_r(plot, base, width_mm = 183, height_mm = 105)
  add_log("Figure_therapy_class_death_mode_landscape", "generated", "therapy_class gate passed sufficiently for supplementary figure", base)
}

only_figure <- Sys.getenv("ONCOLOGY_FIGURE_ONLY", unset = "")
run_selected <- function(key) only_figure == "" || only_figure == key

if (run_selected("death_mode_hotness")) {
  if (gate_passes("death_mode_hotness")) make_hotness() else add_log("Figure_oncology_death_mode_hotness", "rejected", "death_mode_hotness gate did not pass")
}
if (run_selected("tumor_specificity")) {
  if (gate_passes("tumor_specificity")) make_tumor_landscape() else add_log("Figure_tumor_specific_death_landscape", "rejected", "tumor_specificity gate did not pass")
}
if (run_selected("drug_death_triads")) {
  if (gate_passes("drug_death_triads")) make_triads() else add_log("Figure_oncology_drug_death_triads", "rejected", "drug_death_triads gate did not pass")
}
if (run_selected("sentence_relation_audit")) {
  if (gate_passes("sentence_relation_audit")) make_sentence_audit() else add_log("Figure_sentence_level_relation_audit", "rejected", "sentence_relation_audit gate did not pass")
}
if (run_selected("same_drug_multi_death")) {
  if (gate_passes("same_drug_multi_death")) make_same_drug() else add_log("Figure_same_drug_multi_death_signals", "rejected", "same_drug_multi_death gate did not pass")
}
if (run_selected("therapy_class_landscape")) {
  if (gate_passes("therapy_class_landscape")) make_therapy_class() else add_log("Figure_therapy_class_death_mode_landscape", "rejected", "therapy_class_landscape gate did not pass")
}

log_suffix <- ifelse(only_figure == "", "", paste0("_", only_figure))
write_csv_safe(figure_log, file.path(fig_dir, paste0("oncology_figure_export_log", log_suffix, ".csv")))

rejected <- figure_log %>%
  filter(status != "generated") %>%
  transmute(figure_or_module = figure, reason_for_rejection_or_demote = reason, recommended_output = "table_only_or_supplementary_table")
write_csv_safe(rejected, file.path(ROOT, "15_QC_reports", paste0("rejected_oncology_figures_from_plotting", log_suffix, ".csv")))
