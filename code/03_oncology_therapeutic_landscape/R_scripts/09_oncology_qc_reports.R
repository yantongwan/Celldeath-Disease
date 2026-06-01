script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

hot <- read_csv_required(file.path(ROOT, "03_oncology_death_mode_hotness/oncology_death_mode_hotness.csv"))
td <- read_csv_required(file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_selectivity.csv"))
poly <- read_csv_required(file.path(ROOT, "07_same_drug_multi_death_signals/drug_polydeath_signal_table.csv"))
gates <- read_csv_required(file.path(ROOT, "15_QC_reports/oncology_analysis_gate_decision_table.csv"))
fig_log <- read_csv_required(file.path(ROOT, "10_publication_grade_figures/oncology_figure_export_log.csv"))

top_hot <- hot %>%
  arrange(heat_rank) %>%
  mutate(line = paste0(heat_rank, ". ", death_mode, " (PMIDs=", oncology_unique_pmids,
                       "; tumor breadth=", tumor_family_breadth,
                       "; heat=", round(heat_index, 2), ")")) %>%
  pull(line)

top_td <- td %>%
  filter(pair_count >= 3) %>%
  arrange(desc(log2_tumor_death_selectivity)) %>%
  slice_head(n = 10) %>%
  mutate(line = paste0("- ", tumor_family, " / ", death_mode,
                       ": log2 selectivity=", round(log2_tumor_death_selectivity, 2),
                       "; pair_count=", pair_count,
                       "; PMIDs=", unique_pmid_count)) %>%
  pull(line)

top_drugs <- poly %>%
  filter(evidence_grade == "displayable_literature_signal") %>%
  arrange(desc(polydeath_index)) %>%
  slice_head(n = 10) %>%
  mutate(line = paste0("- ", drug_display, " (", therapy_class, "): death-mode breadth=", death_mode_breadth,
                       "; tumor breadth=", tumor_family_breadth,
                       "; high-confidence relation candidates=", high_confidence_relation_count,
                       "; dominant context=", dominant_death_mode)) %>%
  pull(line)

summary_lines <- c(
  "# Oncology Therapeutic Death Landscape Summary",
  "",
  "## Scope",
  paste0("- Oncology scope: ", gates$effect_size_summary[gates$module == "oncology_scope"], "."),
  "- Source scope: v5 selected death-mode 2000-2025 audited atlas records.",
  "- Drug extraction tier: title/abstract lexicon and rule-based sentence relation candidates; MeSH chemical fields unavailable.",
  "",
  "## Oncology Death-Mode Hotness Ranking",
  top_hot,
  "",
  "## Top Tumor-Specific Death-Mode Differences",
  top_td,
  "",
  "## Same-Drug Multi-Death Literature Signals",
  top_drugs,
  "",
  "## Gate Decisions",
  paste0("- ", gates$module, ": ", gates$passes_gate_yes_no, " (", gates$effect_size_summary, ")"),
  "",
  "## Recommended Manuscript Integration",
  "- Do not create a new main oncology mechanistic model.",
  "- Add at most one short Results paragraph or place as supplementary oncology diagnostics.",
  "- Tumor-family death-mode profiles are the strongest biology-facing result.",
  "- Drug, triad, same-drug, and therapy-class outputs should remain supplementary until manual validation is completed.",
  "",
  "## Interpretation Boundary",
  "- Report as drug-associated death-mode literature signals, reported induction/sensitization contexts, or hypothesis-generating tumor-drug-death-mode triads.",
  "- Do not convert these diagnostics into validated drug-induction, clinical, or tumor-vulnerability claims."
)
write_text_safe(summary_lines, file.path(ROOT, "15_QC_reports/oncology_final_module_summary.md"))

fig_qc <- c(
  "# Oncology Figure Readability QC",
  "",
  paste0("- Figures generated: ", sum(fig_log$status == "generated"), "."),
  paste0("- Figures rejected/demoted at plotting stage: ", sum(fig_log$status != "generated"), "."),
  "",
  "## Generated Figures",
  paste0("- ", fig_log$figure[fig_log$status == "generated"], ": PDF/PNG/TIFF exported."),
  "",
  "## Visual Boundary Notes",
  "- Hotness figure is component-wise and descriptive; use as supplementary unless the manuscript explicitly needs an oncology context figure.",
  "- Tumor-specific death landscape is the most suitable figure for a concise main-text oncology paragraph.",
  "- Drug triad, same-drug, therapy-class, and sentence-relation figures should remain supplementary because relation evidence is rule-based and not manually validated.",
  "- No figure should be captioned as drug-induced pathway proof."
)
write_text_safe(fig_qc, file.path(ROOT, "15_QC_reports/oncology_figure_readability_QC.md"))

rejected_lines <- c(
  "# Rejected Or Demoted Oncology Figures",
  "",
  "## Demoted From Main Text",
  "- Composite main oncology therapeutic figure: demoted to supplementary diagnostics because drug and relation evidence is rule-based and manual validation remains incomplete.",
  "- Oncology therapeutic ecology network: not generated because network-style presentation could visually overclaim tumor-drug-death-mode mechanisms and would likely be dense.",
  "",
  "## Generated But Recommended Supplementary",
  "- Figure_oncology_death_mode_hotness: descriptive literature heat index; rank sensitivity was moderate rather than definitive.",
  "- Figure_oncology_drug_death_triads: hypothesis-generating triads only.",
  "- Figure_same_drug_multi_death_signals: same-drug recurrence can reflect publication volume or tool-compound reuse.",
  "- Figure_therapy_class_death_mode_landscape: relation-candidate context only.",
  "- Figure_sentence_level_relation_audit: audit support, not biological result.",
  "",
  "## Safe Main-Text Candidate",
  "- Figure_tumor_specific_death_landscape can support a short main-text oncology heterogeneity paragraph if the manuscript has room, but it should not replace the existing main figure hierarchy."
)
write_text_safe(rejected_lines, file.path(ROOT, "15_QC_reports/rejected_oncology_figures.md"))
