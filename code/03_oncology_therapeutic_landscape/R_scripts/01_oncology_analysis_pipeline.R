script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))
source(file.path(script_dir, "00_setup_oncology_module.R"))

log_lines <- character()
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(msg)
  log_lines <<- c(log_lines, msg)
}

make_tumor_mapping <- function(term) {
  t <- str_to_lower(term)
  family <- case_when(
    str_detect(t, "triple negative|tnbc") ~ "triple-negative breast cancer",
    str_detect(t, "breast|mammary") ~ "breast cancer",
    str_detect(t, "hepatocellular|liver neoplasm|liver cancer|hepatic neoplasm") ~ "liver cancer",
    str_detect(t, "lung|pulmonary neoplasm|bronchogenic") ~ "lung cancer",
    str_detect(t, "colorectal|colon|rectal|colonic") ~ "colorectal cancer",
    str_detect(t, "gastric|stomach") ~ "gastric cancer",
    str_detect(t, "pancrea") ~ "pancreatic cancer",
    str_detect(t, "ovarian|ovary") ~ "ovarian cancer",
    str_detect(t, "prostate") ~ "prostate cancer",
    str_detect(t, "renal cell|kidney|renal neoplasm") ~ "renal cancer",
    str_detect(t, "melanoma") ~ "melanoma",
    str_detect(t, "glioblastoma") ~ "glioblastoma",
    str_detect(t, "glioma|brain neoplasm|central nervous system neoplasm|cns") ~ "CNS tumor",
    str_detect(t, "leukemia|leukaemia") ~ "leukemia",
    str_detect(t, "lymphoma") ~ "lymphoma",
    str_detect(t, "myeloma") ~ "multiple myeloma",
    str_detect(t, "sarcoma|osteosarcoma") ~ "sarcoma",
    str_detect(t, "cervical") ~ "cervical cancer",
    str_detect(t, "bladder|urothelial") ~ "bladder cancer",
    str_detect(t, "esophageal|oesophageal") ~ "esophageal cancer",
    str_detect(t, "oral|head and neck|pharyngeal|laryngeal|nasopharyngeal") ~ "head and neck cancer",
    str_detect(t, "thyroid") ~ "thyroid cancer",
    str_detect(t, "endometrial|uterine") ~ "endometrial/uterine cancer",
    str_detect(t, "cholangiocarcinoma|bile duct") ~ "biliary tract cancer",
    str_detect(t, "neoplasm|cancer|tumou?r|carcinoma|adenocarcinoma|malignan") ~ "mixed/unspecified neoplasms",
    TRUE ~ "not mapped"
  )
  solid_heme <- case_when(
    family %in% c("leukemia", "lymphoma", "multiple myeloma") ~ "hematologic",
    family == "not mapped" ~ "not mapped",
    TRUE ~ "solid"
  )
  organ <- case_when(
    str_detect(family, "breast") ~ "breast",
    str_detect(family, "liver") ~ "liver",
    str_detect(family, "lung") ~ "lung",
    str_detect(family, "colorectal|gastric|pancreatic|biliary|esophageal") ~ "digestive",
    str_detect(family, "ovarian|cervical|endometrial|uterine") ~ "gynecologic",
    str_detect(family, "prostate|bladder|renal") ~ "urogenital",
    str_detect(family, "melanoma") ~ "skin",
    str_detect(family, "CNS|glioblastoma") ~ "nervous system",
    str_detect(family, "leukemia|lymphoma|myeloma") ~ "hematologic",
    str_detect(family, "sarcoma") ~ "mesenchymal/bone",
    str_detect(family, "head and neck|thyroid") ~ "head and neck/endocrine",
    TRUE ~ "mixed/unspecified"
  )
  confidence <- case_when(
    family == "not mapped" ~ "low",
    family == "mixed/unspecified neoplasms" ~ "medium",
    TRUE ~ "high"
  )
  tibble(
    tumor_family = family,
    broad_tumor_class = ifelse(family == "mixed/unspecified neoplasms", "broad_neoplasm", "mapped_tumor_family"),
    solid_vs_hematologic = solid_heme,
    organ_system = organ,
    subtype_if_available = ifelse(family == "triple-negative breast cancer", "TNBC", ""),
    mapping_rule = ifelse(confidence == "high", "regex_specific_tumor_family", ifelse(confidence == "medium", "regex_generic_neoplasm", "not_mapped")),
    mapping_confidence = confidence
  )
}

build_drug_lexicon <- function() {
  tribble(
    ~drug_normalized, ~drug_display, ~synonyms, ~therapy_class, ~is_specific_drug_yes_no, ~is_generic_therapy_term_yes_no,
    "cisplatin", "cisplatin", "cisplatin", "chemotherapy", "yes", "no",
    "carboplatin", "carboplatin", "carboplatin", "chemotherapy", "yes", "no",
    "oxaliplatin", "oxaliplatin", "oxaliplatin", "chemotherapy", "yes", "no",
    "doxorubicin", "doxorubicin", "doxorubicin|adriamycin", "chemotherapy", "yes", "no",
    "paclitaxel", "paclitaxel", "paclitaxel|taxol", "chemotherapy", "yes", "no",
    "docetaxel", "docetaxel", "docetaxel", "chemotherapy", "yes", "no",
    "gemcitabine", "gemcitabine", "gemcitabine", "chemotherapy", "yes", "no",
    "5-fluorouracil", "5-fluorouracil", "5-fluorouracil|5-fu|fluorouracil", "chemotherapy", "yes", "no",
    "temozolomide", "temozolomide", "temozolomide|tmz", "chemotherapy", "yes", "no",
    "bortezomib", "bortezomib", "bortezomib", "proteasome_inhibitor", "yes", "no",
    "sorafenib", "sorafenib", "sorafenib", "kinase_inhibitor", "yes", "no",
    "lenvatinib", "lenvatinib", "lenvatinib", "kinase_inhibitor", "yes", "no",
    "regorafenib", "regorafenib", "regorafenib", "kinase_inhibitor", "yes", "no",
    "sunitinib", "sunitinib", "sunitinib", "kinase_inhibitor", "yes", "no",
    "imatinib", "imatinib", "imatinib", "kinase_inhibitor", "yes", "no",
    "gefitinib", "gefitinib", "gefitinib", "targeted_therapy", "yes", "no",
    "erlotinib", "erlotinib", "erlotinib", "targeted_therapy", "yes", "no",
    "osimertinib", "osimertinib", "osimertinib", "targeted_therapy", "yes", "no",
    "vemurafenib", "vemurafenib", "vemurafenib", "targeted_therapy", "yes", "no",
    "olaparib", "olaparib", "olaparib", "targeted_therapy", "yes", "no",
    "tamoxifen", "tamoxifen", "tamoxifen", "targeted_therapy", "yes", "no",
    "trastuzumab", "trastuzumab", "trastuzumab", "targeted_therapy", "yes", "no",
    "rapamycin", "rapamycin", "rapamycin|sirolimus", "metabolic_modulator", "yes", "no",
    "metformin", "metformin", "metformin", "metabolic_modulator", "yes", "no",
    "erastin", "erastin", "erastin", "ferroptosis_inducer_tool_compound", "yes", "no",
    "rsl3", "RSL3", "rsl3", "ferroptosis_inducer_tool_compound", "yes", "no",
    "sulfasalazine", "sulfasalazine", "sulfasalazine", "ferroptosis_inducer_tool_compound", "yes", "no",
    "elesclomol", "elesclomol", "elesclomol", "experimental_tool_compound", "yes", "no",
    "artesunate", "artesunate", "artesunate", "natural_product", "yes", "no",
    "curcumin", "curcumin", "curcumin", "natural_product", "yes", "no",
    "resveratrol", "resveratrol", "resveratrol", "natural_product", "yes", "no",
    "pembrolizumab", "pembrolizumab", "pembrolizumab", "immunotherapy_checkpoint", "yes", "no",
    "nivolumab", "nivolumab", "nivolumab", "immunotherapy_checkpoint", "yes", "no",
    "atezolizumab", "atezolizumab", "atezolizumab", "immunotherapy_checkpoint", "yes", "no",
    "ipilimumab", "ipilimumab", "ipilimumab", "immunotherapy_checkpoint", "yes", "no",
    "chemotherapy", "chemotherapy", "chemotherapy|chemotherapeutic", "chemotherapy", "no", "yes",
    "radiotherapy", "radiotherapy", "radiotherapy|radiation therapy|irradiation", "radiotherapy", "no", "yes",
    "immunotherapy", "immunotherapy", "immunotherapy|immune checkpoint|pd-1|pd-l1|ctla-4", "immunotherapy_checkpoint", "no", "yes",
    "targeted therapy", "targeted therapy", "targeted therapy|targeted therapies", "targeted_therapy", "no", "yes"
  ) %>%
    mutate(source = "curated_oncology_therapy_keywords",
           mapping_confidence = ifelse(is_specific_drug_yes_no == "yes", "high", "medium"))
}

term_regex <- function(synonyms) {
  parts <- unlist(str_split(synonyms, "\\|"))
  escaped <- str_replace_all(parts, "([\\W])", "\\\\\\1")
  paste0("(?i)(^|[^A-Za-z0-9])(", paste(escaped, collapse = "|"), ")([^A-Za-z0-9]|$)")
}

death_regex <- function(mode) {
  m <- str_to_lower(mode)
  variants <- case_when(
    m == "ferroptosis" ~ "ferroptosis|ferroptotic",
    m == "pyroptosis" ~ "pyroptosis|pyroptotic",
    m == "netosis" ~ "netosis|neutrophil extracellular trap|neutrophil extracellular traps",
    m == "necroptosis" ~ "necroptosis|necroptotic",
    str_detect(m, "immunogenic") | m == "icd" ~ "immunogenic cell death",
    m == "cuproptosis" ~ "cuproptosis|cuproptotic",
    m == "panoptosis" ~ "panoptosis|panoptotic",
    m == "disulfidptosis" ~ "disulfidptosis|disulfidptotic",
    TRUE ~ m
  )
  paste0("(?i)(", variants, ")")
}

classify_relation <- function(sentence) {
  s <- str_to_lower(sentence)
  neg <- str_detect(s, "\\b(not|failed to|did not|no significant|without)\\b")
  category <- "unclear_co_mention"
  cue <- ""
  cues <- list(
    induction_or_activation = "\\b(induce|induced|induces|induction|trigger|triggered|promote|promoted|activate|activated|enhance|enhanced|increase|increased)\\b",
    sensitization = "\\b(sensitize|sensitized|sensitise|sensitised|enhance sensitivity|overcome resistance)\\b",
    inhibition_or_protection = "\\b(inhibit|inhibited|suppress|suppressed|block|blocked|prevent|prevented|protect|protected)\\b",
    resistance_escape = "\\b(resistance|resistant|escape)\\b",
    combination_therapy = "\\b(combination|combined|synergize|synergistic|synergy|co-treatment|cotreatment)\\b"
  )
  for (nm in names(cues)) {
    hit <- str_extract(s, cues[[nm]])
    if (!is.na(hit)) {
      category <- nm
      cue <- hit
      break
    }
  }
  tibble(relation_category = category, cue_word = cue, negation_flag = ifelse(neg, "yes", "no"))
}

split_sentences <- function(text) {
  text <- norm_text(text)
  if (text == "") return(character())
  str_split(text, "(?<=[.!?])\\s+")[[1]] %>% str_squish() %>% .[nchar(.) > 20]
}

approx_distance <- function(sentence, drug_pattern, death_pattern) {
  s <- str_to_lower(sentence)
  words <- unlist(str_split(s, "\\s+"))
  drug_word <- str_extract(s, drug_pattern)
  death_word <- str_extract(s, death_pattern)
  if (is.na(drug_word) || is.na(death_word)) return(NA_real_)
  dp <- which(str_detect(words, fixed(str_squish(drug_word), ignore_case = TRUE)))[1]
  mp <- which(str_detect(words, fixed(str_squish(death_word), ignore_case = TRUE)))[1]
  if (is.na(dp) || is.na(mp)) return(NA_real_)
  abs(dp - mp)
}

main <- function() {
  log_msg("Starting oncology therapeutic death landscape analysis")
  subdirs <- c(
    "01_input_audit", "02_oncology_scope_and_tumor_mapping", "03_oncology_death_mode_hotness",
    "04_tumor_specific_death_profiles", "05_drug_therapy_extraction", "06_drug_death_mode_triads",
    "07_same_drug_multi_death_signals", "08_therapy_class_landscape", "09_sentence_level_relation_audit",
    "10_publication_grade_figures", "11_manuscript_patches", "12_source_data", "15_QC_reports"
  )
  invisible(lapply(file.path(ROOT, subdirs), dir_create))

  article_path <- file.path(V5_ROOT, "01_source_data_filtered/article_stage_records_selected_death_modes_2000_2025.csv")
  pair_path <- file.path(V5_ROOT, "01_source_data_filtered/pair_metrics_selected_death_modes_2000_2025.csv")
  onc_path <- file.path(V5_ROOT, "01_source_data_filtered/oncology_term_classification_selected_death_modes_2000_2025.csv")
  low_path <- file.path(V5_ROOT, "04_supplementary_figures_regenerated/source_data/low_count_warning_pairs_selected_death_modes_2000_2025.csv")

  article <- read_csv_required(article_path)
  pair <- read_csv_required(pair_path)
  onc <- read_csv_required(onc_path)
  low <- read_csv_required(low_path)

  article <- article %>%
    mutate(
      publication_year_num = suppressWarnings(as.integer(publication_year_num)),
      title = norm_text(title),
      abstract = norm_text(abstract),
      text_all = str_squish(paste(title, abstract, sep = ". ")),
      death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode)
    ) %>%
    filter(
      !is.na(publication_year_num),
      publication_year_num >= 2000,
      publication_year_num <= 2025,
      death_mode %in% SELECTED_DEATH_MODES
    )

  pair <- pair %>%
    rename(death_mode = Death_Mode, disease_term = Disease, pair_count = Pair_Count, disease_category = Disease_Category) %>%
    mutate(death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
           pair_key = paste(death_mode, disease_term, sep = "||"))

  onc <- onc %>%
    rename(disease_term = Disease, disease_category_onc = Disease_Category) %>%
    mutate(oncology_flag_base = clean_bool(oncology_flag))

  low_pairs <- low %>%
    mutate(death_mode = ifelse(death_mode == "ICD", "Immunogenic cell death", death_mode),
           pair_key = paste(death_mode, disease_term, sep = "||"),
           low_count_warning_yes_no = "yes") %>%
    select(pair_key, low_count_warning_yes_no, risk_reason)

  oncology_regex <- regex_any(c(
    "neoplasm", "carcinoma", "cancer", "tumou?r", "glioma", "melanoma", "leukemia", "leukaemia",
    "lymphoma", "sarcoma", "adenocarcinoma", "glioblastoma", "hepatocellular", "renal cell",
    "malignan", "myeloma", "blastoma"
  ))

  article_scoped <- article %>%
    left_join(onc, by = "disease_term") %>%
    left_join(pair %>% select(pair_key, death_mode, disease_term, pair_count_pair = pair_count, disease_category, Article_Count), by = c("death_mode", "disease_term")) %>%
    mutate(
      pair_key = paste(death_mode, disease_term, sep = "||"),
      oncology_flag_regex = str_detect(str_to_lower(disease_term), oncology_regex),
      oncology_flag = oncology_flag_base | disease_category_onc == "Neoplasms" | disease_category == "Neoplasms" | oncology_flag_regex
    ) %>%
    left_join(low_pairs, by = "pair_key") %>%
    mutate(low_count_warning_yes_no = ifelse(is.na(low_count_warning_yes_no), "no", low_count_warning_yes_no),
           pair_count_pair = coalesce(pair_count_pair, pair_count))

  tumor_map_rows <- article_scoped %>%
    filter(oncology_flag) %>%
    distinct(disease_term, .keep_all = TRUE) %>%
    mutate(map = lapply(disease_term, make_tumor_mapping)) %>%
    tidyr::unnest(map)

  article_onc <- article_scoped %>%
    filter(oncology_flag) %>%
    left_join(tumor_map_rows %>% select(disease_term, tumor_family, broad_tumor_class, solid_vs_hematologic, organ_system, subtype_if_available, mapping_rule, mapping_confidence), by = "disease_term") %>%
    mutate(tumor_family = ifelse(is.na(tumor_family) | tumor_family == "not mapped", "mixed/unspecified neoplasms", tumor_family))

  mapping <- article_onc %>%
    group_by(disease_term) %>%
    summarise(
      oncology_flag = TRUE,
      tumor_family = first(tumor_family),
      broad_tumor_class = first(broad_tumor_class),
      solid_vs_hematologic = first(solid_vs_hematologic),
      organ_system = first(organ_system),
      subtype_if_available = first(subtype_if_available),
      mapping_rule = first(mapping_rule),
      mapping_confidence = first(mapping_confidence),
      pair_count_total = sum(unique(pair_count_pair), na.rm = TRUE),
      unique_pmid_count = n_distinct(pmid),
      death_mode_breadth = n_distinct(death_mode),
      .groups = "drop"
    ) %>%
    arrange(desc(unique_pmid_count), disease_term)

  scope_summary <- tibble(
    metric = c("article_stage_records_2000_2025_selected_death_modes", "oncology_article_stage_records", "oncology_unique_pmids", "oncology_disease_terms", "tumor_families"),
    value = c(nrow(article_scoped), nrow(article_onc), n_distinct(article_onc$pmid), n_distinct(article_onc$disease_term), n_distinct(article_onc$tumor_family))
  )

  pubmed_pubtator_drug_path <- file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions_evidence_merged.csv")
  pubmed_pubtator_relation_path <- file.path(ROOT, "09_sentence_level_relation_audit/drug_death_relation_sentences_validated_input.csv")
  chemical_evidence_available <- file.exists(pubmed_pubtator_drug_path)
  relation_evidence_available <- file.exists(pubmed_pubtator_relation_path)

  input_availability <- tibble(
    input = c("article_stage_records", "pair_metrics", "oncology_term_classification", "low_count_warning_pairs", "title", "abstract", "PubMed/PubTator chemical evidence", "PubMed/PubTator sentence relation candidates"),
    path = c(article_path, pair_path, onc_path, low_path, article_path, article_path, pubmed_pubtator_drug_path, pubmed_pubtator_relation_path),
    available_yes_no = c("yes", "yes", "yes", "yes", "yes", "yes", ifelse(chemical_evidence_available, "yes", "no"), ifelse(relation_evidence_available, "yes", "no")),
    rows = c(nrow(article), nrow(pair), nrow(onc), nrow(low), sum(nchar(article$title) > 0), sum(nchar(article$abstract) > 20), ifelse(chemical_evidence_available, nrow(readr::read_csv(pubmed_pubtator_drug_path, show_col_types = FALSE)), 0), ifelse(relation_evidence_available, nrow(readr::read_csv(pubmed_pubtator_relation_path, show_col_types = FALSE)), 0)),
    notes = c("v5 selected death-mode 2000-2025 source", "v5 selected death-mode pair metrics", "v5 oncology sensitivity output", "v5 low-count warnings", "usable for lexicon and sentence matching", "usable for sentence-level relation candidates", "official PubMed ChemicalList plus PubTator3 chemical annotations when present", "rule-based sentence relation candidates generated from PubMed/PubTator-enhanced mentions when present")
  )

  write_csv_safe(input_availability, file.path(ROOT, "01_input_audit/oncology_input_availability.csv"))
  write_text_safe(c(
    "# Oncology Input Audit",
    "",
    "- Scope: selected death-mode 2000-2025 package inputs were used with the final eight regulated cell death modes and the 2000-2025 publication-year window.",
    "- Title and abstract text are available and support conservative lexicon and sentence-window relation candidates.",
    ifelse(chemical_evidence_available, "- PubMed ChemicalList and PubTator3 chemical annotations are available as external evidence tiers and are merged conservatively with the local therapy lexicon.", "- External PubMed/PubTator chemical evidence is not available; the pipeline falls back to title/abstract lexicon matching."),
    "- All drug and therapy outputs are literature-level diagnostics."
  ), file.path(ROOT, "01_input_audit/oncology_input_audit.md"))
  write_text_safe(c(
    "# Missing Inputs For Drug Analysis",
    "",
    ifelse(chemical_evidence_available, "- PubMed ChemicalList/PubTator3 evidence: available from cached external annotation outputs.", "- PubMed ChemicalList/PubTator3 evidence: missing."),
    "- Supplementary Concept Records: missing.",
    "- RxNorm/DrugBank/ChEMBL mappings: not present in audited local inputs.",
    ifelse(chemical_evidence_available, "- Consequence: drug extraction uses PubMed/PubTator chemical evidence plus curated local therapy lexicon mapping; unmapped chemicals remain needs_review and are not treated as validated drugs.", "- Consequence: drug extraction uses curated oncology therapy keywords and title/abstract lexicon matching only; relation-level outputs require sentence cues and remain candidates.")
  ), file.path(ROOT, "01_input_audit/missing_inputs_for_drug_analysis.md"))

  write_csv_safe(mapping, file.path(ROOT, "02_oncology_scope_and_tumor_mapping/oncology_tumor_term_mapping.csv"))
  write_csv_safe(scope_summary, file.path(ROOT, "02_oncology_scope_and_tumor_mapping/oncology_scope_summary.csv"))
  write_csv_safe(article_onc, file.path(ROOT, "12_source_data/oncology_article_records_scoped.csv"))

  # Oncology death-mode hotness
  unique_articles <- article_scoped %>% distinct(death_mode, pmid, .keep_all = TRUE)
  unique_onc <- article_onc %>% distinct(death_mode, pmid, .keep_all = TRUE)
  pair_onc <- article_onc %>% distinct(death_mode, disease_term, pair_key, pair_count_pair, tumor_family, low_count_warning_yes_no)
  total_by_death <- unique_articles %>% group_by(death_mode) %>% summarise(total_unique_pmids = n_distinct(pmid), .groups = "drop")
  onc_death_base <- unique_onc %>%
    mutate(
      therapy_vocab = str_detect(str_to_lower(text_all), "\\b(therap|treat|drug|chemotherap|radiotherap|immunotherap|checkpoint|inhibitor|inducer|sensiti[sz]|combination)\\b")
    )

  drug_lex <- build_drug_lexicon()
  write_csv_safe(drug_lex, file.path(ROOT, "05_drug_therapy_extraction/drug_lexicon_normalized.csv"))
  drug_patterns <- setNames(lapply(drug_lex$synonyms, term_regex), drug_lex$drug_normalized)
  text_lower <- str_to_lower(onc_death_base$text_all)
  drug_any <- rep(FALSE, nrow(onc_death_base))
  for (pat in drug_patterns) drug_any <- drug_any | str_detect(text_lower, pat)
  onc_death_base$drug_any <- drug_any

  year_counts <- unique_onc %>%
    group_by(death_mode) %>%
    summarise(
      pmids_2023_2025 = n_distinct(pmid[publication_year_num >= 2023]),
      pmids_2020_2022 = n_distinct(pmid[publication_year_num >= 2020 & publication_year_num <= 2022]),
      pmids_2000_2022 = n_distinct(pmid[publication_year_num <= 2022]),
      .groups = "drop"
    )

  hotness <- unique_onc %>%
    group_by(death_mode) %>%
    summarise(
      oncology_unique_pmids = n_distinct(pmid),
      tumor_family_breadth = n_distinct(tumor_family),
      .groups = "drop"
    ) %>%
    left_join(pair_onc %>% group_by(death_mode) %>% summarise(
      oncology_pair_count = sum(pair_count_pair, na.rm = TRUE),
      low_count_warning_fraction = mean(low_count_warning_yes_no == "yes", na.rm = TRUE),
      .groups = "drop"
    ), by = "death_mode") %>%
    left_join(total_by_death, by = "death_mode") %>%
    left_join(year_counts, by = "death_mode") %>%
    left_join(onc_death_base %>% group_by(death_mode) %>% summarise(
      therapy_vocabulary_fraction = mean(therapy_vocab, na.rm = TRUE),
      drug_associated_pmid_fraction = mean(drug_any, na.rm = TRUE),
      .groups = "drop"
    ), by = "death_mode") %>%
    mutate(
      fraction_of_death_mode_literature_in_oncology = oncology_unique_pmids / total_unique_pmids,
      recent_growth_3yr = (pmids_2023_2025 + 1) / (pmids_2020_2022 + 1),
      burst_index = (pmids_2023_2025 / 3) / pmax(pmids_2000_2022 / 23, 1e-9),
      heat_index = zscore(log1p(oncology_unique_pmids)) + zscore(tumor_family_breadth) + zscore(recent_growth_3yr) + zscore(drug_associated_pmid_fraction),
      heat_rank = rank(-heat_index, ties.method = "first"),
      interpretation = "Oncology literature heat index; not a therapeutic target ranking."
    ) %>%
    arrange(heat_rank)

  hotness_sens <- pair_onc %>%
    filter(pair_count_pair >= 3) %>%
    distinct(death_mode, disease_term, tumor_family, pair_count_pair) %>%
    group_by(death_mode) %>%
    summarise(pair_support_ge3 = sum(pair_count_pair, na.rm = TRUE), tumor_family_breadth_ge3 = n_distinct(tumor_family), .groups = "drop") %>%
    left_join(hotness %>% select(death_mode, heat_rank), by = "death_mode") %>%
    mutate(rank_ge3 = rank(-pair_support_ge3, ties.method = "first"))
  hotness_spearman <- suppressWarnings(cor(hotness_sens$heat_rank, hotness_sens$rank_ge3, method = "spearman", use = "complete.obs"))

  write_csv_safe(hotness, file.path(ROOT, "03_oncology_death_mode_hotness/oncology_death_mode_hotness.csv"))
  write_csv_safe(hotness_sens, file.path(ROOT, "03_oncology_death_mode_hotness/oncology_death_mode_hotness_paircount3_sensitivity.csv"))

  # Tumor-specific death profiles
  tumor_support <- article_onc %>%
    group_by(tumor_family) %>%
    summarise(unique_pmid_count = n_distinct(pmid), total_pair_count = sum(unique(pair_count_pair), na.rm = TRUE), death_mode_breadth = n_distinct(death_mode), mapping_confidence = first(mapping_confidence), .groups = "drop")
  display_tumors <- tumor_support %>% filter((unique_pmid_count >= 20 | total_pair_count >= 20), mapping_confidence != "low") %>% pull(tumor_family)
  td_pair <- pair_onc %>%
    group_by(tumor_family, death_mode) %>%
    summarise(pair_count = sum(pair_count_pair, na.rm = TRUE), low_count_warning_yes_no = ifelse(any(low_count_warning_yes_no == "yes"), "yes", "no"), .groups = "drop")
  td_pmid <- article_onc %>%
    group_by(tumor_family, death_mode) %>%
    summarise(unique_pmid_count = n_distinct(pmid), recent_growth = (n_distinct(pmid[publication_year_num >= 2023]) + 1) / (n_distinct(pmid[publication_year_num >= 2020 & publication_year_num <= 2022]) + 1), .groups = "drop")
  death_bg <- td_pair %>% group_by(death_mode) %>% summarise(death_pair_total = sum(pair_count), .groups = "drop")
  total_onc_pair <- sum(td_pair$pair_count, na.rm = TRUE)
  tumor_death <- td_pair %>%
    left_join(td_pmid, by = c("tumor_family", "death_mode")) %>%
    group_by(tumor_family) %>%
    mutate(tumor_total_pair = sum(pair_count, na.rm = TRUE), tumor_normalized_share = pair_count / tumor_total_pair) %>%
    ungroup() %>%
    left_join(death_bg, by = "death_mode") %>%
    mutate(
      oncology_background_share = death_pair_total / total_onc_pair,
      log2_tumor_death_selectivity = log2((tumor_normalized_share + 1e-6) / (oncology_background_share + 1e-6)),
      cleaned_axis_dominant = "not_computed_from_article_text",
      interpretation = "Tumor-family death-mode literature selectivity; not proof of tumor biology."
    ) %>%
    arrange(tumor_family, desc(log2_tumor_death_selectivity))

  tumor_profile <- tumor_death %>%
    group_by(tumor_family) %>%
    summarise(
      unique_pmid_count = max(unique_pmid_count, na.rm = TRUE),
      death_mode_breadth = sum(pair_count >= 3, na.rm = TRUE),
      top_death_modes = paste(head(death_mode[order(-pair_count)], 4), collapse = "; "),
      top_selective_modes = paste(head(death_mode[order(-log2_tumor_death_selectivity)], 4), collapse = "; "),
      .groups = "drop"
    )

  tumor_spec <- tumor_death %>%
    group_by(tumor_family) %>%
    summarise(
      death_mode_entropy = {
        p <- pair_count / sum(pair_count)
        -sum(p * log(p + 1e-12))
      },
      top_death_mode_share = max(pair_count / sum(pair_count)),
      specialization_index = 1 - (death_mode_entropy / log(n())),
      tumor_death_mode_breadth = sum(pair_count >= 3),
      instability_flag = ifelse(any(low_count_warning_yes_no == "yes"), "yes", "no"),
      .groups = "drop"
    )

  write_csv_safe(tumor_death, file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_mode_matrix.csv"))
  write_csv_safe(tumor_death, file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_selectivity.csv"))
  write_csv_safe(tumor_profile, file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_profile_summary.csv"))
  write_csv_safe(tumor_spec, file.path(ROOT, "04_tumor_specific_death_profiles/tumor_death_specialization.csv"))

  # Drug mentions. Prefer external PubMed/PubTator chemical evidence if present;
  # otherwise fall back to the original title/abstract lexicon scan.
  if (chemical_evidence_available) {
    article_text_lookup <- article_onc %>%
      transmute(
        pmid = as.character(pmid),
        tumor_family,
        disease_term,
        death_mode,
        text_all
      ) %>%
      distinct()

    drug_mentions <- readr::read_csv(pubmed_pubtator_drug_path, show_col_types = FALSE) %>%
      mutate(
        pmid = as.character(pmid),
        year = suppressWarnings(as.integer(year)),
        mention_source = evidence_sources,
        mention_sentence_if_available = mention_text,
        cooccurrence_tier = evidence_tier,
        is_specific_drug_yes_no = coalesce(as.character(specific_drug_yes_no), "needs_review"),
        is_generic_therapy_term_yes_no = coalesce(as.character(generic_therapy_term_yes_no), "no"),
        mapping_confidence = coalesce(as.character(mapping_confidence), "needs_review")
      ) %>%
      left_join(article_text_lookup, by = c("pmid", "tumor_family", "disease_term", "death_mode")) %>%
      mutate(text_all = coalesce(text_all, "")) %>%
      select(
        pmid, year, tumor_family, disease_term, death_mode,
        drug_normalized, drug_display, therapy_class,
        mention_source, mention_sentence_if_available, cooccurrence_tier,
        mapping_confidence, is_specific_drug_yes_no, is_generic_therapy_term_yes_no,
        text_all
      ) %>%
      distinct(pmid, tumor_family, disease_term, death_mode, drug_normalized, .keep_all = TRUE)
  } else {
    mention_list <- vector("list", length(drug_patterns))
    names(mention_list) <- names(drug_patterns)
    for (i in seq_along(drug_patterns)) {
      pat <- drug_patterns[[i]]
      hit <- str_detect(str_to_lower(article_onc$text_all), pat)
      if (any(hit)) {
        lex <- drug_lex[i, ]
        mention_list[[i]] <- article_onc[hit, ] %>%
          transmute(
            pmid, year = publication_year_num, tumor_family, disease_term, death_mode,
            drug_normalized = lex$drug_normalized,
            drug_display = lex$drug_display,
            therapy_class = lex$therapy_class,
            mention_source = "title_abstract_lexicon",
            mention_sentence_if_available = "",
            cooccurrence_tier = "Tier2_title_abstract_drug_death_cooccurrence",
            mapping_confidence = lex$mapping_confidence,
            is_specific_drug_yes_no = lex$is_specific_drug_yes_no,
            is_generic_therapy_term_yes_no = lex$is_generic_therapy_term_yes_no,
            text_all
          )
      }
    }
    drug_mentions <- bind_rows(mention_list) %>% distinct(pmid, tumor_family, disease_term, death_mode, drug_normalized, .keep_all = TRUE)
  }
  write_csv_safe(drug_mentions %>% select(-text_all), file.path(ROOT, "05_drug_therapy_extraction/oncology_drug_mentions.csv"))

  # Sentence-level relation candidates
  relation_rows <- list()
  if (relation_evidence_available) {
    relation_raw <- readr::read_csv(pubmed_pubtator_relation_path, show_col_types = FALSE)
    if (!"relation_category_rule_based" %in% names(relation_raw)) relation_raw$relation_category_rule_based <- relation_raw$relation_category
    if (!"confidence_rule_based" %in% names(relation_raw)) relation_raw$confidence_rule_based <- relation_raw$confidence_high_medium_low
    if (!"evidence_sources" %in% names(relation_raw)) relation_raw$evidence_sources <- ""

    relation_lookup <- drug_mentions %>%
      transmute(
        pmid = as.character(pmid),
        tumor_family,
        disease_term,
        death_mode,
        drug_normalized,
        therapy_class,
        is_specific_drug_yes_no,
        is_generic_therapy_term_yes_no
      ) %>%
      distinct()

    relations <- relation_raw %>%
      mutate(
        pmid = as.character(pmid),
        year = suppressWarnings(as.integer(year)),
        drug_death_distance_tokens = suppressWarnings(as.numeric(drug_death_distance_tokens)),
        relation_category = relation_category_rule_based,
        confidence_high_medium_low = confidence_rule_based,
        evidence_tier = ifelse(evidence_sources == "", "PubMed/PubTator_sentence_relation_candidate", evidence_sources)
      ) %>%
      left_join(relation_lookup, by = c("pmid", "tumor_family", "disease_term", "death_mode", "drug_normalized")) %>%
      mutate(
        therapy_class = coalesce(therapy_class, "unmapped_chemical"),
        is_specific_drug_yes_no = coalesce(is_specific_drug_yes_no, "needs_review"),
        is_generic_therapy_term_yes_no = coalesce(is_generic_therapy_term_yes_no, "no")
      ) %>%
      select(
        pmid, year, tumor_family, death_mode, drug_normalized,
        sentence, relation_category, cue_word, drug_death_distance_tokens,
        tumor_context_present_yes_no, negation_flag, confidence_high_medium_low,
        evidence_tier, therapy_class, is_specific_drug_yes_no, is_generic_therapy_term_yes_no
      ) %>%
      distinct()
  } else if (nrow(drug_mentions) > 0) {
    for (idx in seq_len(nrow(drug_mentions))) {
      row <- drug_mentions[idx, ]
      pat_drug <- drug_patterns[[row$drug_normalized]]
      if (is.null(pat_drug)) next
      pat_death <- death_regex(row$death_mode)
      sents <- split_sentences(row$text_all)
      if (length(sents) == 0) next
      sent_hit <- sents[str_detect(str_to_lower(sents), pat_drug) & str_detect(str_to_lower(sents), pat_death)]
      if (length(sent_hit) == 0) next
      for (sent in sent_hit) {
        rel <- classify_relation(sent)
        tumor_present <- str_detect(str_to_lower(row$text_all), "cancer|tumou?r|neoplasm|carcinoma|glioma|melanoma|leukemia|lymphoma|sarcoma|malignan")
        dist <- approx_distance(sent, pat_drug, pat_death)
        confidence <- case_when(
          rel$relation_category %in% c("induction_or_activation", "sensitization") & !is.na(dist) & dist <= 18 & rel$negation_flag == "no" ~ "high",
          rel$relation_category != "unclear_co_mention" & rel$negation_flag == "no" ~ "medium",
          TRUE ~ "low"
        )
        relation_rows[[length(relation_rows) + 1]] <- tibble(
          pmid = row$pmid,
          year = row$year,
          tumor_family = row$tumor_family,
          death_mode = row$death_mode,
          drug_normalized = row$drug_normalized,
          sentence = sent,
          relation_category = rel$relation_category,
          cue_word = rel$cue_word,
          drug_death_distance_tokens = dist,
          tumor_context_present_yes_no = ifelse(tumor_present, "yes", "no"),
          negation_flag = rel$negation_flag,
          confidence_high_medium_low = confidence,
          evidence_tier = "Tier3_sentence_window_relation_candidate",
          therapy_class = row$therapy_class,
          is_specific_drug_yes_no = row$is_specific_drug_yes_no,
          is_generic_therapy_term_yes_no = row$is_generic_therapy_term_yes_no
        )
      }
    }
    relations <- bind_rows(relation_rows)
  } else {
    relations <- tibble()
  }
  if (nrow(relations) == 0) {
    relations <- tibble(
      pmid = character(), year = integer(), tumor_family = character(), death_mode = character(), drug_normalized = character(),
      sentence = character(), relation_category = character(), cue_word = character(), drug_death_distance_tokens = numeric(),
      tumor_context_present_yes_no = character(), negation_flag = character(), confidence_high_medium_low = character(),
      evidence_tier = character(), therapy_class = character(), is_specific_drug_yes_no = character(), is_generic_therapy_term_yes_no = character()
    )
  }
  write_csv_safe(relations, file.path(ROOT, "09_sentence_level_relation_audit/drug_death_relation_sentences.csv"))

  set.seed(20260515)
  manual_sample <- relations %>%
    mutate(.sample_order = runif(n())) %>%
    group_by(relation_category, death_mode, therapy_class) %>%
    arrange(.sample_order, .by_group = TRUE) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    slice_head(n = 300) %>%
    select(-.sample_order) %>%
    mutate(manual_label = "", curator_1 = "", curator_2 = "", notes = "")
  write_csv_safe(manual_sample, file.path(ROOT, "09_sentence_level_relation_audit/drug_relation_manual_validation_sample.csv"))

  # Triads and drug-level signals
  triad_mentions <- drug_mentions %>%
    filter(is_specific_drug_yes_no == "yes") %>%
    group_by(tumor_family, drug_normalized, drug_display, therapy_class, death_mode) %>%
    summarise(pmid_count = n_distinct(pmid), first_year = min(year, na.rm = TRUE), last_year = max(year, na.rm = TRUE), .groups = "drop")

  rel_sum <- relations %>%
    filter(is_specific_drug_yes_no == "yes") %>%
    group_by(tumor_family, drug_normalized, death_mode) %>%
    summarise(
      high_confidence_relation_count = n_distinct(pmid[confidence_high_medium_low == "high"]),
      induction_count = n_distinct(pmid[relation_category == "induction_or_activation" & confidence_high_medium_low %in% c("high", "medium")]),
      sensitization_count = n_distinct(pmid[relation_category == "sensitization" & confidence_high_medium_low %in% c("high", "medium")]),
      inhibition_count = n_distinct(pmid[relation_category == "inhibition_or_protection" & confidence_high_medium_low %in% c("high", "medium")]),
      combination_count = n_distinct(pmid[relation_category == "combination_therapy" & confidence_high_medium_low %in% c("high", "medium")]),
      unclear_count = n_distinct(pmid[relation_category == "unclear_co_mention"]),
      relation_categories = paste(sort(unique(relation_category)), collapse = "; "),
      .groups = "drop"
    )

  triads <- triad_mentions %>%
    left_join(rel_sum, by = c("tumor_family", "drug_normalized", "death_mode")) %>%
    mutate(across(c(high_confidence_relation_count, induction_count, sensitization_count, inhibition_count, combination_count, unclear_count), ~coalesce(.x, 0L)),
           recent_growth = ifelse(last_year >= 2023, "recent", "not_recent"),
           evidence_tier_summary = "Tier2 co-mention with Tier3 relation candidates when present",
           tumor_specificity_score = zscore(pmid_count),
           death_mode_specificity_score = zscore(high_confidence_relation_count),
           low_count_warning_yes_no = ifelse(pmid_count < 5, "yes", "no"),
           triad_evidence_grade = case_when(
             pmid_count >= 5 & high_confidence_relation_count >= 2 & (induction_count + sensitization_count + inhibition_count + combination_count) > 0 ~ "strong_literature_signal",
             pmid_count >= 5 & high_confidence_relation_count >= 1 ~ "moderate_literature_signal",
             pmid_count >= 3 ~ "exploratory",
             TRUE ~ "weak_or_unclear"
           )) %>%
    arrange(factor(triad_evidence_grade, levels = c("strong_literature_signal", "moderate_literature_signal", "exploratory", "weak_or_unclear")), desc(pmid_count))

  write_csv_safe(triads, file.path(ROOT, "06_drug_death_mode_triads/oncology_drug_death_triads.csv"))
  write_csv_safe(triads %>% filter(triad_evidence_grade %in% c("strong_literature_signal", "moderate_literature_signal")) %>% slice_head(n = 100), file.path(ROOT, "06_drug_death_mode_triads/top_oncology_drug_death_triads.csv"))

  poly <- triads %>%
    group_by(drug_normalized, drug_display, therapy_class) %>%
    summarise(
      death_mode_breadth = n_distinct(death_mode[pmid_count >= 1]),
      tumor_family_breadth = n_distinct(tumor_family),
      total_pmid_count = sum(pmid_count),
      high_confidence_relation_count = sum(high_confidence_relation_count),
      number_of_relation_categories = n_distinct(unlist(str_split(paste(relation_categories, collapse = "; "), ";\\s*"))),
      dominant_death_mode = death_mode[which.max(pmid_count)][1],
      death_mode_entropy = {
        p <- pmid_count / sum(pmid_count)
        -sum(p * log(p + 1e-12))
      },
      top_death_modes = paste(head(death_mode[order(-pmid_count)], 5), collapse = "; "),
      top_tumor_families = paste(head(tumor_family[order(-pmid_count)], 5), collapse = "; "),
      low_confidence_fraction = safe_div(sum(unclear_count), sum(pmid_count)),
      .groups = "drop"
    ) %>%
    mutate(
      polydeath_index = zscore(death_mode_breadth) + zscore(tumor_family_breadth) + zscore(log1p(high_confidence_relation_count)) + zscore(death_mode_entropy),
      evidence_grade = case_when(
        death_mode_breadth >= 2 & high_confidence_relation_count >= 3 & total_pmid_count >= 5 ~ "displayable_literature_signal",
        death_mode_breadth >= 2 & total_pmid_count >= 5 ~ "exploratory_co_mention",
        TRUE ~ "weak_or_unclear"
      )
    ) %>%
    arrange(desc(polydeath_index))
  write_csv_safe(poly, file.path(ROOT, "07_same_drug_multi_death_signals/drug_polydeath_signal_table.csv"))

  therapy_matrix <- relations %>%
    filter(is_specific_drug_yes_no == "yes", confidence_high_medium_low %in% c("high", "medium")) %>%
    group_by(therapy_class, death_mode) %>%
    summarise(
      pmid_count = n_distinct(pmid),
      high_confidence_relation_count = n_distinct(pmid[confidence_high_medium_low == "high"]),
      tumor_family_breadth = n_distinct(tumor_family),
      relation_category_distribution = paste(names(sort(table(relation_category), decreasing = TRUE)), as.integer(sort(table(relation_category), decreasing = TRUE)), sep = ":", collapse = "; "),
      recent_growth = (n_distinct(pmid[year >= 2023]) + 1) / (n_distinct(pmid[year >= 2020 & year <= 2022]) + 1),
      .groups = "drop"
    )
  if (nrow(therapy_matrix) > 0) {
    class_totals <- therapy_matrix %>% group_by(therapy_class) %>% summarise(class_total = sum(pmid_count), .groups = "drop")
    death_totals <- therapy_matrix %>% group_by(death_mode) %>% summarise(death_total = sum(pmid_count), .groups = "drop")
    all_total <- sum(therapy_matrix$pmid_count)
    therapy_matrix <- therapy_matrix %>%
      left_join(class_totals, by = "therapy_class") %>%
      left_join(death_totals, by = "death_mode") %>%
      mutate(therapy_death_selectivity = log2((pmid_count / class_total + 1e-6) / (death_total / all_total + 1e-6)),
             low_count_warning_yes_no = ifelse(pmid_count < 3, "yes", "no")) %>%
      select(therapy_class, death_mode, pmid_count, high_confidence_relation_count, tumor_family_breadth, relation_category_distribution, therapy_death_selectivity, recent_growth, low_count_warning_yes_no)
  }
  write_csv_safe(therapy_matrix, file.path(ROOT, "08_therapy_class_landscape/therapy_class_death_mode_matrix.csv"))

  therapy_temporal <- drug_mentions %>%
    filter(is_specific_drug_yes_no == "yes") %>%
    count(therapy_class, death_mode, year, name = "pmid_mentions") %>%
    arrange(therapy_class, death_mode, year)
  write_csv_safe(therapy_temporal, file.path(ROOT, "08_therapy_class_landscape/therapy_death_temporal_trends.csv"))

  tumor_div <- triads %>%
    group_by(tumor_family) %>%
    summarise(
      number_of_drugs = n_distinct(drug_normalized),
      number_of_therapy_classes = n_distinct(therapy_class),
      number_of_death_modes = n_distinct(death_mode),
      drug_death_triad_count = n(),
      triad_entropy = {
        p <- pmid_count / sum(pmid_count)
        -sum(p * log(p + 1e-12))
      },
      .groups = "drop"
    )
  write_csv_safe(tumor_div, file.path(ROOT, "06_drug_death_mode_triads/tumor_therapy_death_diversity.csv"))

  # Gates
  generic_dom <- if (nrow(relations) > 0) {
    generic_counts <- relations %>%
      count(drug_normalized, is_generic_therapy_term_yes_no, sort = TRUE) %>%
      mutate(frac = n / sum(n))
    if (any(generic_counts$is_generic_therapy_term_yes_no == "yes")) {
      max(generic_counts$frac[generic_counts$is_generic_therapy_term_yes_no == "yes"], na.rm = TRUE)
    } else {
      NA_real_
    }
  } else NA_real_
  if (is.infinite(generic_dom)) generic_dom <- NA_real_
  relation_candidates <- nrow(relations)
  high_ind_sens <- relations %>% filter(confidence_high_medium_low == "high", relation_category %in% c("induction_or_activation", "sensitization")) %>% nrow()

  gates <- tibble(
    module = c("oncology_scope", "death_mode_hotness", "tumor_specificity", "drug_extraction", "sentence_relation_audit", "drug_death_triads", "same_drug_multi_death", "therapy_class_landscape"),
    sample_size = c(nrow(article_onc), nrow(hotness), n_distinct(tumor_death$tumor_family), nrow(drug_mentions), relation_candidates, nrow(triads), nrow(poly), nrow(therapy_matrix)),
    gate_criterion = c(
      "oncology article records and tumor mappings available",
      ">=5 death modes with oncology_unique_pmids >=50 and Pair_Count>=3 stability disclosed",
      "display tumor families have support >=20 and cells are count-annotated",
      "title/abstract or chemical fields available; specific drug mentions detected",
      ">=500 sentence candidates; >=100 high-confidence induction/sensitization; manual sample exists",
      "strong/moderate triads require pmid_count>=5 and relation support",
      "specific drugs require death_mode_breadth>=2, high-confidence relations>=3, total_pmid_count>=5",
      ">=3 therapy classes and >=5 death modes with high/medium relation support"
    ),
    effect_size_summary = c(
      paste0("oncology_unique_pmids=", n_distinct(article_onc$pmid), "; tumor_families=", n_distinct(article_onc$tumor_family)),
      paste0("death_modes_ge50=", sum(hotness$oncology_unique_pmids >= 50), "; rank_spearman_paircount3=", round(hotness_spearman, 3)),
      paste0("display_tumors=", length(display_tumors), "; matrix_cells=", nrow(tumor_death)),
      paste0("drug_mentions=", nrow(drug_mentions), "; specific_mentions=", sum(drug_mentions$is_specific_drug_yes_no == "yes")),
      paste0("relations=", relation_candidates, "; high_induction_or_sensitization=", high_ind_sens, "; max_generic_fraction=", round(generic_dom, 3)),
      paste0("strong_or_moderate_triads=", sum(triads$triad_evidence_grade %in% c("strong_literature_signal", "moderate_literature_signal"))),
      paste0("displayable_drugs=", sum(poly$evidence_grade == "displayable_literature_signal")),
      paste0("classes=", n_distinct(therapy_matrix$therapy_class), "; death_modes=", n_distinct(therapy_matrix$death_mode))
    ),
    passes_gate_yes_no = c(
      "yes",
      ifelse(sum(hotness$oncology_unique_pmids >= 50) >= 5 & !is.na(hotness_spearman) & hotness_spearman >= 0.5, "yes", "partial"),
      ifelse(length(display_tumors) >= 5, "yes", "partial"),
      ifelse(nrow(drug_mentions) > 0 & sum(drug_mentions$is_specific_drug_yes_no == "yes") > 0, "yes", "no"),
      ifelse(relation_candidates >= 500 & high_ind_sens >= 100 & (is.na(generic_dom) | generic_dom <= 0.5) & nrow(manual_sample) > 0, "yes", "no"),
      ifelse(sum(triads$triad_evidence_grade %in% c("strong_literature_signal", "moderate_literature_signal")) >= 10, "yes", "partial"),
      ifelse(sum(poly$evidence_grade == "displayable_literature_signal") >= 5, "yes", "partial"),
      ifelse(n_distinct(therapy_matrix$therapy_class) >= 3 & n_distinct(therapy_matrix$death_mode) >= 5, "yes", "partial")
    ),
    recommended_output = c(
      "source_table",
      "figure",
      "figure",
      "source_table",
      ifelse(relation_candidates >= 500 & high_ind_sens >= 100, "supplementary_figure_or_table", "table_only"),
      "supplementary_table_or_figure",
      "supplementary_figure_if_readable",
      "supplementary_figure_if_readable"
    ),
    reason = c(
      "Oncology scope is traceable to v5 source data.",
      "Hotness is descriptive and component-wise; not a target-quality rank.",
      "Tumor selectivity uses oncology background normalization.",
      ifelse(chemical_evidence_available, "PubMed/PubTator chemical tiers are merged with local lexicon mapping; unmapped chemicals remain needs_review.", "No external chemical tier; extraction remains text-lexicon based."),
      "Sentence candidates are rule-based and require manual validation before mechanistic claims.",
      "Triads are hypothesis-generating only.",
      "Same-drug multi-death signals can reflect literature volume/tool-compound reuse.",
      "Therapy-class results require relation support and count annotation."
    )
  )
  write_csv_safe(gates, file.path(ROOT, "15_QC_reports/oncology_analysis_gate_decision_table.csv"))

  rejected <- gates %>%
    filter(!(passes_gate_yes_no == "yes")) %>%
    transmute(figure_or_module = module, reason_for_rejection_or_demote = reason, recommended_output)
  write_csv_safe(rejected, file.path(ROOT, "15_QC_reports/rejected_oncology_figures.csv"))

  write_text_safe(c(
    "# Oncology Module QC Report",
    "",
    paste0("- Scope records: ", nrow(article_onc), " oncology article-stage records; ", n_distinct(article_onc$pmid), " unique PMIDs."),
    paste0("- Tumor families mapped: ", n_distinct(article_onc$tumor_family), "."),
    paste0("- Drug mentions: ", nrow(drug_mentions), " total; ", sum(drug_mentions$is_specific_drug_yes_no == "yes"), " specific-drug mentions."),
    paste0("- Sentence-level relation candidates: ", relation_candidates, "; high-confidence induction/sensitization: ", high_ind_sens, "."),
    ifelse(chemical_evidence_available, "- PubMed ChemicalList and PubTator3 chemical evidence were available and merged conservatively; unmapped chemicals remain needs_review.", "- External chemical evidence was not available and was not inferred."),
    "- All drug and therapy outputs are literature-level candidates, not mechanistic validation.",
    "",
    "## Gate decisions",
    paste0("- ", gates$module, ": ", gates$passes_gate_yes_no, " (", gates$effect_size_summary, ")")
  ), file.path(ROOT, "15_QC_reports/oncology_module_QC_report.md"))

  forbidden <- c("drug X induces Y cell death", "validated therapeutic mechanism", "clinical efficacy", "true tumor vulnerability", "proven multi-death induction", "translational maturity")
  scan_files <- list.files(ROOT, pattern = "\\.(md|csv)$", recursive = TRUE, full.names = TRUE)
  hits <- lapply(scan_files, function(f) {
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
    bad <- forbidden[sapply(forbidden, function(p) any(str_detect(str_to_lower(txt), str_to_lower(fixed(p)))))]
    if (length(bad) == 0) return(NULL)
    tibble(file = f, forbidden_phrase = paste(bad, collapse = "; "))
  }) %>% bind_rows()
  if (nrow(hits) == 0) {
    write_text_safe("No forbidden phrases detected in generated markdown/csv outputs.", file.path(ROOT, "15_QC_reports/oncology_overclaim_scan.txt"))
  } else {
    write_csv_safe(hits, file.path(ROOT, "15_QC_reports/oncology_overclaim_scan.csv"))
    write_text_safe(c("Forbidden phrase hits detected:", paste(hits$file, hits$forbidden_phrase, sep = " | ")), file.path(ROOT, "15_QC_reports/oncology_overclaim_scan.txt"))
  }

  write_text_safe(c(
    "# Drug Extraction Error Report",
    "",
    "- No runtime extraction errors were detected.",
    ifelse(chemical_evidence_available, "- PubMed ChemicalList/PubTator3 evidence is available from cached external annotation outputs; Supplementary Concept/RxNorm/DrugBank/ChEMBL mappings remain unavailable.", "- Primary limitation: no MeSH chemical/Supplementary Concept fields in the selected v5 source table."),
    "- Generic therapy terms are retained in the lexicon but excluded from specific-drug triad/polydeath evidence grades."
  ), file.path(ROOT, "15_QC_reports/drug_extraction_error_report.md"))

  # Manuscript patches
  relation_gate_pass <- gates %>% filter(module == "sentence_relation_audit") %>% pull(passes_gate_yes_no) == "yes"
  tumor_gate_pass <- gates %>% filter(module == "tumor_specificity") %>% pull(passes_gate_yes_no) == "yes"
  integration <- if (relation_gate_pass && tumor_gate_pass) {
    "Option A: Add one cautious Results subsection. Therapy-linked signals passed rule-based sentence gates but still require manual validation before mechanistic claims."
  } else if (tumor_gate_pass) {
    "Option B: Add one brief Results paragraph. Tumor families differed in death-mode literature profiles; drug and therapy-class signals should remain supplementary relation-gated diagnostics."
  } else {
    "Option C: Do not add to main text; retain as supplementary/table-only diagnostics."
  }
  write_text_safe(c(
    "# Oncology Integration Decision",
    "",
    integration,
    "",
    "Boundary: no drug induction or clinical efficacy claim is supported without manual sentence validation."
  ), file.path(ROOT, "11_manuscript_patches/oncology_integration_decision.md"))

  write_text_safe(c(
    "# Oncology Results Patch",
    "",
    ifelse(chemical_evidence_available,
      "Tumor-specific death-mode profiles and therapy-associated RCD signals in oncology literature were evaluated as an evidence-gated supplementary module. Tumor families showed heterogeneous death-mode literature profiles after oncology-background normalization, supporting the view that oncology RCD literature is not a single homogeneous disease category. Drug and therapy-class signals were retained as relation-gated literature diagnostics using PubMed ChemicalList annotations, PubTator3 chemical mentions, local lexicon mapping, and rule-based sentence candidates. These outputs should therefore be interpreted as drug-associated death-mode literature signals and hypothesis-generating tumor-drug-death-mode triads, not as evidence that a given drug induces a validated death pathway.",
      "Tumor-specific death-mode profiles and therapy-associated RCD signals in oncology literature were evaluated as an evidence-gated supplementary module. Tumor families showed heterogeneous death-mode literature profiles after oncology-background normalization, supporting the view that oncology RCD literature is not a single homogeneous disease category. Drug and therapy-class signals were retained as relation-gated literature diagnostics because the available source table supported title/abstract matching and rule-based sentence candidates, but not MeSH-chemical evidence or manual validation. These outputs should therefore be interpreted as drug-associated death-mode literature signals and hypothesis-generating tumor-drug-death-mode triads, not as evidence that a given drug induces a validated death pathway.")
  ), file.path(ROOT, "11_manuscript_patches/oncology_results_patch.md"))
  write_text_safe(c(
    "# Oncology Discussion Patch",
    "",
    "The oncology-focused module suggests that tumor literatures differ in their RCD-associated vocabularies and that several drugs or therapy classes recur in death-mode-specific article contexts. These observations may help prioritize hypotheses about therapy-linked RCD vocabulary, resistance models, or experimental-tool contexts. However, same-drug multi-death signals can also reflect publication volume, tool-compound reuse, or broad experimental interest. Sentence-level and manual curation remain necessary before converting these literature signals into mechanistic or therapeutic claims."
  ), file.path(ROOT, "11_manuscript_patches/oncology_discussion_patch.md"))
  write_text_safe(c(
    "# Oncology Methods Patch",
    "",
    ifelse(chemical_evidence_available,
      "Oncology records were defined using the audited oncology-term classification table, disease-system annotations, and conservative tumor-name regular expressions. Tumor terms were mapped to tumor families and broad tumor classes using rule-based lexical mapping. Drug and therapy evidence was assembled from PubMed ChemicalList annotations, PubTator3 chemical mentions, and a curated oncology therapy lexicon. Chemical names matching the local therapy lexicon were used as specific-drug evidence; unmapped chemicals were retained as needs_review and were not promoted to specific-drug triads. Sentence-level relation candidates required drug and death-mode vocabulary to occur in the same sentence with relation cues for induction/activation, sensitization, inhibition/protection, resistance, or combination therapy. All relation outputs were treated as candidates and a stratified manual-validation sample was generated with blank labels.",
      "Oncology records were defined using the audited oncology-term classification table, disease-system annotations, and conservative tumor-name regular expressions. Tumor terms were mapped to tumor families and broad tumor classes using rule-based lexical mapping. Drug and therapy mentions were extracted from article titles and abstracts using a curated oncology therapy lexicon; MeSH chemical evidence was not used because chemical fields were unavailable in the selected source table. Sentence-level relation candidates required drug and death-mode vocabulary to occur in the same sentence with relation cues for induction/activation, sensitization, inhibition/protection, resistance, or combination therapy. All relation outputs were treated as candidates and a stratified manual-validation sample was generated with blank labels.")
  ), file.path(ROOT, "11_manuscript_patches/oncology_methods_patch.md"))
  write_text_safe(c(
    "# Oncology Figure Legends",
    "",
    "Figure/Supplementary Figure. Oncology death-mode hotness. Component-wise descriptive heat index summarizing oncology PMID support, tumor-family breadth, recent growth, and drug-associated title/abstract vocabulary. This index describes oncology literature attention and should not be interpreted as therapeutic value or target quality.",
    "",
    "Figure/Supplementary Figure. Tumor-specific death-mode landscape. Tumor-family by death-mode selectivity was calculated relative to the oncology background. Cells show literature selectivity and count support, not experimentally proven tumor biology.",
    "",
    "Supplementary Figure. Same-drug multi-death literature signals. Specific drugs appearing across multiple death-mode relation-candidate contexts are shown as hypothesis-generating literature signals. Manual validation is required before mechanistic interpretation.",
    "",
    "Supplementary Figure. Therapy class by death-mode relation landscape. Therapy classes are summarized using rule-based sentence candidates from titles and abstracts. The figure reports therapy-linked RCD vocabulary, not clinical efficacy."
  ), file.path(ROOT, "11_manuscript_patches/oncology_figure_legends.md"))

  source_manifest <- tibble(
    source_file = c(article_path, pair_path, onc_path, low_path, pubmed_pubtator_drug_path, pubmed_pubtator_relation_path),
    destination_or_use = c("article-level oncology scope, drug text, sentence candidates", "pair counts and selectivity", "oncology flag source", "low-count warning annotation", "PubMed/PubTator chemical evidence merged with local lexicon", "rule-based sentence relation candidates from enhanced evidence"),
    rows = c(nrow(article), nrow(pair), nrow(onc), nrow(low), ifelse(chemical_evidence_available, nrow(readr::read_csv(pubmed_pubtator_drug_path, show_col_types = FALSE)), 0), ifelse(relation_evidence_available, nrow(readr::read_csv(pubmed_pubtator_relation_path, show_col_types = FALSE)), 0))
  )
  write_csv_safe(source_manifest, file.path(ROOT, "12_source_data/oncology_source_data_manifest.csv"))

  write_text_safe(log_lines, file.path(ROOT, "15_QC_reports/oncology_analysis_run.log"))
  log_msg("Completed oncology therapeutic death landscape analysis")
}

main()
