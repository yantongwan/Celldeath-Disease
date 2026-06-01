---
output:
  word_document: default
  html_document: default
---
# Nature-style methods for all selected CellDeath_Figure PNGs

This file summarizes the methods used for the selected PNG figures in the `CellDeath_Figure` repository: `Figure 1.png` to `Figure 5.png` and `Figure S1.png` to `Figure S7.png`. macOS AppleDouble sidecar files named `._Figure*.png` were excluded. The methods are organized into five modules, following a Nature-style emphasis on transparent input data, prespecified filtering, metric definitions, reproducibility and interpretation boundaries.

All figure-level signals are literature-level signals derived from PubMed records, PubMed/PubTator annotations, title/abstract/keyword co-mentions, or structured intermediate tables. They should not be interpreted as causal biological evidence, clinical efficacy evidence, or experimental pathway validation.

## Figure-to-method mapping

| Figure PNG | Main method module | Main data layer |
|---|---|---|
| `Figure 1.png` | Module 1 | Temporal and disease-pair PubMed atlas |
| `Figure 2.png` | Modules 1 and 2 | Literature breadth, disease hubs, country/journal context and radar summaries |
| `Figure 3.png` | Module 2 | Evidence-stage audit |
| `Figure 4.png` | Module 3 | Gene/protein co-mention layer |
| `Figure 5.png` | Module 4 | Oncology therapeutic and approved-drug literature landscape |
| `Figure S1.png` | Module 1 | Normalized temporal trajectories |
| `Figure S2.png` | Modules 1 and 2 | Low-count pair and threshold-robustness analyses |
| `Figure S3.png` | Module 2 | Disease-system, non-oncology hub and publication ecosystem summaries |
| `Figure S4.png` | Module 2 | Country and journal specialization heatmaps |
| `Figure S5.png` | Module 3 | Gene matching workflow, STRING/network asset and oncology fraction scatter |
| `Figure S6.png` | Module 3 | VOSviewer disease-gene networks |
| `Figure S7.png` | Module 4 | Refined approved-drug landscape |

## Module 1. PubMed atlas construction, death-mode screening and temporal metrics

### Death-mode screening and final figure scope

An upstream exploratory term screen evaluated a broader set of regulated cell death concepts. The public figure package is restricted to the final eight selected concepts used in the manuscript figures.

The selected PNG set uses the final eight-concept selected death-mode 2000-2025 figure scope:

`Ferroptosis`, `Pyroptosis`, `NETosis`, `Necroptosis`, `Immunogenic cell death`, `Cuproptosis`, `PANoptosis` and `Disulfidptosis`.

The final atlas was filtered to publication years 2000-2025, and PubMed was not re-queried during final figure regeneration. After filtering, the atlas contained 8 death-mode concepts, 1,621 disease-system-mapped disease terms with at least one filtered pair, 4,793 positive death-mode disease pairs, 51,061 pair-PMID article-stage rows and 20,645 unique PubMed records.

Primary source files:

- `table/00_manifest/selected_death_modes_manifest.csv`
- `raw/01_core_pubmed_atlas/article_stage_records_selected_death_modes_2000_2025.csv`
- `raw/01_core_pubmed_atlas/original_pair_pmids_specific_v4_balanced_all.csv`

### Death-mode x disease pair construction

For each retained death mode `m` and disease term `d`, records were represented as death-mode-disease-PMID links. The pair support metric was:

```text
Pair_Count(m,d) = number of unique PMIDs linked to death mode m and disease term d
```

`Pair_PMID_Links` is the number of row-level links and can exceed the number of unique PMIDs when the same article contributes to more than one mapped context. Unless otherwise stated, plotted article counts use unique PMIDs.

For each death mode:

```text
Death_Article_N(m) = number of unique PMIDs linked to death mode m
Disease_Breadth(m) = number of disease terms d with Pair_Count(m,d) > 0
Total_Pair_Count(m) = sum_d Pair_Count(m,d)
Positive_Pairs(m) = number of retained disease pairs for death mode m
```

For death-mode disease selectivity:

```text
Share_Within_Death_Cutoff(m,d) = Pair_Count(m,d) / sum_d Pair_Count(m,d)
Disease_Background_Share(d) = sum_m Pair_Count(m,d) / sum_{m,d} Pair_Count(m,d)
Log2_Selectivity_Cutoff(m,d) =
  log2((Share_Within_Death_Cutoff(m,d) + 1e-9) /
       (Disease_Background_Share(d) + 1e-9))
```

This selectivity score is a literature concentration metric, not a disease-specific biological effect size.

### Temporal metrics

Annual publication curves were computed from unique PMIDs per death mode and publication year:

```text
Annual_Publication_Count(m,y) = number of unique PMIDs linked to death mode m in year y
```

For each death mode:

```text
Introduction_Year(m) = earliest year y with Annual_Publication_Count(m,y) > 0
Peak_Year(m) = year y with the maximum Annual_Publication_Count(m,y)
Takeoff_Year(m) = earliest year y with
  Annual_Publication_Count(m,y) >= max(3, 0.1 * max_y Annual_Publication_Count(m,y))
```

For normalized temporal plots:

```text
Normalized_Annual_Count(m,y) =
  Annual_Publication_Count(m,y) / max_y Annual_Publication_Count(m,y)
```

This normalization compares adoption trajectories rather than absolute publication volume.

### Low-count and threshold analyses

Low-count panels summarized the abundance of death-mode-disease pairs with limited PMID support. Threshold-robustness analyses recomputed metrics after filtering pairs by `Pair_Count >= 1`, `>= 3`, `>= 5` and `>= 10`. Rank stability was assessed by Spearman correlation, Kendall correlation and top-N Jaccard overlap, where:

```text
TopN_Jaccard(A,B) = |TopN_A intersect TopN_B| / |TopN_A union TopN_B|
```

The low-count plots are used to disclose sensitivity to sparse literature pairs, not to exclude sparse emerging disease signals by default.

## Module 2. Literature hubs, evidence-stage audit and publication ecosystem metrics

### Disease hub score

Disease hub scores were computed at the disease-term level from pair support, death-mode breadth and two pair-level literature scores. Pair-level mechanism and specificity values were aggregated as `Pair_Count`-weighted means. For disease term `d`:

```text
pair_count_d = sum_m Pair_Count(m,d)
death_mode_breadth_d = number of death modes linked to d
mean_cleaned_mcs_proxy_d = weighted mean of Mechanistic_Convergence_Score, weights = Pair_Count
mean_dss_d = weighted mean of Disease_Specificity_Score, weights = Pair_Count
```

The plotted hub score was:

```text
hub_score_d =
  z(log1p(pair_count_d)) +
  z(death_mode_breadth_d) +
  z(mean_cleaned_mcs_proxy_d) +
  z(mean_dss_d)
```

where `z(x) = (x - mean(x)) / sd(x)` within the scored table. If the standard deviation was zero, the z-score component was set to zero. The hub score ranks disease terms as literature hubs and does not imply that the disease is a biological master regulator.

Two hub tables were used: an original hub ranking and a no-neoplasm sensitivity ranking. Oncology terms were identified by disease-system assignment to `Neoplasms` or by term regex including `neoplasm`, `cancer`, `carcinoma`, `tumor/tumour`, `melanoma`, `glioma`, `leukemia`, `lymphoma`, `sarcoma`, `glioblastoma`, `adenocarcinoma`, `hepatocellular` and `myeloma`.

### Evidence-stage audit

Evidence-stage panels were based on PMID-stage records, not unique studies. A single PMID could contribute multiple records when linked to multiple disease or death-mode contexts.

The audited evidence-stage denominator was:

```text
mode_stage_records(m) = number of PMID-stage records for death mode m
stage_share(m,s) = records assigned to stage s for death mode m / mode_stage_records(m)
```

The audit used PubMed PublicationType support to separate text-level clinical vocabulary from publication-type-supported trial evidence. In the final audited table:

- 51,061 PMID-stage records were audited.
- 1,043 original guideline-like assignments were marked as guideline-ineligible in the article corpus.
- 0 valid guideline evidence records were retained.
- 350 records were retained as PublicationType-supported trials.
- 1,926 records from 792 unique PMIDs were reclassified from human observational to preclinical animal when the human signal was only caused by soft terms such as `survival analysis`, `serum`, `plasma` or `biomarker`, no hard human-observational or publication-type support was present, and explicit preclinical-animal support was present.

This audit was designed to prevent clinical-sounding vocabulary from being overinterpreted as mature clinical evidence.

### Country, journal and publication ecosystem metrics

Country and journal summaries used enriched PMID metadata merged to article-stage records. Country signals were derived from affiliation metadata, not author nationality. Journal signals indicate where atlas literature is concentrated, not journal quality.

For country and journal specialization heatmaps, counts were based on distinct PMID-entity-death-mode records. For entity `e`, death mode `m` and `K` death modes:

```text
n(e,m) = unique PMIDs for entity e and death mode m
entity_share(e,m) = (n(e,m) + 0.5) / (sum_m n(e,m) + 0.5*K)
expected_share(m) = (sum_e n(e,m) + 0.5) / (sum_{e,m} n(e,m) + 0.5*K)
log2_RCA(e,m) = log2(entity_share(e,m) / expected_share(m))
```

The full RCA matrix was clipped to `[-4, 4]`; heatmap visualization used a symmetric scale centered at zero, typically displayed over `[-3, 3]`. Positive values indicate a relative concentration of a death mode within the entity compared with the atlas background.

## Module 3. Gene/protein co-mention layer and VOSviewer network extensions

### Gene universe and primary co-mention matching

The gene layer used a human-mouse shared gene universe and title/abstract/keyword co-mention matching. Primary matches required an accepted matching rule such as unambiguous approved symbols/names or high-risk ambiguous symbols with context support. The primary match table contains the fields `pmid`, `year`, `death_mode`, `disease_term`, `disease_system`, `oncology_flag`, `shared_gene_id`, `shared_display_symbol`, `human_symbol`, `mouse_symbol`, `matched_term`, `matched_term_type`, `match_field`, `ambiguity_status`, `primary_inclusion`, `sensitivity_inclusion`, `match_rule`, `context_cue_present_yes_no` and `low_count_pair_warning_yes_no`.

The matching summary for the organized figure package included 640,889 raw matches, 116,674 primary matches, 300,732 sensitivity matches and 5,492 unique primary shared genes across 16,692 unique primary PMIDs. These values describe the co-mention layer used for figure generation and should not be read as direct gene-function evidence.

Primary source files:

- `raw/02_gene_v3_pubmed_context/pair_pmid_text_context.csv`
- `raw/02_gene_v3_pubmed_context/gene_comention_matches_raw.csv`
- `table/10_Figure_S5_gene_workflow_string_scatter/gene_comention_matches_primary.csv`
- `table/10_Figure_S5_gene_workflow_string_scatter/gene_comention_match_summary.csv`

### Gene-level metrics

For a shared gene `g`:

```text
pmid_count_g = number of unique primary-match PMIDs linked to g
death_mode_breadth_g = number of death modes linked to g
disease_term_breadth_g = number of disease terms linked to g
disease_system_breadth_g = number of disease systems linked to g
multi_death_signal_score_g =
  death_mode_breadth_g * log10(pmid_count_g + 1)
```

Broad multi-death genes were ranked by `multi_death_signal_score`, then by death-mode breadth, PMID count and gene symbol. Death-mode-selective genes were defined by `death_mode_breadth = 1`. Gene x death-mode heatmaps used:

```text
heatmap_value(g,m) = log10(PMID_count(g,m) + 1)
```

Oncology-enriched genes were selected when:

```text
oncology_fraction_g >= 0.75
```

Non-oncology-enriched genes were selected when:

```text
non_oncology_fraction_g >= 0.75
```

Mixed genes were those with `0.25 < oncology_fraction_g < 0.75`. The oncology fraction is the fraction of the gene's primary co-mention signal that occurs in oncology-flagged disease contexts.

### STRING and VOSviewer-style networks

STRING and VOSviewer-style panels are visual network summaries of literature-derived co-mention outputs. They do not validate physical protein-protein interactions unless explicitly imported from STRING as a visual asset.

Eight-mode shared gene networks used genes with `death_mode_breadth >= 8`. Gene network nodes were shared genes and edges were weighted by the number of death modes shared by a gene pair:

```text
edge_weight(g1,g2) = number of shared death modes between g1 and g2
```

Disease-gene VOSviewer networks were built separately for non-oncology and oncology disease layers. Disease nodes were retained if they linked to more than five distinct shared genes. Node weights included:

```text
weight<Genes> = number of distinct shared genes linked to the disease node
weight<PMIDs> = number of distinct PMIDs linked to the disease node
weight<DeathModes> = number of distinct death modes linked to the disease node
```

Disease-network edges were weighted by the number of shared genes between disease nodes:

```text
edge_weight(d1,d2) = number of shared genes between disease d1 and disease d2
```

Louvain clustering was applied to weighted graph edges, and layouts were generated with Fruchterman-Reingold layouts before VOSviewer-compatible map and network files were exported. The non-oncology VOSviewer panel was further filtered to nodes with `weight<PMIDs> > 30`, reducing the non-oncology disease network from 862 to 155 nodes and from 308,485 to 11,935 edges, followed by reclustering into five clusters. The oncology disease network contained 156 nodes and 10,396 edges in the exported summary.

## Module 4. Oncology therapeutic and refined approved-drug literature landscape

### Oncology scope and tumor-family mapping

The oncology layer was generated from the selected death-mode 2000-2025 article-stage and pair tables. A record was flagged as oncology if the disease term was classified as `Neoplasms`, if an oncology term regex matched the disease term, or if the upstream oncology classification indicated an oncology context. Tumor-family labels were assigned using rule-based regexes covering breast, liver, lung, colorectal, gastric, pancreatic, ovarian, prostate, renal, melanoma, glioblastoma/CNS, leukemia, lymphoma, myeloma, sarcoma, cervical, bladder, esophageal, head and neck, thyroid, endometrial/uterine, biliary tract and mixed/unspecified neoplasms.

Key oncology metrics were:

```text
oncology_unique_pmids_m = unique oncology PMIDs for death mode m
total_unique_pmids_m = unique PMIDs for death mode m in the full eight-concept atlas
fraction_of_death_mode_literature_in_oncology_m =
  oncology_unique_pmids_m / total_unique_pmids_m
tumor_family_breadth_m = number of tumor families linked to death mode m
recent_growth_3yr_m =
  (PMIDs_2023_2025_m + 1) / (PMIDs_2020_2022_m + 1)
```

### Drug annotation, relation candidates and approved-drug refinement

Drug evidence combined three sources:

1. PubMed EFetch `ChemicalList` annotations.
2. PubTator3 chemical annotations.
3. Local title/abstract drug lexicon fallback.

Network calls were opt-in and cached. The merged drug table retained evidence-source and evidence-tier fields. Sentence-level drug-death relation candidates required a drug mention and a death-mode term in the same sentence. Rule-based relation cues were grouped as `induction_or_activation`, `sensitization`, `inhibition_or_protection`, `resistance_escape`, `combination_therapy` or `unclear_co_mention`. A relation was high confidence when the cue was induction or sensitization, no negation was detected and the estimated drug-death token distance was at most 18; non-unclear non-negated cues were medium confidence; all other candidates were low confidence. These sentence candidates require manual review for mechanism-level claims.

PubMed/PubTator candidates were mapped to a reference universe containing DrugCentral FDA-approved drugs, Drugs@FDA active ingredients and product names, COCONUT natural product candidates and the local curated oncology therapy lexicon. Salt-stripped normalized names were used as a secondary matching layer. For the refined approved-drug figures, only entities mapped to the approved clinical drug reference class were retained.

Common ions, nutrients, metabolites, solvents, excipients and broad physiology terms were excluded from plotted refined approved-drug signals and retained only in an audit table. The refined approved-drug layer contained:

- 14,220 approved clinical drug raw contexts before common-substance exclusion.
- 523 approved clinical drug normalized entities before common-substance exclusion.
- 6,586 excluded common-substance contexts.
- 54 excluded common-substance normalized entities.
- 7,634 refined approved clinical drug contexts used in figures.
- 469 refined approved clinical drug normalized entities used in figures.
- 2,635 refined approved clinical drug PMIDs.

### Oncology heat index and selectivity metrics

The approved-drug-calibrated oncology heat index was:

```text
approved_drug_pmid_fraction_m =
  approved_drug_pmid_count_m / oncology_unique_pmids_m

heat_index_latest_approved_drug_refined_m =
  z(log1p(oncology_unique_pmids_m)) +
  z(tumor_family_breadth_m) +
  z(recent_growth_3yr_m) +
  z(approved_drug_pmid_fraction_m)
```

This index is a component-wise literature index, not a therapeutic target ranking.

Tumor-family death-mode selectivity was computed relative to the oncology background:

```text
tumor_normalized_share(t,m) =
  pair_count(t,m) / sum_m pair_count(t,m)

oncology_background_share(m) =
  death_pair_total_m / sum_{t,m} pair_count(t,m)

log2_tumor_death_selectivity(t,m) =
  log2((tumor_normalized_share(t,m) + 1e-6) /
       (oncology_background_share(m) + 1e-6))
```

Visualization clipped the fill scale to `[-3, 3]`.

Therapy-class death-mode selectivity was:

```text
therapy_death_selectivity(c,m) =
  log2((pmid_count(c,m) / class_total_c + 1e-6) /
       (death_total_m / all_total + 1e-6))
```

The plotted value was clipped to `[-3, 3]`. Points in therapy-class or drug-mode panels indicate induction/sensitization sentence-candidate PMID support when present.

### Drug-mode, tumor-drug-death and same-drug metrics

For a tumor family `t`, drug `r` and death mode `m`:

```text
triad_pmid_count(t,r,m) = number of unique PMIDs for the tumor-drug-death triad
context_count(t,r,m) = number of retained PubMed/PubTator contexts
disease_count(t,r,m) = number of disease terms represented
induction_sensitization_pmids(t,r,m) =
  number of PMIDs with high- or medium-confidence induction/sensitization sentence candidates
```

Approved-drug x death-mode heatmaps used:

```text
approved_drug_heatmap_value(r,m) = log1p(pmid_count(r,m))
```

The same-approved-drug multi-death score used in the refined approved-drug supplementary layer was:

```text
same_drug_multi_death_score_r =
  z(death_mode_breadth_pmid_r) +
  z(log1p(total_pmid_count_r)) +
  z(death_mode_breadth_induction_sensitization_r) +
  z(log1p(total_induction_sensitization_count_r))
```

This score prioritizes approved drugs that recur across multiple death-mode oncology literature contexts, especially when induction/sensitization sentence candidates are present. It does not establish that the drug induces the death mode or has clinical efficacy in the tumor context.

## Module 5. Figure rendering, source traceability and reporting boundaries

### Rendering workflow

Python figures were generated with matplotlib, seaborn, pandas, numpy and networkx. R figures were generated with ggplot2, dplyr, tidyr, patchwork, scales, igraph, clusterProfiler, org.Hs.eg.db and optional ReactomePA where available. VOSviewer-compatible map and network files were exported for network visualization. Publication-grade outputs were saved as PDF, PNG and TIFF, generally at 400-600 dpi, with LZW compression for TIFF exports where used.

The organized figure package keeps the raw, source-data and code layers together:

- `raw/01_core_pubmed_atlas`: PMID-stage atlas records and original pair-PMID links.
- `raw/02_gene_v3_pubmed_context`: gene co-mention raw matches and PMID text context.
- `raw/03_oncology_pubmed_pubtator`: oncology scoped records and PubMed/PubTator drug contexts.
- `table/01_Figure_1_temporal_adoption` to `table/12_Figure_S7_approved_drug_landscape`: figure-specific source tables.
- `code/01_selected_death_modes_v5`: selected death-mode literature atlas scripts.
- `code/02_gene_v3`: gene/protein co-mention and VOSviewer extension scripts.
- `code/03_oncology_therapeutic_landscape`: oncology therapeutic and approved-drug scripts.

### General statistical and reporting conventions

All PMID counts are unique PubMed records unless explicitly labeled as article-stage records, pair-PMID links, context counts or sentence-candidate counts. Composite scores were used for ranking and visualization, not for formal hypothesis testing. z-score components were calculated within the relevant table. Small pseudo-counts were used only to stabilize ratios or log2 selectivity calculations when zero counts were possible.

No randomization, blinding or experimental sample-size calculation was applicable because the study used public bibliographic and annotation data rather than prospectively collected biological samples. The appropriate reporting boundary is therefore bibliometric reproducibility: source records, filtering criteria, derived metrics, code, and figure-specific source data should be provided with the manuscript or as supplementary/source data files.

### Interpretation boundary for the manuscript

Recommended Nature-style wording:

```text
All analyses quantify literature-level associations among regulated cell death terms,
disease terms, genes, tumor contexts and drug mentions. Co-mention, PubMed/PubTator
chemical annotation and sentence-window relation candidates were used to construct
a reproducible atlas of the published literature. These metrics should be interpreted
as evidence-stage, bibliometric and hypothesis-generating signals, not as proof of
causal mechanisms, therapeutic activity or clinical efficacy.
```
