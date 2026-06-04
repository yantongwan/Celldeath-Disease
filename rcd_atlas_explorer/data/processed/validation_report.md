# Validation Report

## Source Files Loaded

| logical_table | rows | source_file | column_mapping |
|---|---:|---|---|
| articles | 22076 | article_metadata_emerging_article_nonapoptosis_dedup_2025.csv | `{"abstract": "Abstract", "country": "JournalCountry", "doi": "DOI", "journal": "Journal", "pmid": "PMID", "publication_types": "PublicationTypes", "publication_year": "PubYear", "title": "Title"}` |
| death_mode_disease_articles | 53889 | pair_pmids_emerging_article_nonapoptosis_dedup_2025.csv | `{"death_mode": "Death_Mode", "disease_term": "Disease", "pair_count": "Pair_Count", "pmid": "PMID"}` |
| death_mode_disease_pairs | 5130 | death_disease_pairs_emerging_article_nonapoptosis_positive_2025.csv | `{"death_mode": "Death_Mode", "disease_term": "Disease", "literature_hub_score": "Disease_Specificity_Score", "pair_count": "Pair_Count", "selectivity_score": "Death_Selectivity_Score"}` |
| gene_mentions | 116674 | gene_comention_matches_primary.csv | `{"context_sentence": "match_sentence_or_snippet", "death_mode": "death_mode", "disease_system": "disease_system", "disease_term": "disease_term", "gene_symbol": "shared_display_symbol", "is_oncology": "oncology_flag", "matched_text": "matched_term", "pmid": "pmid", "source_field": "match_field"}` |
| gene_death_disease_articles | 116674 | gene_comention_matches_primary.csv | `{"context_sentence": "match_sentence_or_snippet", "death_mode": "death_mode", "disease_system": "disease_system", "disease_term": "disease_term", "gene_alias": "matched_term", "gene_symbol": "shared_display_symbol", "is_oncology": "oncology_flag", "pmid": "pmid", "source_field": "match_field"}` |
| drug_mentions | 137861 | oncology_drug_mentions.csv | `{"death_mode": "death_mode", "disease_term": "disease_term", "normalized_drug_name": "drug_normalized", "pmid": "pmid", "reference_class_label": "therapy_class", "relation_confidence": "mapping_confidence", "relation_sentence": "mention_sentence_if_available", "source": "mention_source", "tumor_family": "tumor_family"}` |
| drug_mentions | 7634 | pubmed_pubtator_approved_clinical_drug_refined_contexts.csv | `{"death_mode": "death_mode", "disease_term": "disease_term", "normalized_drug_name": "drug_normalized", "pmid": "pmid", "reference_class_label": "reference_class_label", "relation_confidence": "mapping_confidence", "relation_sentence": "mention_sentence_if_available", "source": "mention_source", "tumor_family": "tumor_family"}` |
| drug_death_disease_articles | 137861 | oncology_drug_mentions.csv | `{"death_mode": "death_mode", "disease_term": "disease_term", "drug_name": "drug_display", "normalized_drug_name": "drug_normalized", "pmid": "pmid", "reference_class_label": "therapy_class", "relation_confidence": "mapping_confidence", "relation_sentence": "mention_sentence_if_available", "tumor_family": "tumor_family"}` |
| drug_death_disease_articles | 7634 | pubmed_pubtator_approved_clinical_drug_refined_contexts.csv | `{"death_mode": "death_mode", "disease_term": "disease_term", "drug_name": "drug_display", "normalized_drug_name": "drug_normalized", "pmid": "pmid", "reference_class_label": "reference_class_label", "relation_confidence": "mapping_confidence", "relation_sentence": "mention_sentence_if_available", "tumor_family": "tumor_family"}` |
| disease_dictionary | 1649 | disease_system_mapping.csv | `{"disease_system": "disease_system", "disease_term": "disease_term"}` |

## Final Table Row Counts

| table | rows | columns |
|---|---:|---:|
| articles | 22076 | 14 |
| death_mode_disease_articles | 51061 | 7 |
| death_mode_disease_pairs | 4793 | 11 |
| drug_death_disease_articles | 137861 | 13 |
| drug_mentions | 68753 | 12 |
| gene_death_disease_articles | 116674 | 11 |
| gene_mentions | 47125 | 7 |
| synonyms | 12 | 3 |

## Warnings and Notes

| severity | message |
|---|---|
| warning | articles: 11 publication years outside 2000-2025 |
| info | articles: 0 rows missing title |
| info | articles: 0 rows missing journal |
| info | articles: 0 rows missing publication_year |
| info | gene_death_disease_articles: 46 ASCL4/artifact rows flagged |
| info | drug_death_disease_articles: 0 negated rows retained |

## Interpretation Boundary

Drug-death-disease associations, gene co-mentions and disease-pair counts are literature-level signals. They do not establish biological causality, pathway activation, therapeutic efficacy or clinical validity.
