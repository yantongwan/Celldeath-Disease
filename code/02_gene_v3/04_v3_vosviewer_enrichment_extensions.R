#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else Sys.getenv("CELLDEATH_ATLAS_ROOT", unset = getwd())
module_root <- file.path(project_root, "reviewer_revision", "16_gene_protein_death_disease_layer_v3")
out_root <- file.path(module_root, "17_vosviewer_enrichment_extensions")

dirs <- c(
  "00_v2style_figures/pdf", "00_v2style_figures/png", "00_v2style_figures/tiff", "00_v2style_figures/source_data",
  "01_eight_mode_genes", "02_enrichment",
  "03_vosviewer_gene_network", "04_vosviewer_nononcology_disease_network",
  "05_vosviewer_oncology_disease_network", "06_system_broad_top20", "15_QC_reports"
)
for (d in dirs) dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

theme_ext <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#1F2933"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#1F2933"),
      plot.title = element_text(face = "bold", size = base_size + 1.2),
      plot.subtitle = element_text(size = base_size, colour = "#4B5563"),
      legend.position = "right"
    )
}

save_pub <- function(plot, name, width = 7.2, height = 8.6) {
  ggsave(file.path(out_root, "00_v2style_figures/pdf", paste0(name, ".pdf")), plot, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(file.path(out_root, "00_v2style_figures/png", paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 400)
  ggsave(file.path(out_root, "00_v2style_figures/tiff", paste0(name, ".tiff")), plot, width = width, height = height, units = "in", dpi = 400, compression = "lzw")
}

broad <- read_required(file.path(module_root, "07_gene_death_disease_comention", "broad_multi_death_gene_signals_primary.csv"))
death <- read_required(file.path(module_root, "07_gene_death_disease_comention", "gene_death_summary_primary.csv"))
primary <- read_required(file.path(module_root, "06_gene_comention_matching", "gene_comention_matches_primary.csv"))
universe <- read_required(file.path(module_root, "02_shared_gene_universe", "shared_human_mouse_gene_universe.csv"))

broad <- broad %>%
  mutate(
    pmid_count = as.numeric(pmid_count),
    death_mode_breadth = as.numeric(death_mode_breadth),
    disease_term_breadth = as.numeric(disease_term_breadth),
    disease_system_breadth = as.numeric(disease_system_breadth),
    multi_death_signal_score = death_mode_breadth * log10(pmid_count + 1),
    source_species = "human+mouse shared"
  )

# -------------------------------------------------------------------------
# v2-style S21A: Top 40 broad multi-death signals.
# -------------------------------------------------------------------------
top40 <- broad %>%
  arrange(desc(multi_death_signal_score), desc(death_mode_breadth), desc(pmid_count), shared_display_symbol) %>%
  slice_head(n = 40)

mode_levels <- sort(unique(top40$death_mode_breadth))
mode_palette <- setNames(
  grDevices::colorRampPalette(c("#D8ECF4", "#5DA5C8", "#173B57"))(length(mode_levels)),
  as.character(mode_levels)
)
top40 <- top40 %>%
  mutate(
    death_mode_breadth_factor = factor(death_mode_breadth, levels = mode_levels),
    fill_colour = mode_palette[as.character(death_mode_breadth)]
  )
write.csv(top40, file.path(out_root, "00_v2style_figures/source_data", "S21A_v3_top40_broad_multi_death_signals_source_data.csv"), row.names = FALSE)

p_s21a <- ggplot(top40, aes(x = factor(shared_display_symbol, levels = rev(shared_display_symbol)), y = pmid_count, fill = death_mode_breadth_factor)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = paste0(death_mode_breadth, " modes")), hjust = -0.05, size = 2.0, colour = "#263238") +
  coord_flip() +
  scale_fill_manual(
    values = mode_palette,
    breaks = as.character(mode_levels),
    labels = paste0(mode_levels, " modes")
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "A Top 40 broad multi-death signals",
    subtitle = "Ortholog-aware display; ranked by breadth x log10(PMID support)",
    x = NULL, y = "PMIDs", fill = "Death-mode breadth"
  ) +
  theme_ext(base_size = 7)
save_pub(p_s21a, "Supplementary_Figure_S21A_v3_top40_broad_multi_death_signals_v2style")

# -------------------------------------------------------------------------
# Eight-mode gene set and enrichment input.
# -------------------------------------------------------------------------
eight_genes <- broad %>%
  filter(death_mode_breadth >= 8) %>%
  left_join(universe %>% dplyr::select(shared_gene_id, human_entrez_id, human_symbol, human_name, mouse_symbol, mouse_entrez_id, hgnc_id, mgi_id), by = "shared_gene_id", suffix = c("", "_universe")) %>%
  mutate(human_entrez_id = ifelse(is.na(human_entrez_id) | human_entrez_id == "", sub("^shared_human_entrez:", "", shared_gene_id), human_entrez_id)) %>%
  arrange(desc(multi_death_signal_score), desc(pmid_count), shared_display_symbol)

write.csv(eight_genes, file.path(out_root, "01_eight_mode_genes", "v3_primary_eight_death_mode_shared_genes.csv"), row.names = FALSE)
write.csv(eight_genes %>% dplyr::select(shared_display_symbol, human_entrez_id, human_symbol, mouse_symbol, pmid_count, death_mode_breadth, death_modes),
          file.path(out_root, "01_eight_mode_genes", "v3_primary_eight_death_mode_enrichment_gene_list.csv"), row.names = FALSE)

gene_ids <- unique(na.omit(as.character(eight_genes$human_entrez_id)))
gene_ids <- gene_ids[gene_ids != ""]
universe_ids <- unique(na.omit(as.character(universe$human_entrez_id)))
universe_ids <- universe_ids[universe_ids != ""]

safe_enrich <- function(expr, output_name) {
  result <- tryCatch(expr, error = function(e) e)
  path <- file.path(out_root, "02_enrichment", output_name)
  if (inherits(result, "error")) {
    write.csv(data.frame(error = conditionMessage(result)), path, row.names = FALSE)
    return(data.frame(source = output_name, error = conditionMessage(result)))
  }
  df <- as.data.frame(result)
  write.csv(df, path, row.names = FALSE)
  if (nrow(df) == 0) return(data.frame(source = output_name, rows = 0, significant_FDR005 = 0))
  data.frame(source = output_name, rows = nrow(df), significant_FDR005 = sum(df$p.adjust < 0.05, na.rm = TRUE))
}

enrichment_summary <- bind_rows(
  safe_enrich(enrichGO(gene = gene_ids, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE), "eight_mode_human_ortholog_GO_BP.csv"),
  safe_enrich(enrichGO(gene = gene_ids, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE), "eight_mode_human_ortholog_GO_CC.csv"),
  safe_enrich(enrichGO(gene = gene_ids, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE), "eight_mode_human_ortholog_GO_MF.csv"),
  safe_enrich(enrichKEGG(gene = gene_ids, universe = universe_ids, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3), "eight_mode_human_ortholog_KEGG.csv"),
  if (has_reactome) safe_enrich(ReactomePA::enrichPathway(gene = gene_ids, universe = universe_ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, readable = TRUE), "eight_mode_human_ortholog_Reactome.csv") else data.frame(source = "eight_mode_human_ortholog_Reactome.csv", error = "ReactomePA unavailable")
)
write.csv(enrichment_summary, file.path(out_root, "02_enrichment", "eight_mode_enrichment_summary.csv"), row.names = FALSE)

top_terms <- list.files(file.path(out_root, "02_enrichment"), pattern = "eight_mode_human_ortholog_.*\\.csv$", full.names = TRUE) %>%
  lapply(function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"Description" %in% names(df)) return(NULL)
    df$source <- basename(path)
    df
  }) %>%
  bind_rows() %>%
  arrange(p.adjust, pvalue) %>%
  group_by(source) %>%
  slice_head(n = 20) %>%
  ungroup()
write.csv(top_terms, file.path(out_root, "02_enrichment", "eight_mode_top_enrichment_terms_combined.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# VOSviewer gene network: nodes are eight-mode genes; edges weighted by shared death modes.
# -------------------------------------------------------------------------
gene_death_sets <- death %>%
  filter(shared_gene_id %in% eight_genes$shared_gene_id) %>%
  group_by(shared_gene_id, shared_display_symbol) %>%
  summarise(
    death_modes = list(sort(unique(death_mode))),
    .groups = "drop"
  )

node_df <- eight_genes %>%
  transmute(
    id = row_number(),
    shared_gene_id,
    label = shared_display_symbol,
    pmid_count = as.numeric(pmid_count),
    death_mode_breadth = as.numeric(death_mode_breadth),
    disease_term_breadth = as.numeric(disease_term_breadth),
    disease_system_breadth = as.numeric(disease_system_breadth)
  )

edge_rows <- list()
if (nrow(gene_death_sets) >= 2) {
  for (i in seq_len(nrow(gene_death_sets) - 1)) {
    for (j in (i + 1):nrow(gene_death_sets)) {
      shared_modes <- intersect(gene_death_sets$death_modes[[i]], gene_death_sets$death_modes[[j]])
      if (length(shared_modes) > 0) {
        edge_rows[[length(edge_rows) + 1]] <- data.frame(
          source_gene_id = gene_death_sets$shared_gene_id[i],
          target_gene_id = gene_death_sets$shared_gene_id[j],
          source_label = gene_death_sets$shared_display_symbol[i],
          target_label = gene_death_sets$shared_display_symbol[j],
          weight = length(shared_modes),
          shared_death_modes = paste(shared_modes, collapse = "|"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
gene_edges <- bind_rows(edge_rows)
gene_edges <- gene_edges %>%
  left_join(node_df %>% dplyr::select(shared_gene_id, source = id), by = c("source_gene_id" = "shared_gene_id")) %>%
  left_join(node_df %>% dplyr::select(shared_gene_id, target = id), by = c("target_gene_id" = "shared_gene_id")) %>%
  dplyr::select(source, target, weight, source_label, target_label, shared_death_modes)

g_gene <- graph_from_data_frame(gene_edges %>% dplyr::select(source, target, weight), directed = FALSE, vertices = node_df %>% dplyr::select(id, label))
clusters_gene <- cluster_louvain(g_gene, weights = E(g_gene)$weight)
layout_gene <- layout_with_fr(g_gene, weights = E(g_gene)$weight, niter = 1000)
node_df$x <- layout_gene[, 1]
node_df$y <- layout_gene[, 2]
node_df$cluster <- membership(clusters_gene)[as.character(node_df$id)]

gene_map <- node_df %>%
  transmute(
    id, label, x, y, cluster,
    `weight<PMIDs>` = pmid_count,
    `weight<DeathModes>` = death_mode_breadth,
    `weight<DiseaseTerms>` = disease_term_breadth,
    `weight<DiseaseSystems>` = disease_system_breadth,
    description = paste0(label, "; shared human-mouse gene; primary PMID count=", pmid_count)
  )
write_tsv(gene_map, file.path(out_root, "03_vosviewer_gene_network", "v3_eight_mode_gene_network_map.txt"))
write_tsv(gene_edges %>% dplyr::select(source, target, weight), file.path(out_root, "03_vosviewer_gene_network", "v3_eight_mode_gene_network_network.txt"))
write.csv(gene_map, file.path(out_root, "03_vosviewer_gene_network", "v3_eight_mode_gene_network_map.csv"), row.names = FALSE)
write.csv(gene_edges, file.path(out_root, "03_vosviewer_gene_network", "v3_eight_mode_gene_network_network.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Disease networks.
# -------------------------------------------------------------------------
make_disease_network <- function(flag_value, out_dir, label) {
  dat <- primary %>%
    filter(oncology_flag == flag_value) %>%
    distinct(disease_term, disease_system, death_mode, pmid, shared_gene_id, shared_display_symbol)

  disease_stats <- dat %>%
    group_by(disease_term) %>%
    summarise(
      disease_system = paste(sort(unique(disease_system)), collapse = "|"),
      gene_count = n_distinct(shared_gene_id),
      pmid_count = n_distinct(pmid),
      death_mode_breadth = n_distinct(death_mode),
      genes = paste(sort(unique(shared_display_symbol)), collapse = "|"),
      .groups = "drop"
    ) %>%
    filter(gene_count > 5) %>%
    arrange(desc(gene_count), disease_term)

  disease_gene <- dat %>%
    semi_join(disease_stats, by = "disease_term") %>%
    distinct(disease_term, shared_gene_id, shared_display_symbol)

  node <- disease_stats %>%
    mutate(id = row_number()) %>%
    dplyr::select(id, disease_term, disease_system, gene_count, pmid_count, death_mode_breadth, genes)

  split_genes <- split(disease_gene$shared_gene_id, disease_gene$disease_term)
  split_symbols <- split(disease_gene$shared_display_symbol, disease_gene$disease_term)
  edge_list <- list()
  disease_names <- node$disease_term
  if (length(disease_names) >= 2) {
    for (i in seq_len(length(disease_names) - 1)) {
      for (j in (i + 1):length(disease_names)) {
        d1 <- disease_names[i]
        d2 <- disease_names[j]
        shared_ids <- intersect(split_genes[[d1]], split_genes[[d2]])
        if (length(shared_ids) > 0) {
          shared_symbols <- intersect(split_symbols[[d1]], split_symbols[[d2]])
          edge_list[[length(edge_list) + 1]] <- data.frame(
            source_disease = d1,
            target_disease = d2,
            weight = length(shared_ids),
            shared_gene_symbols = paste(sort(unique(shared_symbols)), collapse = "|"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  edges <- bind_rows(edge_list) %>%
    left_join(node %>% dplyr::select(disease_term, source = id), by = c("source_disease" = "disease_term")) %>%
    left_join(node %>% dplyr::select(disease_term, target = id), by = c("target_disease" = "disease_term")) %>%
    dplyr::select(source, target, weight, source_disease, target_disease, shared_gene_symbols)

  if (nrow(edges) > 0 && nrow(node) > 1) {
    graph <- graph_from_data_frame(edges %>% dplyr::select(source, target, weight), directed = FALSE, vertices = node %>% dplyr::select(id, disease_term))
    clusters <- cluster_louvain(graph, weights = E(graph)$weight)
    layout <- layout_with_fr(graph, weights = E(graph)$weight, niter = 1000)
    node$x <- layout[, 1]
    node$y <- layout[, 2]
    node$cluster <- membership(clusters)[as.character(node$id)]
  } else {
    node$x <- seq_len(nrow(node))
    node$y <- 0
    node$cluster <- 1
  }

  map <- node %>%
    transmute(
      id,
      label = disease_term,
      x, y, cluster,
      `weight<Genes>` = gene_count,
      `weight<PMIDs>` = pmid_count,
      `weight<DeathModes>` = death_mode_breadth,
      disease_system,
      description = paste0(label, "; disease genes=", gene_count, "; ", label)
    )

  out_path <- file.path(out_root, out_dir)
  write_tsv(map, file.path(out_path, paste0("v3_", label, "_disease_gene_network_map.txt")))
  write_tsv(edges %>% dplyr::select(source, target, weight), file.path(out_path, paste0("v3_", label, "_disease_gene_network_network.txt")))
  write.csv(map, file.path(out_path, paste0("v3_", label, "_disease_gene_network_map.csv")), row.names = FALSE)
  write.csv(edges, file.path(out_path, paste0("v3_", label, "_disease_gene_network_network.csv")), row.names = FALSE)
  write.csv(disease_stats, file.path(out_path, paste0("v3_", label, "_disease_gene_counts.csv")), row.names = FALSE)
  data.frame(layer = label, nodes = nrow(map), edges = nrow(edges), min_gene_count = 6)
}

network_summary <- bind_rows(
  make_disease_network("False", "04_vosviewer_nononcology_disease_network", "nononcology"),
  make_disease_network("True", "05_vosviewer_oncology_disease_network", "oncology")
)
write.csv(network_summary, file.path(out_root, "15_QC_reports", "v3_vosviewer_network_summary.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Disease-system Top20 broad multi-death genes.
# -------------------------------------------------------------------------
system_gene <- primary %>%
  distinct(disease_system, pmid, death_mode, disease_term, shared_gene_id, shared_display_symbol, human_symbol, mouse_symbol, low_count_pair_warning_yes_no) %>%
  group_by(disease_system, shared_gene_id, shared_display_symbol, human_symbol, mouse_symbol) %>%
  summarise(
    pmid_count = n_distinct(pmid),
    death_mode_breadth = n_distinct(death_mode),
    disease_term_breadth = n_distinct(disease_term),
    low_count_warning_fraction = mean(low_count_pair_warning_yes_no == "yes", na.rm = TRUE),
    death_modes = paste(sort(unique(death_mode)), collapse = "|"),
    .groups = "drop"
  ) %>%
  mutate(
    broad_score = death_mode_breadth * log10(pmid_count + 1)
  ) %>%
  arrange(disease_system, desc(broad_score), desc(death_mode_breadth), desc(pmid_count), shared_display_symbol) %>%
  group_by(disease_system) %>%
  mutate(system_rank = row_number()) %>%
  filter(system_rank <= 20) %>%
  ungroup()
write.csv(system_gene, file.path(out_root, "06_system_broad_top20", "v3_system_disease_broad_multi_death_gene_signals_top20.csv"), row.names = FALSE)

report <- c(
  "# v3 VOSviewer and Enrichment Extension Summary",
  "",
  paste0("- Eight-mode shared genes: ", nrow(eight_genes)),
  paste0("- Human Entrez IDs used for enrichment: ", length(gene_ids)),
  paste0("- VOSviewer eight-mode gene network nodes: ", nrow(gene_map)),
  paste0("- VOSviewer eight-mode gene network edges: ", nrow(gene_edges)),
  paste0("- Non-oncology disease network nodes/edges: ", network_summary$nodes[network_summary$layer == "nononcology"], "/", network_summary$edges[network_summary$layer == "nononcology"]),
  paste0("- Oncology disease network nodes/edges: ", network_summary$nodes[network_summary$layer == "oncology"], "/", network_summary$edges[network_summary$layer == "oncology"]),
  paste0("- Disease-system Top20 rows: ", nrow(system_gene)),
  "",
  "Interpretation boundary: all outputs are title/abstract/keyword co-mention outputs. Enrichment uses the eight-death-mode shared gene set as a hypothesis-generating literature-derived input."
)
writeLines(report, file.path(out_root, "15_QC_reports", "v3_vosviewer_enrichment_extension_summary.md"))
cat(paste(report, collapse = "\n"), "\n")
