# R Package Requirements

Core CRAN packages:

```r
install.packages(c(
  "dplyr",
  "ggplot2",
  "igraph",
  "patchwork",
  "readr",
  "scales",
  "stringr",
  "tidyr"
))
```

Bioconductor packages used by the gene/network extension module:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA"
))
```

`ReactomePA` is optional in the current scripts; the enrichment extension checks whether it is installed before using it.
