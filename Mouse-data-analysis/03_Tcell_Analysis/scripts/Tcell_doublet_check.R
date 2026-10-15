#' ---------------------------------------------------------------------------
#' Module 03 — doublet check across the final αβ T-cell clusters (scDblFinder)
#'
#' Supports the annotation of cluster 6 ("Trm type mixed CD8+ CD4", section 4.8c
#' of Tcell_ab_analysis_reclustering.R): a CD4/CD8-mixed cluster could in
#' principle be CD4–CD8 doublets, so per-cell doublet scores were computed with
#' scDblFinder and summarised per cluster. Cluster 6 should show the LOWEST
#' doublet fraction of any cluster if it is a genuine mixed population.
#'
#' scDblFinder is run on raw counts with each sample (orig.ident) processed
#' separately (`samples =`), because doublets can only form within a capture.
#'
#' Input : ../objects/tcell_TCRab_annotated_final04_clean.rds  (Tcell_remove_MAIT.R)
#' Output: ../objects/cl6_doublet_check.csv
#'
#' Run from 03_Tcell_Analysis/scripts/:  Rscript Tcell_doublet_check.R
#' ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(dplyr)
})

set.seed(123)

obj <- readRDS("../objects/tcell_TCRab_annotated_final04_clean.rds")
obj <- JoinLayers(obj)
cat("Loaded:", ncol(obj), "cells\n")

sce <- SingleCellExperiment(
  assays  = list(counts = GetAssayData(obj, assay = "RNA", layer = "counts")),
  colData = obj@meta.data[, c("orig.ident", "celltype_merged")]
)

sce <- scDblFinder(sce, samples = "orig.ident", BPPARAM = SerialParam(RNGseed = 123))

dbl <- data.frame(
  celltype_merged = sce$celltype_merged,
  is_doublet      = sce$scDblFinder.class == "doublet",
  score           = sce$scDblFinder.score
) %>%
  group_by(celltype_merged) %>%
  summarise(n            = n(),
            n_doublet    = sum(is_doublet),
            pct_doublet  = round(100 * mean(is_doublet), 1),
            median_score = round(median(score), 3),
            .groups = "drop") %>%
  arrange(desc(pct_doublet))

print(as.data.frame(dbl), row.names = FALSE)
write.csv(dbl, "../objects/cl6_doublet_check.csv", row.names = FALSE)
cat("Saved: ../objects/cl6_doublet_check.csv\n")
