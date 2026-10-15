#' ---
#' title: "Tcell_ab_analysis_reclustering_v2"
#' author: "Apoorva Sharma"
#' ---
#'
#' **αβ T Cell Reclustering + Annotation + Figures**
#'   - Input: tcell_TCRab_only03.rds (from Tcell_analysis_reclustering.R)
#'   - Reclusters at PC=8, res=0.8 with TCR gene exclusion, the resolution was based on how much granularity/cell types we were expecting.
#'   - Numbered review outputs saved first for annotation review
#'   - Checkpoint RDS saved after clustering (skip recluster on re-runs: set start_from_clustered=TRUE)
#'   - MAIT/Th17 separation (Zbtb16 GMM, cluster 7)
#'   - Tph annotation confirmation by sub-clustering (cluster 2)
#'   - Cluster 6 re-annotated (first we thought it was CD8 Trm) "Trm type mixed CD8+ CD4" (mixed CD4/CD8, section 4.8c)
#'   - Final figures: UMAP, UMAP by condition, pie charts, heatmap, violins
#'
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(readr)
  library(svglite)
  library(pheatmap)
  library(grid)
  library(ggrepel)
  library(mclust)
  library(ggalluvial)
  library(cluster)
  library(RColorBrewer)
})

# Case-insensitive gene name resolver — used throughout for Seurat rowname lookup
resolve_gene <- function(gn, rn) {
  hits <- grep(paste0("^", gn, "$"), rn, ignore.case = TRUE, value = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

## ============================================================================
## CHECKPOINT FLAG — set TRUE to skip reclustering and load from saved RDS
## ============================================================================
start_from_clustered <- FALSE

## ============================================================================
## STAGE 4a — RECLUSTER (skipped if start_from_clustered = TRUE)
## ============================================================================

tcr_pat <- "(^Trav|^Trbv|^Trgv|^Trdv|^Trac$|^Trbc1$|^Trbc2$|^Traj|^Trbj)"

if (!start_from_clustered) {

## ---- 4.1-load ---------------------------------------------------------------
  t_cells_TCRab <- readRDS("../objects/tcell_TCRab_only03.rds")
  t_cells_TCRab <- JoinLayers(t_cells_TCRab)
  DefaultAssay(t_cells_TCRab) <- "RNA"
  cat("Loaded:", ncol(t_cells_TCRab), "αβ TCR cells\n")

## ---- 4.2-normalize ----------------------------------------------------------
  set.seed(123)
  t_cells_TCRab <- NormalizeData(t_cells_TCRab, verbose = FALSE)
  t_cells_TCRab <- FindVariableFeatures(t_cells_TCRab, nfeatures = 3000, verbose = FALSE)

  tcr_in_hvg <- grep(tcr_pat, VariableFeatures(t_cells_TCRab), value = TRUE, ignore.case = TRUE)
  VariableFeatures(t_cells_TCRab) <- setdiff(VariableFeatures(t_cells_TCRab), tcr_in_hvg)
  cat("TCR genes removed from HVGs:", length(tcr_in_hvg), "\n")

  t_cells_TCRab <- ScaleData(t_cells_TCRab,
                              vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)

## ---- 4.3-pca ----------------------------------------------------------------
  t_cells_TCRab <- RunPCA(t_cells_TCRab, npcs = 50, verbose = FALSE, seed.use = 123)

  p_elbow_ab <- ElbowPlot(t_cells_TCRab, ndims = 50) +
    ggtitle("Elbow plot — αβ T cells (TCR genes excluded, 50 PCs)") +
    geom_vline(xintercept = 8, linetype = "dashed", color = "red") +
    annotate("text", x = 9, y = Inf, label = "PC=8 selected",
             hjust = 0, vjust = 1.5, size = 3.5, color = "red") +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank())

  svglite("../figures/elbow_plot_TCRab.svg", width = 8, height = 5)
  print(p_elbow_ab)
  dev.off()
  cat("Saved: ../figures/elbow_plot_TCRab.svg\n")

## ---- 4.4-cluster ------------------------------------------------------------
  pcs_use <- 1:8
  t_cells_TCRab <- FindNeighbors(t_cells_TCRab, reduction = "pca", dims = pcs_use, k.param = 20)
  t_cells_TCRab <- FindClusters(t_cells_TCRab, resolution = 0.8, algorithm = 1,
                                 cluster.name = "RNA_snn_res.0.8")
  t_cells_TCRab <- RunUMAP(t_cells_TCRab, reduction = "pca", dims = pcs_use,
                            n.neighbors = 30, min.dist = 0.15,
                            reduction.name = "umap.pca_TCRab", seed.use = 123)
  Idents(t_cells_TCRab) <- "RNA_snn_res.0.8"
  cat("Clusters (res 0.8):", length(unique(Idents(t_cells_TCRab))), "\n")

## ---- 4.5-markers ------------------------------------------------------------
  markers_r08 <- FindAllMarkers(
    t_cells_TCRab, only.pos = TRUE, test.use = "wilcox",
    logfc.threshold = 0.25, min.pct = 0.1
  )
  top500_r08 <- markers_r08 %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
    dplyr::slice_head(n = 500) %>%
    dplyr::ungroup()
  write.csv(top500_r08, "../objects/markers_TCRab_res0.8_top500.csv", row.names = FALSE)
  cat("Saved: ../objects/markers_TCRab_res0.8_top500.csv\n")

  markers_r08_strict <- FindAllMarkers(
    t_cells_TCRab, only.pos = TRUE, test.use = "wilcox",
    logfc.threshold = 0.5, min.pct = 0.25, min.diff.pct = 0.1
  )
  top500_strict <- markers_r08_strict %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(avg_log2FC, n = 500) %>%
    dplyr::ungroup()
  write.csv(top500_strict, "../objects/markers_TCRab_res0.8_top500_strict.csv", row.names = FALSE)
  cat("Saved: ../objects/markers_TCRab_res0.8_top500_strict.csv\n")

  saveRDS(t_cells_TCRab, "../objects/tcell_TCRab_clustered04a.rds")
  cat("Checkpoint saved: tcell_TCRab_clustered04a.rds\n")

} else {

  t_cells_TCRab <- readRDS("../objects/tcell_TCRab_clustered04a.rds")
  t_cells_TCRab <- JoinLayers(t_cells_TCRab)
  Idents(t_cells_TCRab) <- "RNA_snn_res.0.8"
  cat("Loaded clustered checkpoint:", ncol(t_cells_TCRab), "cells,",
      length(unique(Idents(t_cells_TCRab))), "clusters\n")

}

reduction_used <- "umap.pca_TCRab"

## ============================================================================
## USER CONFIGURATION — Edit annotations and palette after reviewing the
##   numbered figures below. Clusters 1 and 4 both map to "Teff CD8-like" (identical
##   top markers). Clusters 2 and 7 carry provisional labels refined in STAGE 4b:
##   cluster 7 → MAIT contamination removed (GMM); cluster 2 → identity inconclusive,
##   retained as "Tfh/Tph-like" (see section 4.8 conclusion note).
## ============================================================================

# ---- annotations ------------------------------------------------------------
label_short_r08 <- c(
  `0`  = "Th2-like",
  `1`  = "Teff CD8-like",
  `2`  = "Tfh/Tph-like",   # provisional — sub-clustering + marker test inconclusive; retained as Tfh/Tph-like
  `3`  = "Tregs-like",
  `4`  = "Teff CD8-like",
  `5`  = "Naive T-like",
  `6`  = "Trm-like CD8",   # provisional — resolved to "Trm type mixed CD8+ CD4" in section 4.8c (CD4/CD8 mixed)
  `7`  = "Th17/MAIT-like", # provisional — MAIT contamination removed by GMM below
  `8`  = "Trm/Exh-like CD8",
  `9`  = "Teff/Tem-like CD8",
  `10` = "Prolif T"
)

# ---- palette ----------------------------------------------------------------
# Tableau 10 base colors.
# Provisional names (numbered review heatmap) and confirmed final names both
# These colors were set after the annotation.
# mapped to the same color so celltype_merged is covered at every pipeline stage.
pal_tcrab <- c(
  "Teff CD8-like"     = "#4E79A7",
  "Th2-like"          = "#F28E2B",
  "Tfh/Tph-like"      = "#E15759",   # provisional + confirmed (identity inconclusive; see 4.8 notes)
  "Tregs-like"        = "#76B7B2",
  "Trm-like CD8"      = "#59A14F",             # provisional (numbered review); resolved to mixed in section 4.8c
  "Trm type mixed CD8+ CD4" = "#59A14F",       # confirmed final (cluster 6 — mixed CD4/CD8, see section 4.8c)
  "Naive T-like"      = "#EDC948",
  "Th17/MAIT-like"    = "#B07AA1",   # provisional
  "Th17-like"         = "#B07AA1",   # confirmed final
  "Trm/Exh-like CD8"  = "#FF9DA7",
  "Teff/Tem-like CD8" = "#9C755F",
  "Prolif T"          = "#BAB0AC"
)

## ============================================================================
## NUMBERED REVIEW FIGURES — for manual annotation review
## ============================================================================

## ---- cluster-palette --------------------------------------------------------
# Numbered UMAP: one distinct color per cluster (default ggplot2 hue, numerical order).
# No annotation mapping here — clusters 1 and 4 get different colors in these review plots.
cl_ord_num <- as.character(sort(as.integer(unique(as.character(Idents(t_cells_TCRab))))))
cl_pal     <- setNames(scales::hue_pal()(length(cl_ord_num)), cl_ord_num)

## ---- numbered-umap ----------------------------------------------------------
umap_num_emb <- t_cells_TCRab@reductions[["umap.pca_TCRab"]]@cell.embeddings
umap_num_df  <- data.frame(
  UMAP_1    = umap_num_emb[, 1],
  UMAP_2    = umap_num_emb[, 2] * -1,
  cluster   = as.character(Idents(t_cells_TCRab)),
  condition = t_cells_TCRab$condition
)
label_num_df <- umap_num_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap_numbered <- ggplot(umap_num_df, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.4, alpha = 0.7) +
  geom_label(data = label_num_df, aes(label = cluster),
             size = 3.5, fontface = "bold", fill = "white", alpha = 0.85,
             show.legend = FALSE) +
  scale_color_manual(values = cl_pal) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = "none") +
  labs(title = sprintf("αβ T cells — %d clusters (res 0.8, PC=8)", length(cl_ord_num)))

svglite("../figures/UMAP_TCRab_numbered.svg", width = 7, height = 6)
print(p_umap_numbered)
dev.off()
cat("Saved: ../figures/UMAP_TCRab_numbered.svg\n")

## ---- numbered-umap-by-condition ---------------------------------------------
cond_order <- c("BAS", "PD", "REC", "RCU")
cond_panels <- lapply(cond_order, function(cond) {
  fg  <- umap_num_df[umap_num_df$condition == cond, ]
  lbl <- fg %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
  ggplot(fg, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.35, alpha = 0.8) +
    geom_label(data = lbl, aes(UMAP_1, UMAP_2, label = cluster, color = cluster),
               size = 3, fontface = "bold", fill = "white", alpha = 0.85, show.legend = FALSE) +
    scale_color_manual(values = cl_pal) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(), aspect.ratio = 1, legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = cond, x = "UMAP 1", y = "UMAP 2")
})
svglite("../figures/UMAP_TCRab_numbered_byCondition.svg", width = 22, height = 6)
print(patchwork::wrap_plots(cond_panels, nrow = 1))
dev.off()
cat("Saved: ../figures/UMAP_TCRab_numbered_byCondition.svg\n")

## ---- numbered-heatmap -------------------------------------------------------
# Genes as ROWS, clusters as COLUMNS; cluster column strip uses annotation colors
modules_hm <- list(
  "CD8 cytotoxic"  = c("Cd8a","Cd8b1","Gzmb","Gzmk","Prf1","Ccl5","Cxcr3","Nkg7","Fasl","Cst7","Klrk1","Eomes","Stat4"),
  "CD4 Th2/TFH"    = c("Cd4","Gata3","Il4","Il13","Cd40lg","Icos","Il21","Stat6","Gpr183","Bcl6","Cxcr5"),
  "Treg"           = c("Foxp3","Il2ra","Ctla4","Il10","Tnfrsf18","Tnfrsf4","Ikzf2","Tigit","Maf","Batf"),
  "TRM/activation" = c("Itgae","Itga1","Cd69","Cxcr6","Lag3","Havcr2","Pdcd1","Nr4a2","Nr4a3"),
  "Naive/TCM"      = c("Ccr7","Sell","S1pr1","Klf2","Tcf7","Lef1","Foxp1","Bcl2"),
  "Th17/MAIT17"    = c("Rorc","Rora","Il17a","Il17f","Il23r","Ccr6","Zbtb16","Il1r1","Tnfsf11","Trav1"),
  "Th1 axis"       = c("Ifng","Tbx21","Il12rb2"),
  "Cycling"        = c("Mki67","Cdk1","Ccnb1","Top2a","Birc5","Ube2c","Aurkb","Plk1","Nusap1","Cenpf")
)

hm_genes_ord <- unique(unlist(modules_hm))
hm_genes_ord <- hm_genes_ord[hm_genes_ord %in% rownames(t_cells_TCRab)]

Idents(t_cells_TCRab) <- "RNA_snn_res.0.8"
avg_cl <- AverageExpression(t_cells_TCRab, assays = "RNA", slot = "data",
                             features = hm_genes_ord, group.by = "RNA_snn_res.0.8")[[1]]
colnames(avg_cl) <- sub("^g", "", colnames(avg_cl))
col_ord <- as.character(sort(as.integer(colnames(avg_cl))))
avg_cl  <- avg_cl[hm_genes_ord, col_ord]

z_cl          <- t(scale(t(avg_cl)))
z_cl[is.na(z_cl)] <- 0; z_cl[z_cl > 2] <- 2; z_cl[z_cl < -1.5] <- -1.5

# Column annotation: predicted cell type for each cluster (no module strip)
col_ann_num <- data.frame(
  Annotation = factor(label_short_r08[col_ord], levels = unique(label_short_r08)),
  row.names  = col_ord, check.names = FALSE
)
ann_colors_num <- list(Annotation = pal_tcrab[levels(col_ann_num$Annotation)])
ann_colors_num$Annotation <- ann_colors_num$Annotation[!is.na(ann_colors_num$Annotation)]

ph_cl <- pheatmap::pheatmap(
  mat               = z_cl,          # genes as rows, clusters as columns
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  border_color      = NA,
  annotation_col    = col_ann_num,
  annotation_colors = ann_colors_num,
  fontsize_row      = 7,
  fontsize_col      = 10,
  angle_col         = 0,
  breaks            = seq(-1.5, 2, length.out = 101),
  color             = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100),
  main              = "αβ T cells — numbered clusters (avg exp, z-scored by gene)",
  silent            = TRUE
)

svglite("../figures/heatmap_TCRab_numbered.svg", width = 14, height = 12)
grid::grid.newpage()
grid::grid.draw(ph_cl$gtable)
dev.off()
cat("Saved: ../figures/heatmap_TCRab_numbered.svg\n")

cat("\n>>> Review UMAP_TCRab_numbered.svg + heatmap_TCRab_numbered.svg + marker CSVs.\n")
cat("    Update label_short_r08 in the USER CONFIGURATION section above, then re-run.\n\n")

## ============================================================================
## NOTE: Clusters 2 and 7 carry provisional labels in label_short_r08 above.
## Their identities are confirmed by the refinement steps below:
##   Cluster 7 "Th17/MAIT-like" → GMM separates MAIT cells; non-MAIT cells confirmed as "Th17-like"
##   Cluster 2 "Tfh/Tph-like" → sub-clustering + Tbx21/Batf test inconclusive; retained as "Tfh/Tph-like"
## ============================================================================

## ============================================================================
## STAGE 4b — ANNOTATION + CLUSTER REFINEMENT
## ============================================================================

## ---- 4.6-annotate -----------------------------------------------------------
t_cells_TCRab$celltype_r0.8_short <- unname(label_short_r08[as.character(Idents(t_cells_TCRab))])
t_cells_TCRab$celltype_merged     <- t_cells_TCRab$celltype_r0.8_short
cat("Initial cell types:", paste(sort(unique(t_cells_TCRab$celltype_merged)), collapse = ", "), "\n")

## ---- 4.7-mait-gmm -----------------------------------------------------------
# Cluster 7 = Th17 — identify and REMOVE innate-like contaminants (Zbtb16-high AND Trav1+)
# Innate-like cells are excluded from celltype_merged (set to NA) and not shown in final figures.
# Supplementary figures: Zbtb16 distribution, Innate-like vs Th17 by condition, TCR pairing.
cat("[MAIT GMM — cluster 7]\n")

cells_cl7  <- rownames(t_cells_TCRab@meta.data)[as.character(Idents(t_cells_TCRab)) == "7"]
zbtb16_nm  <- grep("^Zbtb16$", rownames(t_cells_TCRab), ignore.case = TRUE, value = TRUE)[1]
trav1_nm   <- grep("^Trav1$",  rownames(t_cells_TCRab), ignore.case = TRUE, value = TRUE)[1]
t_cells_TCRab$mait_label <- NA_character_

if (!is.na(zbtb16_nm) && length(cells_cl7) > 10) {
  zbtb16_vec <- FetchData(t_cells_TCRab, vars = zbtb16_nm, layer = "data")[cells_cl7, zbtb16_nm]

  thresh <- tryCatch({
    fit   <- mclust::Mclust(zbtb16_vec, G = 2, verbose = FALSE)
    mu    <- sort(fit$parameters$mean)
    sigma <- sqrt(fit$parameters$variance$sigmasq)
    if (length(sigma) == 1) sigma <- rep(sigma, 2)
    f_diff <- function(x)
      dnorm(x, mu[1], sigma[1]) * fit$parameters$pro[which(fit$parameters$mean == mu[1])] -
      dnorm(x, mu[2], sigma[2]) * fit$parameters$pro[which(fit$parameters$mean == mu[2])]
    uniroot(f_diff, interval = c(mu[1], mu[2]))$root
  }, error = function(e) {
    message("GMM failed (", conditionMessage(e), ") — histogram valley fallback")
    h     <- hist(zbtb16_vec, breaks = 30, plot = FALSE)
    lo    <- quantile(zbtb16_vec, 0.05); hi <- quantile(zbtb16_vec, 0.95)
    inner <- h$mids > lo & h$mids < hi
    if (any(inner)) h$mids[inner][which.min(h$density[inner])] else median(zbtb16_vec)
  })
  cat("  Zbtb16 threshold:", round(thresh, 4), "\n")

  zbtb16_pos <- zbtb16_vec > thresh
  trav1_pos  <- rep(FALSE, length(cells_cl7))
  if (!is.na(trav1_nm))
    trav1_pos <- FetchData(t_cells_TCRab, vars = trav1_nm, layer = "data")[cells_cl7, trav1_nm] > 0

  mait_cells <- cells_cl7[zbtb16_pos & trav1_pos]
  th17_cells <- setdiff(cells_cl7, mait_cells)

  t_cells_TCRab@meta.data[cells_cl7, "mait_label"]  <- "Th17 (conventional)"
  t_cells_TCRab@meta.data[mait_cells, "mait_label"] <- "Innate-like (MAIT/NKT-like)"
  # Innate-like cells excluded from final annotation (NA = invisible in all final figures)
  t_cells_TCRab@meta.data[mait_cells, "celltype_merged"] <- NA_character_
  # Th17 identity confirmed → relabel from provisional "Th17/MAIT" to final "Th17"
  t_cells_TCRab@meta.data[th17_cells, "celltype_merged"] <- "Th17-like"
  cat(sprintf("  Innate-like (removed): %d  |  Th17-like identity confirmed: %d cells\n",
              length(mait_cells), length(th17_cells)))

  # ---- Zbtb16 distribution histogram -----------------------------------------
  p_zbtb16 <- ggplot(data.frame(Zbtb16 = zbtb16_vec), aes(x = Zbtb16)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = "#6baed6", color = "white", alpha = 0.85) +
    geom_vline(xintercept = thresh, linetype = "dashed", color = "#d73027", linewidth = 0.9) +
    annotate("text", x = thresh + 0.03, y = Inf,
             label = sprintf("thresh = %.3f", thresh),
             hjust = 0, vjust = 1.8, size = 3.8, color = "#d73027") +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(title = "Cluster 7: Zbtb16 expression distribution",
         subtitle = "Dashed line = bimodal GMM threshold",
         x = "Zbtb16 (normalised expression)", y = "Density")
  svglite("../figures/cl7_Zbtb16_distribution.svg", width = 7, height = 5)
  print(p_zbtb16)
  dev.off()
  cat("Saved: ../figures/cl7_Zbtb16_distribution.svg\n")

  # ---- Innate-like vs Th17 by condition bar chart -----------------------------
  bar_df <- t_cells_TCRab@meta.data[cells_cl7, ] %>%
    dplyr::count(condition, mait_label, name = "n") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(prop = n / sum(n),
                  condition = factor(condition, levels = c("BAS","PD","REC","RCU"))) %>%
    dplyr::ungroup()
  p_cl7_bar <- ggplot(bar_df, aes(x = condition, y = prop,
                                   fill = factor(mait_label,
                                     levels = c("Innate-like (MAIT/NKT-like)","Th17 (conventional)")))) +
    geom_col(color = "white", width = 0.65) +
    scale_fill_manual(values = c("Innate-like (MAIT/NKT-like)" = "#d73027",
                                  "Th17 (conventional)"         = "#4E79A7"),
                      name = "Subtype") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.02), expand = c(0, 0)) +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank(), legend.title = element_blank()) +
    labs(title = "Cluster 7: Innate-like vs Th17 by Condition",
         x = "Condition", y = "Proportion of cluster 7 cells")
  svglite("../figures/cl7_innate_vs_Th17_byCondition.svg", width = 6, height = 5)
  print(p_cl7_bar)
  dev.off()
  cat("Saved: ../figures/cl7_innate_vs_Th17_byCondition.svg\n")

  # ---- Contamination fraction supplementary table -----------------------------
  cont_table <- t_cells_TCRab@meta.data[cells_cl7, ] %>%
    dplyr::count(condition, mait_label, name = "n") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(total_cl7 = sum(n), pct = round(100 * n / total_cl7, 1)) %>%
    dplyr::ungroup()
  write.csv(cont_table, "../objects/MAIT_contamination_per_condition.csv", row.names = FALSE)
  cat("Saved: ../objects/MAIT_contamination_per_condition.csv\n")
  print(cont_table)

  # ---- TCR pairing analysis for innate-like cells -----------------------------
  cat("\n[cl7_innate TCR pairing]\n")
  if (length(mait_cells) > 0) {
    mait_meta <- t_cells_TCRab@meta.data[mait_cells, ]
    tcr_col   <- intersect(c("CTgene","CTstrict","CTaa"), colnames(mait_meta))[1]
    cat("  Total innate-like cells :", length(mait_cells), "\n")

    if (!is.na(tcr_col)) {
      has_tcr   <- !is.na(mait_meta[[tcr_col]]) & mait_meta[[tcr_col]] != ""
      cat(sprintf("  Cells with %s data  : %d (%.1f%%)\n",
                  tcr_col, sum(has_tcr), 100 * mean(has_tcr)))
      cat(sprintf("  Cells without TCR data  : %d (%.1f%%)\n",
                  sum(!has_tcr), 100 * mean(!has_tcr)))

      # TCR recovery by condition
      cat("\nTCR recovery by condition:\n")
      recovery <- mait_meta %>%
        dplyr::mutate(has_tcr = !is.na(.data[[tcr_col]]) & .data[[tcr_col]] != "") %>%
        dplyr::count(condition, has_tcr) %>%
        tidyr::pivot_wider(names_from = has_tcr, values_from = n, values_fill = 0)
      print(as.data.frame(recovery))

      # Trav1 detection in CTgene
      cells_tcr    <- mait_cells[has_tcr]
      ctgene_vals  <- mait_meta[has_tcr, tcr_col]
      trav1_ctgene <- grepl("TRAV1[^0-9-]|TRAV1$|Trav1[^0-9-]|Trav1$", ctgene_vals, ignore.case = TRUE)
      traj33_ctgene <- grepl("TRAJ33|Traj33", ctgene_vals, ignore.case = TRUE)
      canonical    <- trav1_ctgene & traj33_ctgene
      cat(sprintf("\nTrav1 detected in %s: %d / %d cells with TCR data (%.1f%%)\n",
                  tcr_col, sum(trav1_ctgene), sum(has_tcr), 100 * mean(trav1_ctgene)))
      cat(sprintf("Canonical MAIT TCR (Trav1 + Traj33): %d cells\n", sum(canonical)))

      # Expression vs CTgene agreement
      if (!is.na(trav1_nm)) {
        trav1_expr <- FetchData(t_cells_TCRab, vars = trav1_nm, layer = "data")[mait_cells, trav1_nm] > 0
        trav1_expr_tcr <- trav1_expr[has_tcr]
        cat("\nExpression-based vs CTgene Trav1 agreement:\n")
        print(table(Trav1_expression = trav1_expr_tcr, Trav1_CTgene = trav1_ctgene))
      }

      # Sample entries
      cat(sprintf("\n[cl7_innate TCR pairing] Sample %s entries (first 10 with data):\n", tcr_col))
      sample_entries <- head(ctgene_vals[trav1_ctgene], 10)
      for (i in seq_along(sample_entries))
        cat(sprintf("  %d. %s\n", i, sample_entries[i]))

      # Save TCR summary CSV
      tcr_summary <- data.frame(
        cell          = cells_tcr,
        condition     = mait_meta[has_tcr, "condition"],
        CTgene        = ctgene_vals,
        Trav1_expr    = if (!is.na(trav1_nm)) FetchData(t_cells_TCRab, vars = trav1_nm,
                          layer = "data")[cells_tcr, trav1_nm] > 0 else NA,
        Trav1_CTgene  = trav1_ctgene,
        Traj33_CTgene = traj33_ctgene,
        canonical_MAIT = canonical
      )
      write.csv(tcr_summary, "../objects/cl7_innate_TCR_summary.csv", row.names = FALSE)
      cat("Saved: ../objects/cl7_innate_TCR_summary.csv\n")
    } else {
      cat("  No CTgene/CTstrict column found in metadata\n")
    }
    saveRDS(subset(t_cells_TCRab, cells = mait_cells), "../objects/cl7_innate_MAIT.rds")
    cat("Saved: ../objects/cl7_innate_MAIT.rds\n")
  }
} else {
  cat("  Zbtb16 not found or cluster 7 too small — skipping GMM\n")
}

## ---- 4.8-tph-subcluster -----------------------------------------------------
# Cluster 2 = Th1/Tph-like — test for genuine Th1/Tfh split vs co-expression
cat("[Tph sub-clustering — cluster 2]\n")
cells_cl2 <- rownames(t_cells_TCRab@meta.data)[as.character(Idents(t_cells_TCRab)) == "2"]
cat("  Cluster 2 cells:", length(cells_cl2), "\n")

if (length(cells_cl2) >= 30) {
  cl2 <- subset(t_cells_TCRab, cells = cells_cl2)
  cl2 <- JoinLayers(cl2)
  DefaultAssay(cl2) <- "RNA"
  cl2 <- NormalizeData(cl2, verbose = FALSE)
  cl2 <- FindVariableFeatures(cl2, nfeatures = 2000, verbose = FALSE)
  VariableFeatures(cl2) <- setdiff(VariableFeatures(cl2),
                                   grep(tcr_pat, VariableFeatures(cl2), value = TRUE, ignore.case = TRUE))
  cl2 <- ScaleData(cl2, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
  cl2 <- RunPCA(cl2, npcs = 30, verbose = FALSE, seed.use = 123)
  cl2 <- RunUMAP(cl2, reduction = "pca", dims = 1:10,
                  n.neighbors = 15, min.dist = 0.2, reduction.name = "umap.cl2", seed.use = 123)
  cl2 <- FindNeighbors(cl2, reduction = "pca", dims = 1:10, k.param = 15)

  for (res in c(0.2, 0.4, 0.6, 0.8)) {
    col_nm <- paste0("cl2_res", gsub("\\.", "", sprintf("%.1f", res)))
    cl2 <- FindClusters(cl2, resolution = res, algorithm = 1, cluster.name = col_nm)
    cat("  Resolution", res, ":", length(unique(cl2@meta.data[[col_nm]])), "sub-clusters\n")
  }

  pca_cl2  <- cl2@reductions$pca@cell.embeddings[, 1:10]
  dist_cl2 <- dist(pca_cl2)
  compute_sil <- function(labels, d) {
    if (length(unique(labels)) < 2) return(0)
    mean(cluster::silhouette(as.integer(factor(labels)), d)[, 3])
  }
  sil_scores <- c("0.2" = compute_sil(cl2$cl2_res02, dist_cl2),
                  "0.4" = compute_sil(cl2$cl2_res04, dist_cl2),
                  "0.6" = compute_sil(cl2$cl2_res06, dist_cl2),
                  "0.8" = compute_sil(cl2$cl2_res08, dist_cl2))
  cat("  Silhouette:", paste(names(sil_scores), round(sil_scores, 3), sep = "=", collapse = ", "), "\n")
  best_res <- names(which.max(sil_scores))
  best_col <- paste0("cl2_res", gsub("\\.", "", sprintf("%.1f", as.numeric(best_res))))
  Idents(cl2) <- best_col
  cat("  Best resolution:", best_res, "\n")

  tbx21_nm <- grep("^Tbx21$", rownames(cl2), ignore.case = TRUE, value = TRUE)[1]
  batf_nm  <- grep("^Batf$",  rownames(cl2), ignore.case = TRUE, value = TRUE)[1]
  tph_split <- FALSE
  peak_tbx <- NA_character_; peak_bat <- NA_character_
  if (!is.na(tbx21_nm) && !is.na(batf_nm)) {
    tb_expr <- FetchData(cl2, vars = c(tbx21_nm, batf_nm), layer = "data")
    cl2_ids <- as.character(Idents(cl2))
    sub_ids <- sort(unique(cl2_ids))
    if (length(sub_ids) >= 2) {
      mean_tbx <- tapply(tb_expr[[tbx21_nm]], cl2_ids, mean)
      mean_bat <- tapply(tb_expr[[batf_nm]],  cl2_ids, mean)
      peak_tbx <- names(which.max(mean_tbx))
      peak_bat <- names(which.max(mean_bat))
      cond_a   <- peak_tbx != peak_bat
      cond_b   <- (mean_tbx[peak_tbx] / max(mean_tbx[setdiff(sub_ids, peak_tbx)]) >= 2.0) &&
                  (mean_bat[peak_bat] / max(mean_bat[setdiff(sub_ids, peak_bat)]) >= 2.0)
      cond_c   <- (mean_tbx[peak_bat] <= median(tb_expr[[tbx21_nm]])) &&
                  (mean_bat[peak_tbx] <= median(tb_expr[[batf_nm]]))
      tph_split <- cond_a && cond_b && cond_c
      cat("  Tbx21/Batf split — a:", cond_a, "b:", cond_b, "c:", cond_c, "=> split:", tph_split, "\n")
    }
  }

  if (tph_split) {
    cl2$cl2_subtype <- ifelse(as.character(Idents(cl2)) == peak_tbx, "Th1-like", "Tfh-like")
    subtype_vec <- setNames(cl2$cl2_subtype, colnames(cl2))
    t_cells_TCRab@meta.data[cells_cl2, "celltype_merged"] <- subtype_vec[cells_cl2]
    cat("  Cluster 2 → Th1-like / Tfh-like split\n")
  } else {
    # Tbx21/Batf co-express across all sub-clusters → no algorithmic basis for a Th1/Tfh split.
    # Identity inconclusive; retained as "Tfh/Tph-like" (see annotation conclusion below).
    cl2$cl2_subtype <- "Tfh/Tph-like"
    t_cells_TCRab@meta.data[cells_cl2, "celltype_merged"] <- "Tfh/Tph-like"
    cat("  Cluster 2 retained as Tfh/Tph-like:", length(cells_cl2), "cells\n")
  }

  # ---- 4.8 annotation conclusion -----------------------------------------------
  # Final label: "Tfh/Tph-like"
  #
  # Three independent lines of evidence were evaluated:
  #
  # 1. Sub-clustering structure (algorithmic):
  #    - Silhouette scores at res 0.2–0.8 all near zero (best = 0.094 at res 0.4 without res 0.2;
  #      0.324 at res 0.2 with only 2 sub-clusters). No resolution produces clearly separated lobes.
  #    - No real sub-structure exists; a genuinely mixed Th1+Tph population would separate
  #      cleanly at some resolution with high silhouette scores (>0.3).
  #
  # 2. Tbx21 / Batf marker test:
  #    - Conditions a (peaks differ), b (2× enrichment), c (mutual exclusion) all FALSE.
  #    - Tbx21 and Batf co-express uniformly across sub-clusters.
  #    - No algorithmic basis to split the cluster into Th1 (Tbx21-hi) and Tph (Batf-hi).
  #
  # 3. Literature and heatmap review (inconclusive):
  #    - CXCR5 low → supports Tph over classical TFH (TFH are CXCR5-hi)
  #    - BCL6 intermediate → usually high in TFH, but in inflamed tissue Tph-like cells
  #      can also upregulate BCL6 transiently; inconclusive here
  #    - Top-expressed genes include: Cd40lg (B cell help active), Batf (Tfh/Tph program),
  #      Il21 (B cell help), Rora (activated tissue CD4), Cxcr6 (peripheral/tissue homing,
  #      supports Tph over GC-TFH which are Cxcr5-hi/Cxcr6-lo)
  #
  # Summary: This cluster exhibits a mixed B-helper transcriptional program (PDCD1+, IL21+,
  # ICOS+, CD40LG+, CXCR6+, BCL6 intermediate) with low CXCR5 and CXCL13, consistent with
  # a Tfh/Tph-like phenotype that cannot be confidently resolved into either canonical identity.

  # ---- Save per-sub-cluster Tbx21/Batf means ----------------------------------
  if (!is.na(tbx21_nm) && !is.na(batf_nm) && exists("mean_tbx") && exists("mean_bat")) {
    subcl_means <- data.frame(
      sub_cluster  = names(mean_tbx),
      mean_Tbx21   = round(as.numeric(mean_tbx), 4),
      mean_Batf    = round(as.numeric(mean_bat[names(mean_tbx)]), 4),
      silhouette   = round(sil_scores[best_res], 4),
      best_res     = best_res,
      split_called = tph_split,
      cond_a_peaks_differ     = if (exists("cond_a")) cond_a else NA,
      cond_b_2x_enriched      = if (exists("cond_b")) cond_b else NA,
      cond_c_mutually_excl    = if (exists("cond_c")) cond_c else NA,
      row.names = NULL
    )
    write.csv(subcl_means, "../objects/cl2_subcluster_Tbx21_Batf_means.csv", row.names = FALSE)
    cat("Saved: ../objects/cl2_subcluster_Tbx21_Batf_means.csv\n")
    print(subcl_means[, c("sub_cluster","mean_Tbx21","mean_Batf")])
  }

  # ---- Multi-resolution UMAP (0.2 / 0.4 / 0.6 / 0.8) -------------------------
  cl2_emb_df <- as.data.frame(cl2@reductions[["umap.cl2"]]@cell.embeddings)
  colnames(cl2_emb_df) <- c("UMAP_1", "UMAP_2")
  res_cols <- c("cl2_res02", "cl2_res04", "cl2_res06", "cl2_res08")
  res_labs <- c("Resolution = 0.2", "Resolution = 0.4", "Resolution = 0.6", "Resolution = 0.8")

  p_res_list <- mapply(function(col, lab) {
    df <- cbind(cl2_emb_df, cluster = as.factor(cl2@meta.data[[col]]))
    n_cl <- nlevels(df$cluster)
    ggplot(df, aes(UMAP_1, UMAP_2, color = cluster)) +
      geom_point(size = 0.6, alpha = 0.8) +
      ggrepel::geom_label_repel(
        data = df %>% dplyr::group_by(cluster) %>%
               dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop"),
        aes(label = cluster), size = 3, label.size = 0.2, fill = "white", alpha = 0.8,
        show.legend = FALSE
      ) +
      scale_color_manual(values = scales::hue_pal()(n_cl), guide = "none") +
      theme_classic(base_size = 11) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 9), axis.text = element_text(size = 8)) +
      labs(title = lab, x = "umapcl2_1", y = "umapcl2_2")
  }, res_cols, res_labs, SIMPLIFY = FALSE)

  p_res_grid <- patchwork::wrap_plots(p_res_list, nrow = 1) +
    patchwork::plot_annotation(
      title = "Cluster 2 sub-clustering at multiple resolutions",
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    )
  svglite("../figures/cl2_subcluster_resolution_sweep.svg", width = 20, height = 5)
  print(p_res_grid)
  dev.off()
  cat("Saved: ../figures/cl2_subcluster_resolution_sweep.svg\n")

  # ---- Tfh/Tph-like marker feature plots on cl2 UMAP --------------------------
  # 9 markers (incl. Cxcr3 — helps resolve Tph vs Th1 axis) in a 3×3 grid on cl2 sub-UMAP
  tph_panel_raw <- c("Tbx21","Batf","Cxcr5","Bcl6","Ifng","Pdcd1","Cxcr6","Maf","Cxcr3")
  tph_panel_nms <- sapply(tph_panel_raw, resolve_gene, rownames(cl2))
  tph_panel_nms <- tph_panel_nms[!is.na(tph_panel_nms)]

  cl2_umap_emb <- cl2@reductions[["umap.cl2"]]@cell.embeddings
  cl2_umap_df  <- data.frame(
    UMAP_1 = cl2_umap_emb[, 1],
    UMAP_2 = cl2_umap_emb[, 2],
    row.names = rownames(cl2_umap_emb)
  )
  feat_expr <- FetchData(cl2, vars = unname(tph_panel_nms), layer = "data")

  # ---- Tph marker feature plots on cl2 UMAP ----------------------------------
  # Expression-colored scatter: cells sorted by expression (zero-expr rendered first
  # as background, high-expr on top for visual foreground). Title shows % expressing.
  dp_list <- lapply(unname(tph_panel_nms), function(gn) {
    plot_df  <- cbind(cl2_umap_df, expr = feat_expr[[gn]])
    frac_pos <- sum(plot_df$expr > 0) / nrow(plot_df)
    # Sort so zero-expression cells render first (background), high-expression on top
    plot_df  <- plot_df[order(plot_df$expr), ]

    ggplot(plot_df, aes(UMAP_1, UMAP_2, color = expr)) +
      geom_point(size = 1.2, alpha = 0.9) +
      scale_color_gradientn(
        colors = c("grey90", "#FEE5D9", "#FC9272", "#DE2D26", "#67000D"),
        name   = "Expr",
        guide  = guide_colorbar(barwidth = 0.4, barheight = 3,
                                title.theme = element_text(size = 7),
                                label.theme = element_text(size = 6))
      ) +
      theme_classic(base_size = 10) +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 7), axis.text = element_text(size = 6),
            plot.title = element_text(size = 11, face = "bold")) +
      labs(title = sprintf("%s  (%.0f%% expr)", gn, frac_pos * 100),
           x = "umapcl2_1", y = "umapcl2_2")
  })
  svglite("../figures/cl2_Tph_featureplot.svg", width = 18, height = 18)
  print(patchwork::wrap_plots(dp_list, nrow = 3, ncol = 3))
  dev.off()
  cat("Saved: ../figures/cl2_Tph_featureplot.svg\n")

  # ---- Tph marker expression summary table ------------------------------------
  all_expr_full <- FetchData(t_cells_TCRab, vars = unname(tph_panel_nms), layer = "data")
  cl2_mask      <- rownames(t_cells_TCRab@meta.data) %in% cells_cl2
  mean_cl2   <- colMeans(all_expr_full[ cl2_mask, unname(tph_panel_nms), drop = FALSE])
  mean_rest  <- colMeans(all_expr_full[!cl2_mask, unname(tph_panel_nms), drop = FALSE])
  pct_cl2    <- colMeans(all_expr_full[ cl2_mask, unname(tph_panel_nms), drop = FALSE] > 0) * 100
  pct_rest   <- colMeans(all_expr_full[!cl2_mask, unname(tph_panel_nms), drop = FALSE] > 0) * 100
  tph_summary <- data.frame(
    gene        = unname(tph_panel_nms),
    mean_cl2    = round(mean_cl2,  4),
    mean_rest   = round(mean_rest, 4),
    log2FC      = round(log2((mean_cl2 + 0.01) / (mean_rest + 0.01)), 3),
    pct_cl2     = round(pct_cl2,  1),
    pct_rest    = round(pct_rest, 1),
    row.names   = NULL
  )
  write.csv(tph_summary, "../objects/cl2_Tph_marker_summary.csv", row.names = FALSE)
  cat("Saved: ../objects/cl2_Tph_marker_summary.csv\n")
  print(tph_summary)

  # ---- Tph marker violin + cell dots ------------------------------------------
  # Confirmed Tph vs all other αβ T cells; dots rasterized to keep file size small
  tph_vln_df <- all_expr_full %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(group = factor(
      ifelse(cell %in% cells_cl2, "Tfh/Tph-like", "Other T cells"),
      levels = c("Tfh/Tph-like", "Other T cells")
    )) %>%
    tidyr::pivot_longer(cols = all_of(unname(tph_panel_nms)), names_to = "gene", values_to = "expr") %>%
    dplyr::mutate(gene = factor(gene, levels = unname(tph_panel_nms)))

  p_tph_vln <- ggplot(tph_vln_df, aes(x = group, y = expr, fill = group, color = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.55, color = NA) +
    ggrastr::geom_jitter_rast(width = 0.2, size = 0.25, alpha = 0.25, shape = 16,
                              raster.dpi = 150) +
    facet_wrap(~ gene, scales = "free_y", nrow = 3) +
    scale_fill_manual(values  = c("Tfh/Tph-like" = "#E15759", "Other T cells" = "grey60")) +
    scale_color_manual(values = c("Tfh/Tph-like" = "#E15759", "Other T cells" = "grey60")) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          strip.background = element_blank(), strip.text = element_text(face = "bold")) +
    labs(title = "Tfh/Tph-like marker panel — cluster 2 vs other αβ T cells", y = "Normalized expression")
  svglite("../figures/cl2_Tph_violin.svg", width = 12, height = 10)
  print(p_tph_vln)
  dev.off()
  cat("Saved: ../figures/cl2_Tph_violin.svg\n")

  # ---- Extended Tfh/Tph program violin panel (16 genes) -----------------------
  # Second violin panel: broader Tfh/Tph B-helper + naive/costimulation program.
  tph_ext_req <- c("Cxcl13", "Icos", "Il21", "Cd40lg", "Pdcd1", "Tox", "Batf", "Prdm1",
                   "Bcl6", "Cxcr5", "Cxcr6", "Sell", "Ccr7", "Tnfrsf4", "Tnfsf8", "Cd200")
  tph_ext_nms <- sapply(tph_ext_req, resolve_gene, rownames(t_cells_TCRab))
  missing_ext <- tph_ext_req[is.na(tph_ext_nms)]
  if (length(missing_ext))
    cat("  Extended Tph panel — genes not found (dropped):", paste(missing_ext, collapse = ", "), "\n")
  tph_ext_nms <- tph_ext_nms[!is.na(tph_ext_nms)]

  ext_expr <- FetchData(t_cells_TCRab, vars = unname(tph_ext_nms), layer = "data")
  tph_ext_df <- ext_expr %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(group = factor(
      ifelse(cell %in% cells_cl2, "Tfh/Tph-like", "Other T cells"),
      levels = c("Tfh/Tph-like", "Other T cells")
    )) %>%
    tidyr::pivot_longer(cols = all_of(unname(tph_ext_nms)), names_to = "gene", values_to = "expr") %>%
    dplyr::mutate(gene = factor(gene, levels = unname(tph_ext_nms)))

  p_tph_ext <- ggplot(tph_ext_df, aes(x = group, y = expr, fill = group, color = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.55, color = NA) +
    ggrastr::geom_jitter_rast(width = 0.2, size = 0.25, alpha = 0.25, shape = 16,
                              raster.dpi = 150) +
    facet_wrap(~ gene, scales = "free_y", nrow = 4) +
    scale_fill_manual(values  = c("Tfh/Tph-like" = "#E15759", "Other T cells" = "grey60")) +
    scale_color_manual(values = c("Tfh/Tph-like" = "#E15759", "Other T cells" = "grey60")) +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          strip.background = element_blank(), strip.text = element_text(face = "bold")) +
    labs(title = "Extended Tfh/Tph program — cluster 2 vs other αβ T cells", y = "Normalized expression")
  svglite("../figures/cl2_Tph_violin_extended.svg", width = 13, height = 12)
  print(p_tph_ext)
  dev.off()
  cat("Saved: ../figures/cl2_Tph_violin_extended.svg\n")

  # ---- Combined: feature panel (top) + violin panel (bottom) ------------------
  p_feature_panel <- patchwork::wrap_plots(dp_list, nrow = 3, ncol = 3)
  combined_tph <- p_feature_panel / p_tph_vln +
    patchwork::plot_layout(heights = c(1.6, 1)) +
    patchwork::plot_annotation(
      title = "Tfh/Tph-like cluster 2 characterization",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  svglite("../figures/cl2_Tph_combined.svg", width = 18, height = 28)
  print(combined_tph)
  dev.off()
  cat("Saved: ../figures/cl2_Tph_combined.svg\n")

  saveRDS(cl2, "../objects/cl2_Tph_subcluster.rds")
  cat("Saved: ../objects/cl2_Tph_subcluster.rds\n")

  # Alluvial plot (resolution stability)
  alluvial_df <- data.frame(res_0.4 = cl2$cl2_res04, res_0.6 = cl2$cl2_res06,
                             res_0.8 = cl2$cl2_res08) %>%
    dplyr::count(res_0.4, res_0.6, res_0.8, name = "Freq")
  p_alluvial <- ggplot(alluvial_df,
    aes(y = Freq, axis1 = res_0.4, axis2 = res_0.6, axis3 = res_0.8)) +
    geom_alluvium(aes(fill = res_0.4), width = 1/3) +
    geom_stratum(width = 1/3, fill = "grey88", color = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(limits = c("Res 0.4","Res 0.6","Res 0.8"), expand = c(0.15, 0.05)) +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(title = "Cluster 2 sub-clustering stability", y = "Cells", fill = "Sub-cluster (res 0.4)")
  svglite("../figures/UMAP_cl2_alluvial.svg", width = 8, height = 5)
  print(p_alluvial)
  dev.off()
  cat("Saved: ../figures/UMAP_cl2_alluvial.svg\n")

  if (!is.na(tbx21_nm) && !is.na(batf_nm)) {
    scatter_df <- FetchData(cl2, vars = c(tbx21_nm, batf_nm, best_col), layer = "data")
    p_scatter <- ggplot(scatter_df,
      aes(x = .data[[tbx21_nm]], y = .data[[batf_nm]], color = .data[[best_col]])) +
      geom_point(alpha = 0.7, size = 1.5) +
      theme_classic(base_size = 12) +
      theme(panel.grid = element_blank()) +
      labs(title = "Cluster 2: Tbx21 vs Batf", color = "Sub-cluster", x = tbx21_nm, y = batf_nm)
    svglite("../figures/cl2_Tbx21_Batf_scatter.svg", width = 5, height = 5)
    print(p_scatter)
    dev.off()
    cat("Saved: ../figures/cl2_Tbx21_Batf_scatter.svg\n")
  }
} else {
  cat("  Cluster 2 too small —", length(cells_cl2), "cells, skipping sub-clustering\n")
}

## ---- 4.8b-tph-identity-heatmap -----------------------------------------------
# Compact heatmap: Tfh/Tph-like identity genes × all annotated T cell subtypes.
# Shows PD-1hi/CXCR5-low/PRDM1+/IL21+/CXCL13-low pattern supporting Tfh/Tph-like annotation.
tph_id_genes <- c("Pdcd1", "Cxcr5", "Bcl6", "Prdm1", "Il21", "Cxcl13",
                   "Maf", "Slamf6", "Ccr2", "Ccr5", "Cx3cr1")
tph_id_genes <- intersect(tph_id_genes, rownames(t_cells_TCRab))

if (length(tph_id_genes) >= 3) {
  # confirmed cell types only — drop provisional names and NA
  ct_tph_hm <- names(pal_tcrab)[
    names(pal_tcrab) %in% unique(na.omit(t_cells_TCRab$celltype_merged)) &
    !names(pal_tcrab) %in% c("Tph/Th1-like", "Th17/MAIT-like")
  ]

  avg_tph <- sapply(ct_tph_hm, function(ct) {
    idx <- !is.na(t_cells_TCRab$celltype_merged) & t_cells_TCRab$celltype_merged == ct
    if (!any(idx)) return(setNames(rep(NA_real_, length(tph_id_genes)), tph_id_genes))
    colMeans(FetchData(t_cells_TCRab, vars = tph_id_genes,
                       cells = rownames(t_cells_TCRab@meta.data)[idx], layer = "data"))
  })  # rows = genes, cols = cell types

  z_tph <- t(scale(t(avg_tph)))           # z-score each gene across cell types
  z_tph[is.nan(z_tph)] <- 0
  z_tph <- pmax(pmin(z_tph, 2), -1.5)    # same clamp as main heatmap

  col_ann_tph <- data.frame(Celltype = factor(ct_tph_hm, levels = ct_tph_hm),
                              row.names = ct_tph_hm)
  ann_col_tph <- list(Celltype = pal_tcrab[ct_tph_hm])

  ph_tph_id <- pheatmap::pheatmap(
    mat                      = z_tph,
    cluster_rows             = FALSE,
    cluster_cols             = TRUE,
    clustering_distance_cols = "correlation",
    clustering_method        = "ward.D2",
    show_rownames            = TRUE,
    show_colnames            = TRUE,
    border_color             = NA,
    annotation_col           = col_ann_tph,
    annotation_colors        = ann_col_tph,
    annotation_legend        = TRUE,
    fontsize                 = 11,
    fontsize_row             = 10,
    angle_col                = 45,
    breaks                   = seq(-1.5, 2, length.out = 101),
    color                    = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100),
    main                     = "Tfh/Tph-like identity markers (avg exp, z-scored)",
    silent                   = TRUE
  )
  svglite("../figures/heatmap_Tph_identity.svg", width = 9, height = 5)
  grid::grid.draw(ph_tph_id$gtable)
  dev.off()
  cat("Saved: ../figures/heatmap_Tph_identity.svg\n")
} else {
  cat("  Tph identity heatmap skipped — genes not found in dataset\n")
}

## ---- 4.8c-cluster6-mixed-lineage --------------------------------------------
# Cluster 6 (provisional "Trm-like CD8") is a tissue-resident/IEL-like population
# that is only ~50% CD8-marker-positive (weakest of all CD8 clusters) and ~26%
# Cd4+/CD8-. scDblFinder ruled out doublets (1.2% — lowest of any T-cell cluster; Tcell_doublet_check.R,
# see objects/cl6_doublet_check.csv). The evidence below shows it is a genuine
# CD4/CD8 MIXTURE with two functionally distinct arms (cytotoxic CD8 + helper CD4),
# not a clean CD8 cluster and not a single "CD8-helper" population. Because ~24% of
# cells carry no lineage marker (dropout) and cannot be split reliably per-cell, it
# is annotated honestly as "Trm type mixed CD8+ CD4" (mixed = neither CD4 nor CD8,
# so it drops out of all lineage-specific downstream analyses).
cat("[Cluster 6 mixed-lineage assessment]\n")

cells_cl6 <- rownames(t_cells_TCRab@meta.data)[as.character(Idents(t_cells_TCRab)) == "6"]

if (length(cells_cl6) > 10) {
  # ---- raw (>0) mutually-exclusive lineage call --------------------------------
  raw_c6 <- FetchData(t_cells_TCRab, vars = c("Cd4", "Cd8a", "Cd8b1"),
                      cells = cells_cl6, layer = "counts")
  cd8_raw <- raw_c6$Cd8a > 0 | raw_c6$Cd8b1 > 0
  cd4_raw <- raw_c6$Cd4 > 0
  lin_raw <- ifelse( cd8_raw & !cd4_raw, "CD8",
             ifelse( cd4_raw & !cd8_raw, "CD4",
             ifelse( cd8_raw &  cd4_raw, "double-positive", "double-negative")))
  names(lin_raw) <- cells_cl6
  n_cd8_raw <- sum(lin_raw == "CD8"); n_cd4_raw <- sum(lin_raw == "CD4")
  n_dp <- sum(lin_raw == "double-positive"); n_dn <- sum(lin_raw == "double-negative")
  cat(sprintf("  raw >0: CD8=%d CD4=%d double+=%d double-(dropout)=%d\n",
              n_cd8_raw, n_cd4_raw, n_dp, n_dn))

  # ---- module scores (for the two-arm scatter; CD4 uses Cd40lg helper co-marker)
  cd8_set <- list(intersect(c("Cd8a", "Cd8b1"), rownames(t_cells_TCRab)))
  cd4_set <- list(intersect(c("Cd4", "Cd40lg"), rownames(t_cells_TCRab)))
  t_cells_TCRab <- AddModuleScore(t_cells_TCRab, features = cd8_set, name = "CD8score", seed = 123)
  t_cells_TCRab <- AddModuleScore(t_cells_TCRab, features = cd4_set, name = "CD4score", seed = 123)

  # ---- kNN label transfer (majority vote) on a local cluster-6 re-embedding ----
  # Local PCA (TCR genes excluded, same recipe as reclustering), k=30 neighbours.
  # Keep the raw call for marker-resolvable cells; assign the double-negative
  # (dropout) cells a lineage by majority vote of their neighbours' raw call; exclude
  # the few double-positive cells. Resolves the ~24% dropout without guessing per gene.
  # This is categorical label transfer, not expression imputation — no counts altered.
  cl6 <- subset(t_cells_TCRab, cells = cells_cl6)
  cl6 <- JoinLayers(cl6)
  cl6 <- NormalizeData(cl6, verbose = FALSE)
  cl6 <- FindVariableFeatures(cl6, nfeatures = 2000, verbose = FALSE)
  VariableFeatures(cl6) <- setdiff(VariableFeatures(cl6),
                                   grep(tcr_pat, VariableFeatures(cl6), value = TRUE, ignore.case = TRUE))
  cl6 <- ScaleData(cl6, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
  cl6 <- RunPCA(cl6, npcs = 30, seed.use = 123, verbose = FALSE)
  pcs_c6 <- Embeddings(cl6, "pca")[cells_cl6, 1:30]

  set.seed(123)
  knn_c6 <- FNN::get.knn(pcs_c6, k = 30)$nn.index
  rownames(knn_c6) <- cells_cl6
  lin_knn <- lin_raw
  for (i in which(lin_raw == "double-negative")) {
    nbr <- lin_raw[knn_c6[i, ]]
    nbr <- nbr[nbr %in% c("CD8", "CD4")]
    lin_knn[i] <- if (length(nbr)) names(sort(table(nbr), decreasing = TRUE))[1] else NA_character_
  }
  lin_knn[lin_knn == "double-positive"]  <- NA_character_    # DP excluded from composition
  n_cd8_knn <- sum(lin_knn == "CD8", na.rm = TRUE)
  n_cd4_knn <- sum(lin_knn == "CD4", na.rm = TRUE)
  n_res_knn <- n_cd8_knn + n_cd4_knn
  pct_cd8_knn <- round(100 * n_cd8_knn / n_res_knn)
  cat(sprintf("  kNN label transfer (majority vote): CD8=%d CD4=%d (%.0f%% / %.0f%% of %d resolvable)\n",
              n_cd8_knn, n_cd4_knn, pct_cd8_knn, 100 - pct_cd8_knn, n_res_knn))

  # second concordant method: raw-resolvable only (no label transfer)
  pct_cd8_raw <- round(100 * n_cd8_raw / (n_cd8_raw + n_cd4_raw))

  # ---- evidence table ----------------------------------------------------------
  comp_tbl <- data.frame(
    method = c("Raw marker+ (resolvable only)", "kNN label transfer (majority vote)", "double-positive (excluded)"),
    n_CD8  = c(n_cd8_raw, n_cd8_knn, NA),
    n_CD4  = c(n_cd4_raw, n_cd4_knn, NA),
    pct_CD8 = c(pct_cd8_raw, pct_cd8_knn, NA),
    pct_CD4 = c(100 - pct_cd8_raw, 100 - pct_cd8_knn, NA),
    note   = c(sprintf("%d cells with exactly one lineage marker", n_cd8_raw + n_cd4_raw),
               sprintf("%d dropout cells assigned by k=30-neighbour majority vote", n_dn),
               sprintf("%d Cd4+Cd8+ cells excluded", n_dp))
  )
  write.csv(comp_tbl, "../objects/cl6_mixed_lineage_composition.csv", row.names = FALSE)
  cat("Saved: ../objects/cl6_mixed_lineage_composition.csv\n")

  # ---- separability test: does unsupervised sub-clustering split CD4 from CD8? --
  # If lineage were a separable axis, sub-clusters would be lineage-pure. Instead
  # sub-clustering the local re-embedding splits on activation state and every
  # sub-cluster stays CD4/CD8-mixed → the merge is justified, not an artefact.
  cl6 <- FindNeighbors(cl6, reduction = "pca", dims = 1:10, k.param = 15, verbose = FALSE)
  cl6 <- RunUMAP(cl6, reduction = "pca", dims = 1:10, reduction.name = "umap.cl6",
                 seed.use = 123, verbose = FALSE)
  dist_cl6 <- dist(Embeddings(cl6, "pca")[, 1:10])
  sil_c6 <- setNames(numeric(4), c("0.2", "0.4", "0.6", "0.8"))
  for (res in c(0.2, 0.4, 0.6, 0.8)) {
    col_nm <- paste0("cl6_res", gsub("\\.", "", sprintf("%.1f", res)))
    cl6    <- FindClusters(cl6, resolution = res, algorithm = 1, cluster.name = col_nm, verbose = FALSE)
    lab    <- cl6@meta.data[[col_nm]]
    sil_c6[as.character(res)] <- if (length(unique(lab)) > 1)
      mean(cluster::silhouette(as.integer(factor(lab)), dist_cl6)[, 3]) else 0
  }
  best_res_c6 <- names(which.max(sil_c6))
  best_col_c6 <- paste0("cl6_res", gsub("\\.", "", sprintf("%.1f", as.numeric(best_res_c6))))
  cl6$subcluster <- cl6@meta.data[[best_col_c6]]
  cat(sprintf("  cluster-6 sub-clustering silhouette: %s | best res %s (%d sub-clusters)\n",
              paste(names(sil_c6), round(sil_c6, 3), sep = "=", collapse = ", "),
              best_res_c6, length(unique(cl6$subcluster))))

  # lineage composition within each sub-cluster, from RAW marker+ calls (resolvable
  # cells only; dropout excluded, no imputation) — the "not separable" evidence.
  subD <- data.frame(subcluster = as.character(cl6$subcluster),
                     lineage    = lin_raw[colnames(cl6)]) %>%
    dplyr::filter(lineage %in% c("CD8", "CD4")) %>%
    dplyr::count(subcluster, lineage, name = "n") %>%
    dplyr::group_by(subcluster) %>%
    dplyr::mutate(pct = 100 * n / sum(n), n_sub = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sub_lab  = sprintf("sub %s\n(n=%d)", subcluster, n_sub),
                  lineage  = factor(lineage, levels = c("CD4", "CD8")))
  write.csv(subD[, c("subcluster", "lineage", "n", "pct")],
            "../objects/cl6_subcluster_lineage_composition.csv", row.names = FALSE)
  cat("Saved: ../objects/cl6_subcluster_lineage_composition.csv\n")

  # =============================================================================
  # SUPPLEMENTARY FIGURE — justification of the mixed-lineage annotation
  # =============================================================================
  base_c6 <- theme_classic(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.key.size = unit(0.4, "cm"))

  # ---- Panel A: lineage proportion among RESOLVABLE cells only -----------------
  # CD8 vs CD4 of the marker-resolvable cells; dropout (n_dn) and double-positive
  # (n_dp) are not drawn — they are reported in the subtitle and caption instead.
  n_resolv <- n_cd8_raw + n_cd4_raw
  compA <- data.frame(lineage = c("CD8", "CD4"), n = c(n_cd8_raw, n_cd4_raw))
  compA$pct <- 100 * compA$n / sum(compA$n)
  compA$lineage <- factor(compA$lineage, levels = c("CD8", "CD4"))   # CD8 top, CD4 bottom
  pA <- ggplot(compA, aes(x = "Cluster 6", y = pct, fill = lineage)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%s\n%.0f%%\n(n=%d)", lineage, pct, n)),
              position = position_stack(vjust = 0.5), size = 2.9, color = "white",
              fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = c("CD8" = "#4E79A7", "CD4" = "#E15759"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    base_c6 + theme(legend.position = "none",
                    plot.subtitle = element_text(hjust = 0.5, size = 7)) +
    labs(title = "A  Lineage proportion",
         subtitle = sprintf("of %d resolvable cells\n(%d dropout + %d DP unassigned)", n_resolv, n_dn, n_dp),
         x = NULL, y = "% of resolvable cells")

  # ---- Panel B: CD8 vs CD4 module-score scatter, colored by raw marker call ----
  msB <- t_cells_TCRab@meta.data[cells_cl6, c("CD8score1", "CD4score1")]
  msB$call <- factor(lin_raw[rownames(msB)],
    levels = c("CD8", "CD4", "double-negative", "double-positive"),
    labels = c("CD8", "CD4", "unresolved (dropout)", "double-positive"))
  pB <- ggplot(msB, aes(CD8score1, CD4score1, color = call)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_color_manual(values = c("CD8" = "#4E79A7", "CD4" = "#E15759",
                                  "unresolved (dropout)" = "grey70",
                                  "double-positive" = "grey40"), name = "Marker call (raw)") +
    base_c6 + labs(title = "B  Two arms: CD8 vs CD4 module score (per cell)",
                   x = "CD8 module score (Cd8a, Cd8b1)", y = "CD4 module score (Cd4, Cd40lg)")

  # ---- Panel C: program markers, cluster 6 vs clean CD8 & clean CD4 ------------
  ref_cd8 <- "Teff CD8-like"; ref_cd4 <- "Th2-like"
  progC_genes <- intersect(c("Cd4", "Cd8a", "Cd8b1",                 # lineage
                             "Itga1", "Itgae", "Cxcr6", "Id2",       # Trm/IEL
                             "Fos"),                                 # activation
                          rownames(t_cells_TCRab))
  cellsC <- rownames(t_cells_TCRab@meta.data)[
    t_cells_TCRab$celltype_merged %in% c(ref_cd8, ref_cd4) & !is.na(t_cells_TCRab$celltype_merged)]
  grpC <- setNames(t_cells_TCRab$celltype_merged[cellsC], cellsC)
  grpC[cells_cl6] <- "Cluster 6 (mixed)"
  exprC <- FetchData(t_cells_TCRab, vars = progC_genes,
                    cells = c(cells_cl6, cellsC), layer = "data")
  dotC <- exprC %>% tibble::rownames_to_column("cell") %>%
    mutate(group = grpC[cell]) %>%
    tidyr::pivot_longer(all_of(progC_genes), names_to = "gene", values_to = "expr") %>%
    group_by(group, gene) %>%
    summarise(avg = mean(expr), pct = 100 * mean(expr > 0), .groups = "drop") %>%
    group_by(gene) %>% mutate(avg_scaled = as.numeric(scale(avg))) %>% ungroup() %>%
    mutate(gene = factor(gene, levels = progC_genes),
           group = factor(group, levels = c(ref_cd4, "Cluster 6 (mixed)", ref_cd8)))
  pC <- ggplot(dotC, aes(gene, group)) +
    geom_point(aes(size = pct, color = avg_scaled)) +
    scale_size_continuous(range = c(0.5, 6), name = "% expr") +
    scale_color_gradientn(colors = c("#4575b4", "grey90", "#d73027"), name = "avg (z)") +
    base_c6 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "C  Lineage / Trm-IEL / activation program markers", x = NULL, y = NULL)

  # ---- Panel D: sub-clusters stay lineage-mixed (composition bar) --------------
  pD <- ggplot(subD, aes(x = sub_lab, y = pct, fill = lineage)) +
    geom_col(width = 0.55, color = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(round(pct), "%")),
              position = position_stack(vjust = 0.5), size = 2.8, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("CD8" = "#4E79A7", "CD4" = "#E15759"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    base_c6 + theme(plot.title = element_text(size = 9.5), plot.subtitle = element_text(hjust = 0.5, size = 7)) +
    labs(title = sprintf("D  Splits on activation, not lineage\n(max silhouette = %.2f)", max(sil_c6)),
         subtitle = "each sub-cluster stays CD4/CD8-mixed",
         x = NULL, y = "% resolvable")

  # ---- Panel E: sub-clustering feature plots -----------------------------------
  # Activation (IEG) module score segregates by sub-cluster = the split axis;
  # Trm module score is uniform = shared program; Cd8a/Cd4 scatter = lineage not the axis.
  ieg_set <- list(intersect(c("Atf3", "Fosb", "Dusp1", "Nr4a1", "Fos", "Zfp36"), rownames(cl6)))
  trm_set <- list(intersect(c("Itgae", "Itga1", "Cxcr6", "Id2"), rownames(cl6)))
  cl6 <- AddModuleScore(cl6, features = ieg_set, name = "ActScore", seed = 123)
  cl6 <- AddModuleScore(cl6, features = trm_set, name = "TrmScore", seed = 123)

  cl6_emb <- Embeddings(cl6, "umap.cl6")
  d_df <- data.frame(UMAP_1 = cl6_emb[, 1], UMAP_2 = cl6_emb[, 2],
                     subcluster = factor(cl6$subcluster),
                     row.names = colnames(cl6))
  d_expr <- FetchData(cl6, vars = intersect(c("Cd8a", "Cd4"), rownames(cl6)), layer = "data")

  base_umap_d <- base_c6 + theme(plot.title = element_text(size = 9),
                                 axis.title = element_text(size = 7),
                                 axis.text = element_blank(), axis.ticks = element_blank(),
                                 legend.key.size = unit(0.25, "cm"),
                                 legend.text = element_text(size = 6),
                                 legend.title = element_text(size = 7))
  p_e_sub <- ggplot(d_df, aes(UMAP_1, UMAP_2, color = subcluster)) +
    geom_point(size = 0.5, alpha = 0.85) +
    scale_color_manual(values = c("0" = "#7570B3", "1" = "#D95F02"), name = "sub") +
    base_umap_d + labs(title = "E  Sub-clusters", x = "UMAP1", y = "UMAP2")

  feat_plot <- function(vals, ttl, diverging = FALSE) {
    df <- data.frame(UMAP_1 = d_df$UMAP_1, UMAP_2 = d_df$UMAP_2, v = vals)
    df <- df[order(df$v), ]
    g <- ggplot(df, aes(UMAP_1, UMAP_2, color = v)) +
      geom_point(size = 0.5, alpha = 0.85) + base_umap_d +
      labs(title = ttl, x = "UMAP1", y = "UMAP2")
    if (diverging)
      g + scale_color_gradient2(low = "#4575b4", mid = "grey92", high = "#d73027", midpoint = 0, name = NULL)
    else
      g + scale_color_gradientn(colors = c("grey88", "#FEE0D2", "#FC9272", "#DE2D26"), name = NULL)
  }
  p_e_act <- feat_plot(cl6$ActScore1[rownames(d_df)], "Activation score (IEG)", diverging = TRUE)
  p_e_trm <- feat_plot(cl6$TrmScore1[rownames(d_df)], "Trm score (shared)",     diverging = TRUE)
  p_e_cd8 <- feat_plot(d_expr[rownames(d_df), "Cd8a"], "Cd8a")
  p_e_cd4 <- feat_plot(d_expr[rownames(d_df), "Cd4"],  "Cd4")

  cap_c6 <- paste0(
    "Cluster 6 (781 cells; provisional 'Trm-like CD8') carries a shared Trm/IEL program but contains both CD8 and CD4 T cells. ",
    "(A) Among ", n_resolv, " cells with a detectable lineage marker, ", pct_cd8_raw, "% are CD8 and ", 100 - pct_cd8_raw,
    "% CD4; a further ", n_dn, " (", round(100 * n_dn / length(cells_cl6)), "%) had no detectable Cd4/Cd8a/Cd8b1 (dropout) and ",
    n_dp, " were double-positive, both left unassigned (majority-vote imputation of the dropout cells gives a concordant ",
    pct_cd8_knn, "/", 100 - pct_cd8_knn, ", robustness check only). (B) CD8 and CD4 module scores form two arms. ",
    "(C) The cluster expresses both lineages' markers on a shared Trm/IEL program vs clean CD8 (Teff CD8-like) and CD4 (Th2-like) references. ",
    "(D-E) Lineage is NOT a separable axis: unsupervised sub-clustering (max silhouette ", sprintf("%.2f", max(sil_c6)),
    " at res ", best_res_c6, ") splits on immediate-early/stress genes (Atf3, Fosb, Dusp1, Nr4a1, Fos; sub-cluster 1) not on Cd4/Cd8 — ",
    "the activation module score segregates by sub-cluster while the Trm score is uniform and Cd8a/Cd4 scatter across both, ",
    "and every sub-cluster stays CD4/CD8-mixed. Doublets were ruled out (scDblFinder 1.2%, lowest of any T-cell cluster; Tcell_doublet_check.R). ",
    "The dominant transcriptional axis is a shared Trm/IEL + activation program; CD4/CD8 lineage is a minor axis the data does not resolve. ",
    "The cluster is annotated 'Trm type mixed CD8+ CD4' rather than split, and excluded from CD4/CD8-specific analyses.")

  row1 <- (pA | pB | pC | pD) + patchwork::plot_layout(widths = c(0.6, 1.3, 1.3, 0.9))
  row2 <- (p_e_sub | p_e_act | p_e_trm | p_e_cd8 | p_e_cd4)
  supp_c6 <- (row1 / row2) +
    patchwork::plot_layout(heights = c(1, 0.9)) +
    patchwork::plot_annotation(
      title = "Supplementary — justification of the 'Trm type mixed CD8+ CD4' annotation (cluster 6)",
      caption = stringr::str_wrap(cap_c6, 175),
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                    plot.caption = element_text(hjust = 0, size = 8.5)))
  svglite("../figures/Supp_MixedTrm_annotation.svg", width = 18, height = 10)
  print(supp_c6)
  dev.off()
  cat("Saved: ../figures/Supp_MixedTrm_annotation.svg\n")

  # ---- relabel: provisional "Trm-like CD8" -> confirmed mixed label ------------
  t_cells_TCRab@meta.data[cells_cl6, "celltype_merged"] <- "Trm type mixed CD8+ CD4"
  cat(sprintf("  Relabeled %d cells: 'Trm-like CD8' -> 'Trm type mixed CD8+ CD4'\n", length(cells_cl6)))
} else {
  cat("  Cluster 6 too small —", length(cells_cl6), "cells, skipping mixed-lineage assessment\n")
}

## ---- 4.9-save-annotated -----------------------------------------------------
saveRDS(t_cells_TCRab, "../objects/tcell_TCRab_annotated_final04.rds")
cat("Saved: ../objects/tcell_TCRab_annotated_final04.rds\n")

## ============================================================================
## STAGE 4c — FINAL FIGURES
## ============================================================================

## ---- 4.10-palette -----------------------------------------------------------
# Extend palette for any new cell types (e.g., Th1-like/Tfh-like from split)
all_ct     <- unique(t_cells_TCRab$celltype_merged)
missing_ct <- setdiff(all_ct[!is.na(all_ct)], names(pal_tcrab))
if (length(missing_ct) > 0) {
  extra_cols <- setNames(grDevices::hcl.colors(length(missing_ct), "Dynamic"), missing_ct)
  pal_tcrab  <- c(pal_tcrab, extra_cols)
  cat("Added fallback colors for:", paste(missing_ct, collapse = ", "), "\n")
}
Idents(t_cells_TCRab) <- "celltype_merged"

## ---- 4.11-umap-base ---------------------------------------------------------
umap_emb <- t_cells_TCRab@reductions[[reduction_used]]@cell.embeddings
umap_df  <- data.frame(
  UMAP_1    = umap_emb[, 1],
  UMAP_2    = umap_emb[, 2] * -1,
  celltype  = t_cells_TCRab$celltype_merged,
  condition = t_cells_TCRab$condition,
  row.names = rownames(umap_emb)
)
theme_umap <- theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(), legend.title = element_blank(), aspect.ratio = 1,
        axis.line = element_blank(), axis.ticks = element_blank())

## ---- 4.12-umap-annotated ----------------------------------------------------
label_df <- umap_df %>%
  dplyr::filter(!is.na(celltype)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap_annot <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_label_repel(data = label_df, aes(label = celltype),
                   size = 5, fill = "white", alpha = 0.85,
                   show.legend = FALSE, box.padding = 0.4, max.overlaps = 20) +
  scale_color_manual(values = pal_tcrab, na.translate = FALSE) +
  theme_umap + labs(title = NULL) +
  theme(legend.position = "none")

svglite("../figures/UMAP_TCRab_annotated.svg", width = 7, height = 5)
print(p_umap_annot)
dev.off()
cat("Saved: ../figures/UMAP_TCRab_annotated.svg\n")

## ---- 4.13-umap-by-condition -------------------------------------------------
cond_levels <- c("BAS", "PD", "REC", "RCU")
conditions  <- intersect(cond_levels, unique(umap_df$condition))
p_cond_list <- lapply(conditions, function(cond) {
  fg_cond <- umap_df[umap_df$condition == cond & !is.na(umap_df$celltype), ]
  ggplot(fg_cond, aes(UMAP_1, UMAP_2, color = celltype)) +
    geom_point(size = 0.45, alpha = 0.9) +
    scale_color_manual(values = pal_tcrab, na.translate = FALSE) +
    theme_umap + labs(title = cond) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
})
svglite("../figures/UMAP_TCRab_byCondition.svg", width = 18, height = 5)
print(patchwork::wrap_plots(p_cond_list, nrow = 1))
dev.off()
cat("Saved: ../figures/UMAP_TCRab_byCondition.svg\n")

## ---- 4.14-pie-charts --------------------------------------------------------
pie_df <- t_cells_TCRab@meta.data %>%
  dplyr::filter(!is.na(celltype_merged)) %>%
  dplyr::count(condition, celltype_merged, name = "n") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(prop = n / sum(n),
                label = ifelse(prop >= 0.02, paste0(round(prop * 100, 1), "%"), "")) %>%
  dplyr::ungroup()
pie_df$condition <- factor(pie_df$condition, levels = cond_levels)

p_pie_list <- lapply(conditions, function(cond) {
  ggplot(pie_df[pie_df$condition == cond, ], aes(x = 1, y = prop, fill = celltype_merged)) +
    geom_col(width = 1, color = NA) +
    coord_polar(theta = "y") +
    geom_text(aes(x = 1.65, label = label),
              position = position_stack(vjust = 0.5), size = 3.9, color = "black",
              fontface = "bold") +
    xlim(c(0, 2.1)) +
    scale_fill_manual(values = pal_tcrab) +
    theme_void(base_size = 11) +
    labs(title = cond) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.title = element_blank(), legend.text = element_text(size = 12))
})
p_pies <- patchwork::wrap_plots(p_pie_list, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 15, face = "bold")) &
  guides(fill = guide_legend(nrow = 1))
svglite("../figures/Pie_TCRab_byCondition.svg", width = 14, height = 6)
print(p_pies)
dev.off()
cat("Saved: ../figures/Pie_TCRab_byCondition.svg\n")

## ---- 4.15-heatmap-final -----------------------------------------------------
# Cell types as ROWS, genes as COLUMNS
# Row annotation: n_cells bar + Celltype color strip (matching UMAP)
Idents(t_cells_TCRab) <- "celltype_merged"

hm_genes_flat <- unique(unlist(modules_hm))
avg_exp   <- AverageExpression(t_cells_TCRab, assays = "RNA", slot = "data",
                                group.by = "celltype_merged")[[1]]
kept_genes <- hm_genes_flat[hm_genes_flat %in% rownames(avg_exp)]
mat_hm     <- avg_exp[kept_genes, ]

z_hm           <- t(scale(t(mat_hm)))
z_hm[is.na(z_hm)] <- 0; z_hm[z_hm > 2] <- 2; z_hm[z_hm < -1.5] <- -1.5

ct_levels  <- colnames(z_hm)
size_map   <- table(t_cells_TCRab$celltype_merged)

row_ann_final <- data.frame(
  n_cells  = as.integer(size_map[ct_levels]),
  Celltype = factor(ct_levels, levels = ct_levels),
  row.names = ct_levels, check.names = FALSE
)
ann_colors_final <- list(
  n_cells  = colorRampPalette(c("white", "steelblue4"))(100),
  Celltype = pal_tcrab[ct_levels]
)
ann_colors_final$Celltype <- ann_colors_final$Celltype[!is.na(ann_colors_final$Celltype)]

ph_final <- pheatmap::pheatmap(
  mat                      = z_hm,
  cluster_rows             = FALSE,
  cluster_cols             = TRUE,
  clustering_distance_cols = "correlation",
  clustering_method        = "ward.D2",
  show_rownames            = TRUE,
  show_colnames            = TRUE,
  border_color             = NA,
  annotation_col           = row_ann_final,
  annotation_colors        = ann_colors_final,
  annotation_legend        = FALSE,
  annotation_names_col     = FALSE,
  legend                   = FALSE,          # suppress built-in z-score scale; drawn manually below
  fontsize                 = 13,
  fontsize_row             = 10,
  fontsize_col             = 13,
  angle_col                = 45,
  breaks                   = seq(-1.5, 2, length.out = 101),
  color                    = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100),
  main                     = "αβ T cell subtypes (avg exp, z-scored by gene)",
  silent                   = TRUE
)

# Post-process: make all text grobs bold
make_bold_grobs <- function(g) {
  if (inherits(g, "text")) {
    gp_list <- if (is.null(g$gp)) list() else as.list(g$gp)
    gp_list[["font"]]     <- NULL
    gp_list[["fontface"]] <- "bold"
    g$gp <- do.call(grid::gpar, gp_list)
    return(g)
  }
  if (inherits(g, "gTree") && !is.null(g$children)) {
    g$children <- do.call(grid::gList, lapply(g$children, make_bold_grobs))
  }
  g
}
gt_bold <- ph_final$gtable
gt_bold$grobs <- lapply(gt_bold$grobs, make_bold_grobs)

# Un-bold gene labels (row_names), size 11
make_plain_size11 <- function(g) {
  if (inherits(g, "text")) {
    gp_list <- if (is.null(g$gp)) list() else as.list(g$gp)
    gp_list[["font"]]     <- NULL
    gp_list[["fontface"]] <- "plain"
    gp_list[["fontsize"]] <- 11
    g$gp <- do.call(grid::gpar, gp_list)
    return(g)
  }
  if (inherits(g, "gTree") && !is.null(g$children)) {
    g$children <- do.call(grid::gList, lapply(g$children, make_plain_size11))
  }
  g
}
rn_idx <- which(gt_bold$layout$name == "row_names")
if (length(rn_idx) > 0)
  gt_bold$grobs[[rn_idx]] <- make_plain_size11(gt_bold$grobs[[rn_idx]])

# ---- Manual legends: z-score (left) + n_cells (right), side by side -----------
hm_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(100)
n_vals    <- as.integer(size_map[ct_levels])
n_min     <- min(n_vals);  n_max <- max(n_vals)
n_colors  <- colorRampPalette(c("white", "steelblue4"))(100)

# z-score gradient: top = 2 (red), bottom = -1.5 (blue)
zscore_grob <- grid::grobTree(
  grid::rasterGrob(matrix(rev(hm_colors), ncol = 1),
                   x = unit(0.18, "npc"), y = unit(0.50, "npc"),
                   width = unit(0.35, "npc"), height = unit(0.82, "npc"),
                   just = c("left", "center")),
  grid::textGrob("2",    x=unit(0.58,"npc"), y=unit(0.97,"npc"), just=c("left","center"), gp=grid::gpar(fontsize=11, fontface="bold")),
  grid::textGrob("1",    x=unit(0.58,"npc"), y=unit(0.82,"npc"), just=c("left","center"), gp=grid::gpar(fontsize=11, fontface="bold")),
  grid::textGrob("0",    x=unit(0.58,"npc"), y=unit(0.57,"npc"), just=c("left","center"), gp=grid::gpar(fontsize=11, fontface="bold")),
  grid::textGrob("-1",   x=unit(0.58,"npc"), y=unit(0.32,"npc"), just=c("left","center"), gp=grid::gpar(fontsize=11, fontface="bold")),
  grid::textGrob("-1.5", x=unit(0.58,"npc"), y=unit(0.09,"npc"), just=c("left","center"), gp=grid::gpar(fontsize=11, fontface="bold"))
)

# n_cells gradient: top = n_max (dark blue), bottom = n_min (white)
ncells_grob <- grid::grobTree(
  grid::textGrob("n_cells",
                 x=unit(0.05,"npc"), y=unit(0.98,"npc"), just=c("left","top"),
                 gp=grid::gpar(fontsize=12, fontface="bold")),
  grid::rasterGrob(matrix(rev(n_colors), ncol = 1),
                   x = unit(0.05, "npc"), y = unit(0.50, "npc"),
                   width = unit(0.35, "npc"), height = unit(0.74, "npc"),
                   just = c("left", "center")),
  grid::textGrob(formatC(n_max, format="d", big.mark=","),
                 x=unit(0.45,"npc"), y=unit(0.88,"npc"), just=c("left","center"),
                 gp=grid::gpar(fontsize=11, fontface="bold")),
  grid::textGrob(formatC(n_min, format="d", big.mark=","),
                 x=unit(0.45,"npc"), y=unit(0.12,"npc"), just=c("left","center"),
                 gp=grid::gpar(fontsize=11, fontface="bold"))
)

# ---- Canvas: 9.5 × 22 in ------------------------------------------------------
# Pheatmap (no legend) in left 8 in × full 22 in
# Both legend bars side-by-side in right 1.5 in, vertically centred on matrix area
# Matrix occupies roughly y=3 in to y=18 in from bottom in a 22 in canvas
leg_h     <- 5    # inches tall — compact publication legend
leg_y_bot <- 8.5  # center on the canvas (22/2 - 5/2 = 8.5)

svglite("../figures/heatmap_TCRab.svg", width = 9.5, height = 22)
grid::grid.newpage()

# Pheatmap — left 8/9.5 of canvas, full height
grid::pushViewport(grid::viewport(
  x = 0, y = 0, width = 8/9.5, height = 1, just = c("left", "bottom")
))
grid::grid.draw(gt_bold)
grid::popViewport()

# z-score legend — first 0.7 in of the right strip
grid::pushViewport(grid::viewport(
  x = unit(8.1, "in"), y = unit(leg_y_bot, "in"),
  width = unit(0.65, "in"), height = unit(leg_h, "in"),
  just = c("left", "bottom")
))
grid::grid.draw(zscore_grob)
grid::popViewport()

# n_cells legend — next 0.75 in of the right strip
grid::pushViewport(grid::viewport(
  x = unit(8.8, "in"), y = unit(leg_y_bot, "in"),
  width = unit(0.65, "in"), height = unit(leg_h, "in"),
  just = c("left", "bottom")
))
grid::grid.draw(ncells_grob)
grid::popViewport()

dev.off()
cat("Saved: ../figures/heatmap_TCRab.svg\n")

## ---- 4.16-violins -----------------------------------------------------------
# Curated marker set — organized by function, ordered for stacked violin readability:
#   Lineage: Cd4, Cd8a
#   Naive/Central memory: Ccr7, Sell, Tcf7, Bcl2
#   Th2: Gata3
#   Th17: Rorc, Il17a, Il17f
#   Tfh/Tph: Bcl6, Il21, Icos, Cd40lg, Pdcd1
#   Treg: Foxp3, Il2ra, Ctla4
#   TRM: Itgae, Cxcr6, Id2
#   Effector/cytotoxic: Ccl5, Gzmk, Gzmb, Ifng
#   Exhaustion: Havcr2
#   Proliferation: Mki67
vln_genes_raw <- c("Cd4","Cd8a",
                   "Ccr7","Sell","Tcf7","Bcl2",
                   "Gata3",
                   "Rorc","Il17a","Il17f",
                   "Bcl6","Il21","Icos","Cd40lg","Pdcd1",
                   "Foxp3","Il2ra","Ctla4",
                   "Itgae","Cxcr6","Id2",
                   "Ccl5","Gzmk","Gzmb","Ifng",
                   "Havcr2",
                   "Mki67")
vln_genes  <- sapply(vln_genes_raw, resolve_gene, rownames(t_cells_TCRab))
vln_genes  <- vln_genes[!is.na(vln_genes)]
ct_ordered <- names(pal_tcrab)[names(pal_tcrab) %in% unique(t_cells_TCRab$celltype_merged)]

# Slightly darker version of each cell-type color for jittered dots
darken_col <- function(hex, amount = 0.25) {
  m <- col2rgb(hex) / 255
  rgb(pmax(0, m[1] - amount), pmax(0, m[2] - amount), pmax(0, m[3] - amount))
}
pal_tcrab_dark <- setNames(sapply(pal_tcrab, darken_col), names(pal_tcrab))

vln_list <- lapply(vln_genes, function(gn) {
  expr_df <- data.frame(
    expr     = FetchData(t_cells_TCRab, vars = gn, layer = "data")[, 1],
    celltype = t_cells_TCRab$celltype_merged
  )
  expr_df <- expr_df[!is.na(expr_df$celltype), ]
  expr_df$celltype <- factor(expr_df$celltype, levels = ct_ordered)

  ggplot(expr_df, aes(x = celltype, y = expr, fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75, color = NA) +
    ggrastr::geom_jitter_rast(aes(color = celltype),
                              width = 0.25, size = 0.4, alpha = 0.45, shape = 16,
                              raster.dpi = 150) +
    scale_fill_manual(values  = pal_tcrab, name = "Cell type") +
    scale_color_manual(values = pal_tcrab_dark, guide = "none") +
    theme_classic(base_size = 10) +
    theme(panel.grid      = element_blank(),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          axis.title.x    = element_blank(),
          plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
          axis.title.y    = element_text(size = 9),
          axis.text.y     = element_text(size = 8),
          legend.position = "right") +
    labs(y = "Expression", title = gn)
})

# Collect shared legend from first plot; all others drop it
vln_wrap <- patchwork::wrap_plots(vln_list, ncol = 3) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position  = "right",
        legend.text      = element_text(size = 11, face = "bold"),
        legend.title     = element_text(size = 12, face = "bold"),
        legend.key.size  = unit(0.55, "cm"))

svglite("../figures/violin_markers.svg", width = 15, height = 22)
print(vln_wrap)
dev.off()
cat("Saved: ../figures/violin_markers.svg\n")

## ---- 4.17-stacked-violin ----------------------------------------------------
# Combined stacked violin using the same vln_genes as violin_markers.svg above
# (Gzmb, Prf1, Nkg7, Tbx21, Ifng, Rorc, Il17a, Foxp3, Il2ra, Bcl6, Cxcr5,
#  Batf, Icos, Maf, Itgae, Cd69, Pdcd1, Lag3, Ccr7, Sell, Lef1, Mki67, Zbtb16, Trav1)
# Genes as rows (facet_grid, labels on left), all T cell types on x-axis.
stk_expr <- FetchData(t_cells_TCRab, vars = unname(vln_genes), layer = "data")
stk_expr$celltype_merged <- t_cells_TCRab$celltype_merged[rownames(stk_expr)]

stk_long <- stk_expr %>%
  dplyr::filter(!is.na(celltype_merged), celltype_merged %in% ct_ordered) %>%
  dplyr::mutate(celltype_merged = factor(celltype_merged, levels = ct_ordered)) %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(cols = all_of(unname(vln_genes)), names_to = "gene", values_to = "expr") %>%
  dplyr::mutate(gene = factor(gene, levels = rev(unname(vln_genes))))

p_stk_vln <- ggplot(stk_long, aes(x = celltype_merged, y = expr, fill = celltype_merged)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, color = NA) +
  facet_grid(gene ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = pal_tcrab, guide = "none") +
  theme_classic(base_size = 10) +
  theme(
    panel.grid        = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9, face = "bold.italic"),
    strip.placement   = "outside",
    strip.background  = element_blank(),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    axis.title        = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 9),
    panel.spacing.y   = unit(0.05, "lines")
  )

svglite("../figures/violin_markers_stacked.svg", width = 8, height = 14)
print(p_stk_vln)
dev.off()
cat("Saved: ../figures/violin_markers_stacked.svg\n")

## ---- done -------------------------------------------------------------------
cat("\n=== Module 03 complete ===\n")
cat("Objects: tcell_TCRab_clustered04a.rds, tcell_TCRab_annotated_final04.rds,\n")
cat("         cl7_innate_MAIT.rds, cl2_Tph_subcluster.rds\n")
cat("Figures: UMAP_TCRab_numbered.svg, UMAP_TCRab_numbered_byCondition.svg,\n")
cat("         heatmap_TCRab_numbered.svg, UMAP_TCRab_annotated.svg,\n")
cat("         UMAP_TCRab_byCondition.svg, Pie_TCRab_byCondition.svg,\n")
cat("         heatmap_TCRab.svg, violin_markers.svg\n")

if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
writeLines(capture.output(sessionInfo()), "../results/sessionInfo_stage4.txt")
