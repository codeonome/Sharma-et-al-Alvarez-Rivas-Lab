#' ---
#' title: "Tcell_analysis_v2 — Stages 1-3"
#' author: "Apoorva Sharma"
#' ---
#'
#' **T Cell Analysis — Stages 1–3**
#'   - Extract T cells from global object
#'   - Initial clustering for contamination removal
#'   - TCR status classification → keep only αβ TCR cells
#'   - Output: tcell_TCRab_only03.rds → hand off to Tcell_ab_analysis_reclustering.R
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
  library(scales)
})

## ---- stage1-load ------------------------------------------------------------
# 1. Load global object and subset T cells
rna_tcr <- readRDS("../../01_Global_Analysis/objects/global_clustered_annotated.rds")
cat("Major cell types in global object:", levels(rna_tcr@meta.data$major_celltype), "\n")

t_cells <- subset(rna_tcr, subset = major_celltype == "T")
cat("T cells after subsetting:", ncol(t_cells), "\n")

saveRDS(t_cells, "../objects/tcell_only_unclustered01.RDS")
cat("Saved: ../objects/tcell_only_unclustered01.RDS\n")

## ---- stage2-preprocess ------------------------------------------------------
# 2. Re-processing — Normalize, Variable Features, Scale, PCA
set.seed(123)
DefaultAssay(t_cells) <- "RNA"
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 3000)
t_cells <- ScaleData(t_cells, vars.to.regress = "percent.mt", verbose = FALSE)
t_cells <- RunPCA(t_cells, npcs = 50, verbose = FALSE, seed.use = 123)

p_elbow <- ElbowPlot(t_cells, ndims = 50) +
  ggtitle("Elbow plot — T cells (initial PCA, 50 PCs)") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank())

svglite("../figures/elbow_plot_tcell.svg", width = 8, height = 5)
print(p_elbow)
dev.off()
cat("Saved: ../figures/elbow_plot_tcell.svg\n")

## ---- stage2-cluster ---------------------------------------------------------
# 3. Clustering — 20 PCs (break ~15, extra PCs to separate contaminants)
set.seed(123)
t_cells <- FindNeighbors(t_cells, reduction = "pca", dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.3, cluster.name = "pca_clusters_tcells_0.3")
t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:20,
                   reduction.name = "umap.pca_tcell", seed.use = 123)
Idents(t_cells) <- "pca_clusters_tcells_0.3"

tcell_markers <- FindAllMarkers(
  t_cells, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox"
)
top500 <- tcell_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC), p_val_adj, dplyr::desc(pct.1 - pct.2), .by_group = TRUE) %>%
  dplyr::slice_head(n = 500) %>%
  dplyr::ungroup()

write.csv(top500, "../objects/tcell_top500_markers_initial.csv", row.names = FALSE)
cat("Saved: ../objects/tcell_top500_markers_initial.csv\n")

saveRDS(t_cells, "../objects/tcell_only_clustered_unfilter02.rds")
cat("Saved: ../objects/tcell_only_clustered_unfilter02.rds\n")

## ---- initial-cluster-annotations --------------------------------------------
# Manual review of tcell_top500_markers_initial.csv produced these provisional identities
# (res 0.3 clustering on all T cells, 20 PCs). These cluster numbers are only used
# for contamination removal here — they do NOT map to the final αβ T cell clusters
# assigned in Tcell_ab_analysis_reclustering.R.
#
#   Cluster  Provisional identity        Top markers                  Decision
#   -------  --------------------------  ---------------------------  --------
#   0        CD4+ T (Naive/Th2-like)     Il7r, Sell, Tcf7             KEEP
#   1        CD8+ T effector             Gzmk, Nkg7, Ccl5             KEEP
#   2        CD4+ T (Th1/Tph-like)       Ifng, Pdcd1, Cxcr6           KEEP
#   3        CD4+ Tregs                  Foxp3, Il2ra, Ctla4           KEEP
#   4        CD8+ T effector             Gzmk, Nkg7, Ccl5             KEEP
#   5        Naive T                     Sell, Tcf7, Lef1              KEEP
#   6        NK / ILC1-like              Nkg7, Klrd1, Gzmb             REMOVE
#   7        Neutrophil / PMN            S100a8, S100a9, Mpo           REMOVE
#   8        CD8+ Trm                    Cd69, Itgae, Cxcr6            KEEP
#   9        ILC2-like / or Th2 ambiguous       Il1rl1, Gata3; low/no TCR     KEEP*
#   10       cPMN / granulocyte          Csf3r, S100a8, Lcn2           REMOVE
#   11       Proliferating T             Mki67, Top2a, Pcna            KEEP
#   12       Germinal center B           Aicda, Bcl6, Cd19             REMOVE
#   13       Activated B                 Cd79a, Cd19, Ighg              REMOVE
#   14       Pro-B / pre-B               Igll1, Vpreb1, Cd79a          REMOVE
#   15       CD8+ T (Trm/Exh-like)       Cd69, Lag3, Havcr2            KEEP
#
# * Cluster 9 note: TCR pairing check (below) showed ~95.8% of cluster 9 cells
#   have No TCR — consistent with ILC2s or innate-like lymphocytes co-clustering
#   with genuine Th2 T cells. Removing cluster 9 entirely would discard the ~4%
#   that carry a real αβ TCR. Instead, we keep cluster 9 here and rely on the
#   αβ TCR subsetting step (stage 3) to exclude the non-T fraction automatically.

## ---- tcr-classification -----------------------------------------------------
classify_tcr_status <- function(meta) {
  use_cols <- intersect(c("CTgene","CTaa","CTstrict","chains","c_call","tcr_chains"), colnames(meta))
  txt <- if (length(use_cols) == 0) rep("", nrow(meta)) else
    apply(meta[, use_cols, drop = FALSE], 1, function(r) { r[is.na(r)] <- ""; paste(r, collapse = ";") })

  has_a  <- grepl("TRAV|TRAJ|TRAC|\\bTRA\\b", txt, TRUE)
  has_b  <- grepl("TRBV|TRBJ|TRBC|\\bTRB\\b", txt, TRUE)
  has_gd <- grepl("TRGV|TRGJ|\\bTRG\\b|TRDV|TRDJ|\\bTRD\\b", txt, TRUE)

  status <- ifelse(has_gd & !(has_a | has_b), "γδ (TRG/TRD)",
            ifelse(has_a & has_b,              "Paired (TRA+TRB)",
            ifelse(has_a != has_b,             "Single chain (TRA or TRB)",
                                               "No TCR")))
  factor(status, levels = c("Paired (TRA+TRB)", "Single chain (TRA or TRB)", "γδ (TRG/TRD)", "No TCR"))
}

## ---- cluster9-check ---------------------------------------------------------
# 4. Cluster 9 TCR pairing check before contamination removal
#    Th2 cells carry αβ TCR; ILC2s do not — check resolves the ambiguity.
t_cells$tcr_status_all <- classify_tcr_status(t_cells@meta.data)

cluster9_cells <- t_cells@meta.data[t_cells$pca_clusters_tcells_0.3 == "9", ]
cluster9_tcr   <- as.data.frame(table(TCR_Status = cluster9_cells$tcr_status_all)) %>%
  dplyr::mutate(Percentage = round(100 * prop.table(Freq), 2), Cluster = "9")
cat("\n--- Cluster 9 TCR pairing status ---\n")
print(cluster9_tcr)
write.csv(cluster9_tcr, "../objects/cluster9_TCR_pairing_check.csv", row.names = FALSE)
cat("Saved: ../objects/cluster9_TCR_pairing_check.csv\n\n")

## ---- contamination-removal --------------------------------------------------
# 5. Remove contamination clusters
#    v1 removed: 6 (NK/ILC1-like), 7 (PMN), 10 (cPMN), 12 (germ B), 13 (act B), 14 (pro B)
#    Cluster 9 NOT removed, as annotation was not clear —  upon checking-most of them were ICl2, αβ TCR subset below excludes ~96% of these ILC2-like cells
remove_clusters <- c("6", "7", "10", "12", "13", "14")
t_cells_filt <- subset(t_cells, idents = setdiff(levels(Idents(t_cells)), remove_clusters))
Idents(t_cells_filt) <- droplevels(Idents(t_cells_filt))
cat("Remaining T cell clusters:", paste(levels(Idents(t_cells_filt)), collapse = ", "), "\n")
cat("Remaining T cells:", ncol(t_cells_filt), "\n")

## ---- cell-count-bar ---------------------------------------------------------
if (!"condition" %in% colnames(t_cells_filt@meta.data))
  t_cells_filt$condition <- factor(sub("_.*", "", t_cells_filt$orig.ident), levels = c("BAS","PD","REC","RCU"))
t_cells_filt$condition <- factor(t_cells_filt$condition, levels = c("BAS","PD","REC","RCU"))

cluster_counts <- t_cells_filt@meta.data %>%
  dplyr::count(pca_clusters_tcells_0.3, condition, name = "n_cells") %>%
  dplyr::mutate(cluster = factor(pca_clusters_tcells_0.3,
                                  levels = sort(as.numeric(as.character(unique(pca_clusters_tcells_0.3))))))

p_cellcount <- ggplot(cluster_counts, aes(x = cluster, y = n_cells, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("BAS"="#4E79A7","PD"="#F28E2B","REC"="#59A14F","RCU"="#E15759")) +
  labs(title = "Cell counts per cluster per condition", x = "Cluster (res 0.3)", y = "N cells", fill = "Condition") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))

svglite("../figures/cell_counts_per_cluster_per_condition.svg", width = 10, height = 5)
print(p_cellcount)
dev.off()
cat("Saved: ../figures/cell_counts_per_cluster_per_condition.svg\n")

## ---- tcr-status-summary -----------------------------------------------------
# 6. TCR status for filtered T cells
t_cells_filt$tcr_status <- classify_tcr_status(t_cells_filt@meta.data)

overall_tcr <- as.data.frame(table(TCR_Status = t_cells_filt$tcr_status)) %>%
  dplyr::mutate(Percentage = round(100 * prop.table(Freq), 2))
write_csv(overall_tcr, "../objects/tcr_overall_filtered.csv")
cat("Saved: ../objects/tcr_overall_filtered.csv\n")

if (!"condition" %in% colnames(t_cells_filt@meta.data))
  t_cells_filt$condition <- factor(sub("_.*", "", t_cells_filt$orig.ident), levels = c("BAS","PD","REC","RCU"))
t_cells_filt$condition <- factor(t_cells_filt$condition, levels = c("BAS","PD","REC","RCU"))

per_condition <- t_cells_filt@meta.data %>%
  dplyr::count(condition, tcr_status, name = "Count") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(Percentage = round(100 * Count / sum(Count), 2)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(condition, tcr_status)

## ---- tcr-bar-chart ----------------------------------------------------------
pal_tcr <- c(
  "Paired (TRA+TRB)"          = "#5c7fa0",
  "Single chain (TRA or TRB)" = "#8eafc7",
  "γδ (TRG/TRD)"              = "#b8cdd8",
  "No TCR"                    = "#e5eef3"
)

p_tcr_bar <- ggplot(per_condition, aes(x = condition, y = Count, fill = tcr_status)) +
  geom_bar(stat = "identity", color = "white", position = "fill") +
  geom_text(aes(label = ifelse(Percentage >= 2.0, paste0(Percentage, "%"), "")),
            position = position_fill(vjust = 0.5), color = "black", size = 3.5) +
  scale_fill_manual(values = pal_tcr, name = "TCR class", drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "TCR chain pairing per condition", x = "Condition", y = "Proportion") +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))

svglite("../figures/TCR_status_by_condition.svg", width = 7, height = 5)
print(p_tcr_bar)
dev.off()
cat("Saved: ../figures/TCR_status_by_condition.svg\n")

## ---- subset-tcrab -----------------------------------------------------------
# 7. Keep only αβ TCR cells (Paired + Single chain)
keep_tcr   <- c("Paired (TRA+TRB)", "Single chain (TRA or TRB)")
keep_cells <- Cells(t_cells_filt)[t_cells_filt$tcr_status %in% keep_tcr]
t_cells_TCRab <- subset(t_cells_filt, cells = keep_cells)
cat("Kept αβ TCR cells:", length(keep_cells), "/", ncol(t_cells_filt), "\n")

saveRDS(t_cells_TCRab, "../objects/tcell_TCRab_only03.rds")
cat("Saved: ../objects/tcell_TCRab_only03.rds\n")

cat("\n=== Stages 1-3 complete. Run Tcell_ab_analysis_reclustering.R next. ===\n")

if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
writeLines(capture.output(sessionInfo()), "../results/sessionInfo_stages1to3.txt")
