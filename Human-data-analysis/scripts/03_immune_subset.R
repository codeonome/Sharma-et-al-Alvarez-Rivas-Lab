# =====================================================================
# Step 03 -- Immune compartment: subset, re-integrate, recluster, and
#            select T cells
# Input  : objects/02_global_clustered.rds
# Output : objects/03_immune_subset.rds    (48,346 cells; immune_res.<r>)
#          objects/03_tcell_subset.rds     (21,899 cells)
#          results/03_immune_subset/*      global + immune markers
#          figures/03_immune_subset/*
#
# Global res-0.1 annotation (12 clusters; markers in cluster_markers_res0.1.csv):
#   0 vascular endothelium | 1 T/NK | 2 fibroblasts | 3 plasma | 4 myeloid
#   5 mural | 6 epithelium | 7 mast | 8 B | 9 lymphatic endothelium
#   10 cycling (mixed) | 11 neural crest (mixed)
#   -> immune clusters carried forward: 1, 3, 4, 7, 8
#
# Immune res-0.2 annotation (13 clusters):
#   0 CD4 T | 1 CD8 T | 2 plasma A | 3 mono/macrophage | 4 plasma B (CD38hi)
#   5 B | 6 mast | 7 NK (+ gd T) | 8 pDC | 9 endothelial contaminant
#   10 fibroblast contaminant | 11 cDC1 | 12 plasma C (lambda, n = 81)
#   -> T cells: clusters 0 and 1
#
# Reproducibility note: there is no set.seed() before the immune Harmony
# call. Its random state is set by RunPCA() and RunUMAP() (seed.use = 42)
# immediately before it, so those calls must stay in this order.
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(harmony)
})
source("config.R"); source("R/helpers.R")
STEP <- "03_immune_subset"
d <- step_dirs(STEP)
dir.create(file.path(d$figures, "resolution_sweep"), showWarnings = FALSE)

combined <- readRDS(OBJ$global)
global_col <- paste0("global_res.", GLOBAL_RES_FINAL)
Idents(combined) <- global_col
combined$seurat_clusters <- combined[[global_col]][, 1]

# ---------------------------------------------------------------------
# 1. Global markers at the chosen resolution
# ---------------------------------------------------------------------
markers_and_top10(combined,
                  file.path(d$results, sprintf("cluster_markers_res%s.csv", GLOBAL_RES_FINAL)),
                  file.path(d$results, sprintf("top10_markers_res%s.csv", GLOBAL_RES_FINAL)))

# T-cell marker detection in every global cluster (checks for stray T cells
# outside the immune selection)
tcell_markers <- c("CD3D", "CD3E", "CD3G", "TRAC")
df_t <- FetchData(combined, vars = c(tcell_markers, global_col))
colnames(df_t)[ncol(df_t)] <- "cluster"
p_stray <- ggplot(pivot_longer(df_t, all_of(tcell_markers), names_to = "gene", values_to = "expr"),
                  aes(x = cluster, y = expr)) +
  geom_violin(fill = "steelblue", scale = "width", trim = TRUE) +
  facet_wrap(~gene, ncol = 2, scales = "free_y") + theme_classic() +
  ggtitle(sprintf("T cell markers by cluster (res %s)", GLOBAL_RES_FINAL))
ggsave(file.path(d$figures, sprintf("tcell_markers_all_clusters_res%s.png", GLOBAL_RES_FINAL)),
       p_stray, width = 10, height = 7, dpi = 150)

# ---------------------------------------------------------------------
# 2. Subset to immune clusters and reprocess within the compartment
# ---------------------------------------------------------------------
immune_obj <- subset(combined, cells = WhichCells(combined, idents = IMMUNE_CLUSTERS))
rm(combined); gc()
log_step("Immune cells: ", ncol(immune_obj))

DefaultAssay(immune_obj) <- "RNA"
immune_obj[["RNA"]] <- split(immune_obj[["RNA"]], f = immune_obj$sample)
immune_obj <- NormalizeData(immune_obj)
immune_obj <- FindVariableFeatures(immune_obj, nfeatures = N_VARIABLE_FEATURES)
immune_obj <- ScaleData(immune_obj)
immune_obj <- RunPCA(immune_obj, npcs = IMMUNE_NPCS, seed.use = SEED_PCA_UMAP)
plot_elbow(immune_obj, "Immune subset -- elbow plot (PCA, pre-integration)",
           file.path(d$figures, "elbow_immune_pca.svg"))

immune_obj <- RunUMAP(immune_obj, reduction = "pca", dims = IMMUNE_DIMS,
                      reduction.name = "umap.unintegrated.immune", seed.use = SEED_PCA_UMAP)

# ---------------------------------------------------------------------
# 3. Harmony integration (batch = sample) and UMAP
# ---------------------------------------------------------------------
immune_obj <- IntegrateLayers(immune_obj, method = HarmonyIntegration, orig.reduction = "pca",
                              new.reduction = "integrated.harmony.immune")
immune_obj <- JoinLayers(immune_obj)
set.seed(SEED_UMAP_IMM)
immune_obj <- RunUMAP(immune_obj, reduction = "integrated.harmony.immune", dims = IMMUNE_DIMS,
                      reduction.name = "umap.immune", seed.use = SEED_PCA_UMAP)

for (red in c("umap.unintegrated.immune", "umap.immune")) {
  tag <- if (red == "umap.immune") "harmony" else "unintegrated"
  plot_umap(immune_obj, red, "sample", paste("Immune subset --", tag, "by sample"),
            file.path(d$figures, paste0("umap_immune_", tag, "_by_sample")),
            facet = TRUE, facet_ncol = 5, width = 16, height = 16)
  plot_umap(immune_obj, red, "condition", paste("Immune subset --", tag, "by condition"),
            file.path(d$figures, paste0("umap_immune_", tag, "_by_condition")),
            facet = TRUE, facet_ncol = 2, width = 9, height = 5)
}

# ---------------------------------------------------------------------
# 4. Resolution sweep within the immune compartment
# ---------------------------------------------------------------------
immune_obj <- FindNeighbors(immune_obj, reduction = "integrated.harmony.immune", dims = IMMUNE_DIMS)
immune_obj <- FindClusters(immune_obj, resolution = IMMUNE_RES_SWEEP, random.seed = SEED_CLUSTERS)
immune_obj <- rename_res_cols(immune_obj, "immune_res.")
for (cc in res_cols(immune_obj, "immune_res.")) {
  r <- sub("immune_res.", "", cc, fixed = TRUE)
  plot_umap(immune_obj, "umap.immune", cc,
            sprintf("Immune subset -- res %s (%d clusters)", r, length(unique(immune_obj[[cc]][, 1]))),
            file.path(d$figures, "resolution_sweep", paste0("umap_immune_res", r)))
}

immune_obj$immune_clusters <- immune_obj[[paste0("immune_res.", IMMUNE_RES_FINAL)]][, 1]
Idents(immune_obj) <- "immune_clusters"

# ---------------------------------------------------------------------
# 5. Immune markers and canonical lineage panels
# ---------------------------------------------------------------------
markers_and_top10(immune_obj,
                  file.path(d$results, sprintf("immune_subset_markers_res%s.csv", IMMUNE_RES_FINAL)),
                  file.path(d$results, sprintf("immune_subset_top10_markers_res%s.csv", IMMUNE_RES_FINAL)))

marker_panel <- list(
  Tcell   = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "LCK", "ITK"),
  NK      = c("GNLY", "NKG7", "KLRD1", "KLRF1", "PRF1"),
  Bcell   = c("MS4A1", "CD79A", "PAX5"),
  Myeloid = c("LYZ", "C1QA", "FCN1"),
  Mast    = c("TPSAB1", "CPA3"),
  Plasma  = c("MZB1", "DERL3"),
  pDC     = c("LILRA4", "CLEC4C")
)
avail <- lapply(marker_panel, function(g) intersect(g, rownames(immune_obj)))
for (grp in names(avail)) {
  if (length(avail[[grp]]) == 0) next
  plot_marker_panel(immune_obj, "umap.immune", avail[[grp]], paste(grp, "markers"),
                    file.path(d$figures, paste0("immune_subset_markers_", grp)))
}
avg_expr <- AverageExpression(immune_obj, features = unlist(avail, use.names = FALSE),
                              group.by = "immune_clusters", layer = "data")$RNA
write.csv(round(avg_expr, 2), file.path(d$results, "immune_subset_marker_avg_expr.csv"))

saveRDS(immune_obj, OBJ$immune)

# ---------------------------------------------------------------------
# 6. T cells = immune clusters 0 (CD4) and 1 (CD8)
# ---------------------------------------------------------------------
tcell_obj <- subset(immune_obj, idents = TCELL_CLUSTERS)
log_step(sprintf("T cells: %d / %d immune cells (%.1f%%)", ncol(tcell_obj), ncol(immune_obj),
                 100 * ncol(tcell_obj) / ncol(immune_obj)))
saveRDS(tcell_obj, OBJ$tcell)
save_session(STEP, d)
log_step("Done: ", OBJ$immune, ", ", OBJ$tcell)
