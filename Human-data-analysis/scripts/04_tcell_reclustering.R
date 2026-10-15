# =====================================================================
# Step 04 -- T cells: re-integrate and recluster from scratch
# Input  : objects/03_tcell_subset.rds
# Output : objects/04_tcell_reclustered.rds   (tcell_res.<r>; tcell_clusters = res 0.4)
#          results/04_tcell_reclustering/tcell_subset_markers.csv (+ top10, panel averages)
#          figures/04_tcell_reclustering/*
#
# Res-0.4 annotation (11 clusters; final labels are applied in step 06):
#   0 CD8 TEM GZMK/GZMH | 1 CD4 naive/CM | 2 CD8 TEM ZEB2-hi | 3 Th17 | 4 Treg
#   5 CD4 resting/CM TSHZ2+ | 6 CD8 cytotoxic GZMB/GNLY | 7 MAIT
#   8 CD4 Tph-like CXCL13+ | 9 ribosomal-gene-dominated, no lineage markers
#   10 activated/transitional, inside cluster 0 territory
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(harmony)
})
source("config.R"); source("R/helpers.R")
STEP <- "04_tcell_reclustering"
d <- step_dirs(STEP)
dir.create(file.path(d$figures, "resolution_sweep"), showWarnings = FALSE)
set.seed(SEED_GLOBAL)

tcell_obj <- readRDS(OBJ$tcell)
print(table(tcell_obj$sample))

# ---------------------------------------------------------------------
# 1. Variable features and PCA within T cells
# ---------------------------------------------------------------------
DefaultAssay(tcell_obj) <- "RNA"
tcell_obj[["RNA"]] <- split(tcell_obj[["RNA"]], f = tcell_obj$sample)
tcell_obj <- NormalizeData(tcell_obj)
tcell_obj <- FindVariableFeatures(tcell_obj, nfeatures = N_VARIABLE_FEATURES)
tcell_obj <- ScaleData(tcell_obj)
tcell_obj <- RunPCA(tcell_obj, npcs = TCELL_NPCS, seed.use = SEED_PCA_UMAP)
plot_elbow(tcell_obj, "T-cell subset -- elbow plot (PCA)", file.path(d$figures, "elbow_tcell_pca.svg"))

tcell_obj <- RunUMAP(tcell_obj, reduction = "pca", dims = TCELL_DIMS,
                     reduction.name = "umap.unintegrated.tcell", seed.use = SEED_PCA_UMAP)

# ---------------------------------------------------------------------
# 2. Harmony integration (batch = sample) and UMAP
# ---------------------------------------------------------------------
set.seed(SEED_HARMONY_T)
tcell_obj <- IntegrateLayers(tcell_obj, method = HarmonyIntegration, orig.reduction = "pca",
                             new.reduction = "integrated.harmony.tcell")
tcell_obj <- JoinLayers(tcell_obj)
tcell_obj <- RunUMAP(tcell_obj, reduction = "integrated.harmony.tcell", dims = TCELL_DIMS,
                     reduction.name = TCELL_UMAP, seed.use = SEED_PCA_UMAP)

for (red in c("umap.unintegrated.tcell", TCELL_UMAP)) {
  tag <- if (red == TCELL_UMAP) "harmony" else "unintegrated"
  plot_umap(tcell_obj, red, "sample", paste("T cells --", tag, "by sample"),
            file.path(d$figures, paste0("umap_tcell_", tag, "_by_sample")),
            facet = TRUE, facet_ncol = 5, width = 16, height = 16, point_size = 0.5)
  plot_umap(tcell_obj, red, "condition", paste("T cells --", tag, "by condition"),
            file.path(d$figures, paste0("umap_tcell_", tag, "_by_condition")),
            facet = TRUE, facet_ncol = 2, width = 9, height = 5, point_size = 0.5)
}

# ---------------------------------------------------------------------
# 3. Resolution sweep; res 0.4 carried forward as tcell_clusters
# ---------------------------------------------------------------------
tcell_obj <- FindNeighbors(tcell_obj, reduction = "integrated.harmony.tcell", dims = TCELL_DIMS)
tcell_obj <- FindClusters(tcell_obj, resolution = TCELL_RES_SWEEP, random.seed = SEED_CLUSTERS)
tcell_obj <- rename_res_cols(tcell_obj, "tcell_res.")
for (cc in res_cols(tcell_obj, "tcell_res.")) {
  r <- sub("tcell_res.", "", cc, fixed = TRUE)
  plot_umap(tcell_obj, TCELL_UMAP, cc,
            sprintf("T cells -- res %s (%d clusters)", r, length(unique(tcell_obj[[cc]][, 1]))),
            file.path(d$figures, "resolution_sweep", paste0("umap_tcell_res", r)), point_size = 0.5)
}

tcell_obj$tcell_clusters <- tcell_obj[[TCELL_RES_COL]][, 1]
Idents(tcell_obj) <- "tcell_clusters"

# ---------------------------------------------------------------------
# 4. Markers and T-cell subset panels
# ---------------------------------------------------------------------
markers_and_top10(tcell_obj, file.path(d$results, "tcell_subset_markers.csv"),
                  file.path(d$results, "tcell_subset_top10_markers.csv"))

tcell_panel <- list(
  CD4           = c("CD4", "IL7R"),
  CD8           = c("CD8A", "CD8B"),
  Naive         = c("CCR7", "SELL", "TCF7", "LEF1"),
  Treg          = c("FOXP3", "IL2RA", "CTLA4"),
  Cytotoxic     = c("GZMA", "GZMB", "GZMK", "PRF1", "NKG7"),
  GammaDelta    = c("TRDC", "TRGC1", "TRGC2"),
  Proliferating = c("MKI67", "TOP2A")
)
avail <- lapply(tcell_panel, function(g) intersect(g, rownames(tcell_obj)))
for (grp in names(avail)) {
  if (length(avail[[grp]]) == 0) next
  plot_marker_panel(tcell_obj, TCELL_UMAP, avail[[grp]], paste(grp, "markers"),
                    file.path(d$figures, paste0("tcell_subset_markers_", grp)))
}
avg_expr <- AverageExpression(tcell_obj, features = unlist(avail, use.names = FALSE),
                              group.by = "tcell_clusters", layer = "data")$RNA
write.csv(round(avg_expr, 2), file.path(d$results, "tcell_subset_marker_avg_expr.csv"))

saveRDS(tcell_obj, OBJ$tcell_cluster)
save_session(STEP, d)
log_step("Done: ", OBJ$tcell_cluster)
