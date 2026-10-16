# =====================================================================
# Step 00 -- QC, merge and integration of the 21 GSE164241 samples
#
# Author : Tyler Chinsky (original analysis, 2026-08-28)
# Input  : Cell Ranger outputs, one folder per sample, in RAW_10X_DIR
#          (<sample>/outs/filtered_feature_bc_matrix)
# Output : objects/00_integrated.rds        (the pipeline's starting object)
#          results/00_integration/*         QC cell counts, cluster markers
#          figures/00_integration/*         QC and UMAP plots
#
# Status : NOT re-verified. The published integrated object (md5
#          a3926732393aa6a5e7387cf13565016e) is the verified starting point
#          of this pipeline; steps 01-09 reproduce every published result
#          from it. This script mirrors the code that produced that object.
#          Automated cell-type annotation run in that session (10x Cloud,
#          scType, CloudAzimuth) was not used in the analysis and is omitted;
#          cell types were annotated manually from marker genes (steps 03-06).
#
# Repeat biopsies (PD161/PD161b, PD164/PD164b/PD164c) are kept
# as independent samples: each is its own integration batch. `donor`
# records the base subject for tracking only.
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(harmony)
})
source("config.R"); source("R/helpers.R")
STEP <- "00_integration"
d <- step_dirs(STEP)
set.seed(SEED_GLOBAL)

# ---------------------------------------------------------------------
# 1. Locate samples
# ---------------------------------------------------------------------
sample_ids  <- list.files(RAW_10X_DIR, pattern = "^(GM|PD)")
sample_ids  <- sample_ids[dir.exists(file.path(RAW_10X_DIR, sample_ids, "outs"))]
sample_dirs <- setNames(file.path(RAW_10X_DIR, sample_ids, "outs", "filtered_feature_bc_matrix"), sample_ids)
stopifnot(length(sample_dirs) > 0)
missing <- names(sample_dirs)[!dir.exists(sample_dirs)]
if (length(missing)) stop("Missing matrix dirs: ", paste(missing, collapse = ", "))
log_step("Found ", length(sample_dirs), " samples")

# ---------------------------------------------------------------------
# 2. Load and QC each sample
# ---------------------------------------------------------------------

load_sample <- function(sample_name) {
  log_step("Loading ", sample_name)
  obj <- CreateSeuratObject(counts = Read10X(sample_dirs[[sample_name]]), project = sample_name,
                            min.cells = QC_MIN_CELLS_PER_GENE, min.features = QC_MIN_FEATURES)
  obj$sample    <- sample_name
  obj$condition <- ifelse(grepl("^PD", sample_name), PD_LABEL, GM_LABEL)
  obj$donor     <- sub("[a-c]$", "", sample_name)
  obj$replicate <- ifelse(grepl("[a-c]$", sample_name), sub("^.*([a-c])$", "\\1", sample_name), "a")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  ggsave(file.path(d$figures, paste0(sample_name, "_qc_metrics.png")), width = 12, height = 4, dpi = 150)
  ggsave(file.path(d$figures, paste0(sample_name, "_qc_scatter.png")),
         plot = FeatureScatter(obj, "nCount_RNA", "percent.mt") + FeatureScatter(obj, "nCount_RNA", "nFeature_RNA"),
         width = 12, height = 6, dpi = 150)

  n_before <- ncol(obj)
  obj <- subset(obj, subset = nFeature_RNA > QC_MIN_FEATURES & nFeature_RNA < QC_MAX_FEATURES &
                  percent.mt < QC_MAX_PERCENT_MT)
  log_step("  ", sample_name, ": ", n_before, " -> ", ncol(obj), " cells")
  obj
}

seurat_list <- setNames(lapply(names(sample_dirs), load_sample), names(sample_dirs))
seurat_list <- seurat_list[vapply(seurat_list, ncol, double(1)) > QC_MIN_CELLS_SAMPLE]

qc_summary <- data.frame(sample = names(seurat_list), cells = vapply(seurat_list, ncol, double(1)),
                         condition = ifelse(grepl("^PD", names(seurat_list)), "PD", "GM"),
                         donor = sub("[a-c]$", "", names(seurat_list)))
write.csv(qc_summary, file.path(d$results, "qc_cell_counts.csv"), row.names = FALSE)
reps <- qc_summary %>% dplyr::count(donor) %>% dplyr::filter(n > 1)

# ---------------------------------------------------------------------
# 3. Merge; one counts layer per sample = the integration batches
# ---------------------------------------------------------------------
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list),
                  project = "GSE164241")
rm(seurat_list); gc()
if (length(Layers(combined, search = "counts")) == 1) {
  combined[["RNA"]] <- split(combined[["RNA"]], f = combined$sample)
}

# ---------------------------------------------------------------------
# 4. Normalise, variable features, PCA, unintegrated UMAP
# ---------------------------------------------------------------------
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = N_VARIABLE_FEATURES)
ggsave(file.path(d$figures, "combined_variable_features.png"),
       LabelPoints(VariableFeaturePlot(combined), points = head(VariableFeatures(combined), 10), repel = TRUE),
       width = 8, height = 6, dpi = 150)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = INT_NPCS)
ggsave(file.path(d$figures, "combined_elbow_plot.png"), ElbowPlot(combined, ndims = INT_NPCS),
       width = 8, height = 6, dpi = 150)
ggsave(file.path(d$figures, "combined_pca_by_sample.png"),
       DimPlot(combined, reduction = "pca", group.by = "sample"), width = 9, height = 6, dpi = 150)
png(file.path(d$figures, "combined_pca_heatmap.png"), width = 1200, height = 1600)
DimHeatmap(combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
combined <- RunUMAP(combined, dims = INT_UMAP_DIMS_UNINT, reduction = "pca",
                    reduction.name = "umap.unintegrated")

# ---------------------------------------------------------------------
# 5. Integration: CCA and Harmony on the same split object
# ---------------------------------------------------------------------
combined <- IntegrateLayers(combined, method = CCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.cca", verbose = TRUE)
combined <- IntegrateLayers(combined, method = HarmonyIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.harmony", verbose = TRUE)
combined <- JoinLayers(combined)

# ---------------------------------------------------------------------
# 6. Cluster + UMAP for each integration
# ---------------------------------------------------------------------
methods <- c(cca = "integrated.cca", harmony = "integrated.harmony")
for (m in names(methods)) {
  combined <- FindNeighbors(combined, dims = INT_CLUSTER_DIMS, reduction = methods[[m]])
  combined <- FindClusters(combined, resolution = INT_CLUSTER_RES, cluster.name = paste0("clusters_", m))
  combined <- RunUMAP(combined, dims = INT_CLUSTER_DIMS, reduction = methods[[m]],
                      reduction.name = paste0("umap.", m))
  u <- paste0("umap.", m)
  for (g in c("sample", "condition")) {
    ggsave(file.path(d$figures, sprintf("umap_%s_by_%s.png", m, g)),
           DimPlot(combined, reduction = u, group.by = g), width = 10, height = 7, dpi = 150)
  }
  ggsave(file.path(d$figures, sprintf("umap_%s_clusters.png", m)),
         DimPlot(combined, reduction = u, group.by = paste0("clusters_", m), label = TRUE) + NoLegend(),
         width = 9, height = 7, dpi = 150)
}
write.csv(as.data.frame.matrix(table(combined$clusters_cca, combined$clusters_harmony)),
          file.path(d$results, "cluster_concordance_cca_vs_harmony.csv"))

# Harmony is the primary integration used by every downstream step
PRIMARY <- "harmony"
Idents(combined) <- paste0("clusters_", PRIMARY)
combined$seurat_clusters <- combined[[paste0("clusters_", PRIMARY)]][, 1]
UMAP_PRIMARY <- paste0("umap.", PRIMARY)

# ---------------------------------------------------------------------
# 7. Markers
# ---------------------------------------------------------------------
markers_and_top10(combined, file.path(d$results, "cluster_markers.csv"),
                  file.path(d$results, "top10_markers_per_cluster.csv"))

saveRDS(combined, OBJ$integrated)
save_session(STEP, d)
log_step("Done: ", OBJ$integrated)
