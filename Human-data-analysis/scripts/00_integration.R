# =====================================================================
# Step 00 -- QC, merge and integration of the 21 GSE164241 samples
#
# Author : Tyler Chinsky (original analysis, CCR, 2026-08-28)
# Input  : Cell Ranger outputs, one folder per sample, in RAW_10X_DIR
#          (<sample>/outs/filtered_feature_bc_matrix and, if present,
#           <sample>/outs/cell_types/cell_types.csv from 10x Cloud)
# Output : objects/00_integrated.rds        (the pipeline's starting object)
#          results/00_integration/*         QC counts, markers, annotations
#          figures/00_integration/*         QC and UMAP plots
#
# Status : NOT re-verified. The published integrated object (md5
#          a3926732393aa6a5e7387cf13565016e) is the verified starting point
#          of this pipeline; steps 01-09 reproduce every published result
#          from it. This script mirrors the code that produced that object
#          (/projects/academic/carlaalv/Human_integration/seurat_integration.R).
#          SingleR was not run in that session and is omitted here.
#          Annotation (scType, CloudAzimuth) needs internet access.
#
# Repeat biopsies (GM184a, PD161/PD161b, PD164/PD164b/PD164c) are kept
# as independent samples: each is its own integration batch. `donor`
# records the base subject for tracking only.
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(stringr); library(harmony)
  library(HGNChelper)   # scType
  library(AzimuthAPI)   # CloudAzimuth
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
annot_files <- setNames(file.path(RAW_10X_DIR, sample_ids, "outs", "cell_types", "cell_types.csv"), sample_ids)
stopifnot(length(sample_dirs) > 0)
missing <- names(sample_dirs)[!dir.exists(sample_dirs)]
if (length(missing)) stop("Missing matrix dirs: ", paste(missing, collapse = ", "))
log_step("Found ", length(sample_dirs), " samples; ", sum(file.exists(annot_files)), " with 10x cell_types.csv")

# ---------------------------------------------------------------------
# 2. Load, annotate and QC each sample
# ---------------------------------------------------------------------
ANNOT_10X <- c("coarse_cell_type", "fine_cell_type", "cell_count_in_model",
               "coarse_cell_type_unfiltered", "umi_count")

load_sample <- function(sample_name) {
  log_step("Loading ", sample_name)
  obj <- CreateSeuratObject(counts = Read10X(sample_dirs[[sample_name]]), project = sample_name,
                            min.cells = QC_MIN_CELLS_PER_GENE, min.features = QC_MIN_FEATURES)
  obj$sample    <- sample_name
  obj$condition <- ifelse(grepl("^PD", sample_name), PD_LABEL, GM_LABEL)
  obj$donor     <- sub("[a-c]$", "", sample_name)
  obj$replicate <- ifelse(grepl("[a-c]$", sample_name), sub("^.*([a-c])$", "\\1", sample_name), "a")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # 10x Cloud labels must be attached before merge(), while barcodes are bare
  af <- annot_files[[sample_name]]
  if (file.exists(af)) {
    ct <- read.csv(af, stringsAsFactors = FALSE)
    bc_col <- grep("^barcode", names(ct), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(bc_col)) stop("No barcode column in ", af)
    rownames(ct) <- ct[[bc_col]]
    want <- intersect(ANNOT_10X, names(ct))
    ann <- ct[colnames(obj), want, drop = FALSE]
    rownames(ann) <- colnames(obj)
    if (sum(!is.na(ann[[want[1]]])) == 0) warning(sample_name, ": zero barcode overlap with cell_types.csv")
    obj <- AddMetaData(obj, ann)
  } else {
    for (cc in ANNOT_10X) obj[[cc]] <- NA
  }

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

# ---------------------------------------------------------------------
# 8. Annotation: 10x Cloud labels, scType, CloudAzimuth
# ---------------------------------------------------------------------
if (any(!is.na(combined$coarse_cell_type))) {
  for (lv in c("coarse", "fine")) {
    col <- paste0(lv, "_cell_type")
    write.csv(as.data.frame.matrix(table(Idents(combined), combined[[col]][, 1], useNA = "ifany")),
              file.path(d$results, sprintf("cluster_vs_10x_%s.csv", lv)))
  }
}

# scType: one label per cluster from marker enrichment on scaled data
sctype_ok <- tryCatch({
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  TRUE
}, error = function(e) { message("scType source failed: ", conditionMessage(e)); FALSE })
if (sctype_ok) {
  gs_list <- gene_sets_prepare(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", SCTYPE_TISSUE)
  es.max <- sctype_score(scRNAseqData = GetAssayData(combined, assay = "RNA", layer = "scale.data"),
                         scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  cl_res <- do.call("rbind", lapply(unique(Idents(combined)), function(cl) {
    cells <- WhichCells(combined, idents = cl)
    es.cl <- sort(rowSums(es.max[, colnames(es.max) %in% cells, drop = FALSE]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.cl), scores = es.cl, ncells = length(cells)), 10)
  }))
  top_sctype <- cl_res %>% group_by(cluster) %>% slice_max(scores, n = 1)
  top_sctype$type[top_sctype$scores < top_sctype$ncells / 4] <- "Unknown"
  write.csv(cl_res, file.path(d$results, "sctype_all_scores.csv"), row.names = FALSE)
  write.csv(top_sctype, file.path(d$results, "sctype_cluster_calls.csv"), row.names = FALSE)
  combined$scType <- as.character(top_sctype$type[match(as.character(Idents(combined)),
                                                        as.character(top_sctype$cluster))])
} else {
  combined$scType <- NA_character_
}

# CloudAzimuth: hosted hierarchical reference mapping (broad/medium/fine)
AZ_COLS <- c("azimuth_broad", "azimuth_medium", "azimuth_fine", "azimuth_label",
             "full_hierarchical_labels", "final_level_labels", "final_level_confidence",
             "full_consistent_hierarchy")
tryCatch({
  az <- CloudAzimuth(combined)
  for (cc in intersect(AZ_COLS, colnames(az[[]]))) combined[[cc]] <- az[[cc]][, 1]
  rm(az); gc()
  if ("final_level_confidence" %in% colnames(combined@meta.data)) {
    combined$final_level_confidence <- as.numeric(combined$final_level_confidence)
  }
  write.csv(as.data.frame.matrix(table(Idents(combined), combined$azimuth_fine, useNA = "ifany")),
            file.path(d$results, "cluster_vs_azimuth_fine.csv"))
}, error = function(e) message("CloudAzimuth failed: ", conditionMessage(e)))

# Per-cluster majority label from every available annotation, side by side
annot_cols <- intersect(c("coarse_cell_type", "fine_cell_type", "scType",
                          "azimuth_broad", "azimuth_medium", "azimuth_fine"), colnames(combined@meta.data))
annot_cols <- annot_cols[vapply(annot_cols, function(x) any(!is.na(combined[[x]][, 1])), logical(1))]
consensus <- lapply(annot_cols, function(ac) {
  combined@meta.data %>%
    mutate(cluster = as.character(Idents(combined))) %>%
    filter(!is.na(.data[[ac]])) %>%
    dplyr::count(cluster, label = .data[[ac]]) %>%
    group_by(cluster) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(cluster, method = ac, label)
}) %>% bind_rows() %>%
  pivot_wider(names_from = method, values_from = label) %>%
  arrange(suppressWarnings(as.numeric(cluster)))
write.csv(consensus, file.path(d$results, "cluster_annotation_consensus.csv"), row.names = FALSE)

# Composition per sample for each label source (repeat biopsies as separate rows)
for (lc in annot_cols) {
  comp <- combined@meta.data %>%
    filter(!is.na(.data[[lc]])) %>%
    dplyr::count(sample, donor, condition, label = .data[[lc]]) %>%
    group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
  write.csv(comp, file.path(d$results, sprintf("composition_%s.csv", lc)), row.names = FALSE)
}

saveRDS(combined, OBJ$integrated)
save_session(STEP, d)
log_step("Done: ", OBJ$integrated)
