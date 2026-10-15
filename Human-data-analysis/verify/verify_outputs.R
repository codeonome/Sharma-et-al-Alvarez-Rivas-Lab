# =====================================================================
# Verify pipeline outputs against the published results.
#   Rscript verify/verify_outputs.R <published_root> [step ...]
#   (run from pipeline/; PIPELINE_OUT selects the output root as usual)
# <published_root> is the folder holding 01_Global_Analysis/,
# 02_Immune_Subset/, 03_Tcell_Reclustering/, 04_Tcell_Final/.
# Tables: IDENTICAL (bytes) or EQUAL (all values equal within 1e-9).
# Objects: same cells, ARI = 1 for every cluster column, embeddings equal.
# Writes verify/report_<steps>.csv next to the outputs; exits non-zero on
# any mismatch.
# =====================================================================
suppressPackageStartupMessages(library(Seurat))
source("config.R")
a <- commandArgs(trailingOnly = TRUE)
PUB <- a[1]
STEPS <- if (length(a) > 1) a[-1] else sprintf("%02d", 2:9)

M01 <- "01_Global_Analysis"; M02 <- "02_Immune_Subset"; M03 <- "03_Tcell_Reclustering"; M04 <- "04_Tcell_Final"
tables <- rbind(
  data.frame(step = "02", new = "02_global_clustering/cluster_quality_flags.csv", pub = file.path(M01, "results/cluster_quality_flags.csv")),
  data.frame(step = "03", new = file.path("03_immune_subset", c("cluster_markers_res0.1.csv", "top10_markers_res0.1.csv")),
             pub = file.path(M01, "results", c("cluster_markers_res0.1.csv", "top10_markers_res0.1.csv"))),
  data.frame(step = "03", new = file.path("03_immune_subset", c("immune_subset_markers_res0.2.csv", "immune_subset_top10_markers_res0.2.csv", "immune_subset_marker_avg_expr.csv")),
             pub = file.path(M02, "results", c("immune_subset_markers_res0.2.csv", "immune_subset_top10_markers_res0.2.csv", "immune_subset_marker_avg_expr.csv"))),
  data.frame(step = "04", new = file.path("04_tcell_reclustering", c("tcell_subset_markers.csv", "tcell_subset_top10_markers.csv", "tcell_subset_marker_avg_expr.csv")),
             pub = file.path(M03, "results", c("tcell_subset_markers.csv", "tcell_subset_top10_markers.csv", "tcell_subset_marker_avg_expr.csv"))),
  data.frame(step = "05", new = file.path("05_tcell_validation", paste0("validation_cluster", c("3_frac_expressing", "8_frac_expressing", "9_qc_metrics"), ".csv")),
             pub = file.path(M03, "results", paste0("validation_cluster", c("3_frac_expressing", "8_frac_expressing", "9_qc_metrics"), ".csv"))),
  data.frame(step = "06", new = file.path("06_tcell_cleanup_annotation", f06 <- c("sample_to_donor_map.csv", "ig_doublets_by_cluster.csv", "nk_check_by_cluster.csv",
                                                                                "tcell_cleanup_changelog.csv", "cluster_renaming_map.csv")),
             pub = file.path(M04, "results", f06)),
  data.frame(step = "07", new = file.path("07_differential_abundance", f07 <- c(
    "sample_x_cluster_long.csv", "sample_x_cluster_counts.csv", "sample_x_cluster_proportions.csv",
    "donor_x_cluster_long.csv", "donor_x_cluster_proportions.csv", "DA_wilcox_sample_level.csv",
    "DA_wilcox_donor_level.csv", "DA_propeller_sample_level.csv", "DA_propeller_donor_level.csv",
    "DA_propeller_donor_blocked.csv", "intra_donor_correlation.txt", "DA_all_tests_combined.csv",
    "DA_all_tests_wide.csv", "DA_leave_one_donor_out.csv", "DA_leave_one_donor_out_summary.csv",
    "intra_donor_variance_components.csv", "DA_propeller_rho_sensitivity.csv",
    "DA_propeller_rho_sensitivity_wide.csv", "DA_propeller_rho_sensitivity_summary.csv")),
    pub = file.path(M04, "results", f07)),
  data.frame(step = "08", new = file.path("08_figures", f08 <- c("donut_pooled_percentages.csv", "DA_bar_blocked_values.csv")),
             pub = file.path(M04, "results", f08)),
  data.frame(step = "09", new = file.path("09_pseudobulk_de", f09 <- c(
    paste0("DE_PD_vs_GM_", c("allTcells", paste0("C", 1:9)), ".csv"),
    "DE_per_cluster_all_genes.csv", "DE_per_cluster_summary.csv", "DE_highlight_genes_per_cluster.csv")),
    pub = file.path(M04, "results", f09))
)

objects <- list(
  "02" = list(new = OBJ$global, pub = file.path(M01, "objects/integrated_seurat_object_res_sweep.rds"),
              cols = setNames(paste0("RNA_snn_res.", GLOBAL_RES_SWEEP), paste0("global_res.", GLOBAL_RES_SWEEP)),
              reds = character(0)),
  "03" = list(new = OBJ$immune, pub = file.path(M02, "objects/immune_subset_seurat_object.rds"),
              cols = c(setNames(paste0("immune_res.", IMMUNE_RES_SWEEP), paste0("immune_res.", IMMUNE_RES_SWEEP)),
                       immune_clusters = "immune_clusters"),
              reds = c("pca", "integrated.harmony.immune", "umap.immune")),
  "03t" = list(new = OBJ$tcell, pub = file.path(M02, "objects/tcell_seurat_object.rds"),
               cols = c(immune_clusters = "immune_clusters"), reds = character(0)),
  "04" = list(new = OBJ$tcell_cluster, pub = file.path(M03, "objects/tcell_seurat_object_harmony_after_clustering_res0.4.rds"),
              cols = c(setNames(paste0("tcell_res.", TCELL_RES_SWEEP), paste0("tcell_res.", TCELL_RES_SWEEP)),
                       tcell_clusters = "tcell_clusters"),
              reds = c("pca", "integrated.harmony.tcell", TCELL_UMAP)),
  "06" = list(new = OBJ$tcell_final, pub = file.path(M04, "objects/tcell_final_seurat_object.rds"),
              cols = c(cluster_id = "cluster_id", cluster_code = "cluster_code", donor = "donor"),
              reds = TCELL_UMAP)
)

cmp_table <- function(new, pub) {
  if (!file.exists(new)) return(c("MISSING", "output not found"))
  if (unname(tools::md5sum(new)) == unname(tools::md5sum(pub))) return(c("IDENTICAL", ""))
  if (grepl("\\.txt$", new)) return(c("DIFFERENT", "text differs"))
  x <- read.csv(new, check.names = FALSE, stringsAsFactors = FALSE)
  y <- read.csv(pub, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(dim(x), dim(y))) return(c("DIFFERENT", sprintf("dim %s vs %s", paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x"))))
  if (!identical(names(x), names(y))) return(c("DIFFERENT", "column names differ"))
  worst <- 0
  for (j in seq_along(x)) {
    if (is.numeric(x[[j]]) && is.numeric(y[[j]])) {
      dd <- abs(x[[j]] - y[[j]]) / pmax(1, abs(y[[j]]))
      if (any(is.na(x[[j]]) != is.na(y[[j]]))) return(c("DIFFERENT", paste("NA pattern:", names(x)[j])))
      worst <- max(worst, dd, na.rm = TRUE)
    } else if (!identical(as.character(x[[j]]), as.character(y[[j]]))) {
      return(c("DIFFERENT", paste("values differ in column", names(x)[j])))
    }
  }
  if (worst < 1e-9) c("EQUAL", sprintf("max rel. diff %.1e", worst)) else c("DIFFERENT", sprintf("max rel. diff %.1e", worst))
}

cmp_object <- function(spec) {
  if (!file.exists(spec$new)) return(data.frame(item = "object", status = "MISSING", detail = spec$new))
  x <- readRDS(spec$new); y <- readRDS(file.path(PUB, spec$pub))
  out <- list(data.frame(item = "cells", status = if (setequal(colnames(x), colnames(y))) "IDENTICAL" else "DIFFERENT",
                         detail = sprintf("%d vs %d", ncol(x), ncol(y))))
  cells <- intersect(colnames(x), colnames(y))
  for (nc in names(spec$cols)) {
    a1 <- as.character(x@meta.data[cells, nc]); b1 <- as.character(y@meta.data[cells, spec$cols[[nc]]])
    ari <- mclust::adjustedRandIndex(a1, b1)
    st <- if (identical(a1, b1)) "IDENTICAL" else if (ari == 1) "EQUAL" else "DIFFERENT"
    out[[length(out) + 1]] <- data.frame(item = nc, status = st, detail = sprintf("ARI %.6f", ari))
  }
  for (r in spec$reds) {
    e1 <- Embeddings(x, r)[cells, ]; e2 <- Embeddings(y, r)[cells, ]
    md <- max(abs(e1 - e2))
    out[[length(out) + 1]] <- data.frame(item = r, status = if (md == 0) "IDENTICAL" else if (md < 1e-6) "EQUAL" else "DIFFERENT",
                                         detail = sprintf("max |diff| %.2e", md))
  }
  do.call(rbind, out)
}

rows <- list()
for (i in which(tables$step %in% STEPS)) {
  r <- cmp_table(file.path(OUT_ROOT, "results", tables$new[i]), file.path(PUB, tables$pub[i]))
  rows[[length(rows) + 1]] <- data.frame(step = tables$step[i], type = "table", item = basename(tables$new[i]),
                                         status = r[1], detail = r[2])
}
for (k in names(objects)) {
  if (!substr(k, 1, 2) %in% STEPS) next
  o <- cmp_object(objects[[k]])
  rows[[length(rows) + 1]] <- data.frame(step = k, type = "object", item = paste(basename(objects[[k]]$new), o$item, sep = ": "),
                                         status = o$status, detail = o$detail)
}
rep <- do.call(rbind, rows)
options(width = 200); print(rep, row.names = FALSE)
dir.create(file.path(OUT_ROOT, "verify"), showWarnings = FALSE)
write.csv(rep, file.path(OUT_ROOT, "verify", paste0("report_steps_", paste(STEPS, collapse = "-"), ".csv")), row.names = FALSE)
cat(sprintf("\n%d checks: %d IDENTICAL, %d EQUAL, %d DIFFERENT, %d MISSING\n", nrow(rep),
            sum(rep$status == "IDENTICAL"), sum(rep$status == "EQUAL"),
            sum(rep$status == "DIFFERENT"), sum(rep$status == "MISSING")))
if (any(rep$status %in% c("DIFFERENT", "MISSING"))) quit(status = 1)
