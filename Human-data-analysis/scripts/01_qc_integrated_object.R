# =====================================================================
# Step 01 -- Sanity checks on the integrated object
# Input  : objects/00_integrated.rds
# Output : results/01_qc/qc_report.txt, qc_per_sample.csv, qc_checks.csv
# Read-only on the object. Fails (non-zero exit) if a check fails.
# =====================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("config.R"); source("R/helpers.R")
STEP <- "01_qc"
d <- step_dirs(STEP)

combined <- readRDS(OBJ$integrated)
md <- combined@meta.data
sink(file.path(d$results, "qc_report.txt"), split = TRUE)

cat("== Object ==\n"); print(combined)
cat("Reductions:", paste(Reductions(combined), collapse = ", "), "\n")
cat("Counts layers:", length(Layers(combined, search = "counts")), "\n")

cat("\n== Cells per sample / condition ==\n")
print(table(md$sample)); print(table(md$condition))
donor_tab <- md %>% distinct(sample, donor) %>% count(donor)
reps <- donor_tab %>% filter(n > 1)
cat("Donors with repeat biopsies:", paste(reps$donor, collapse = ", "), "\n")

qc_per_sample <- md %>% group_by(sample) %>%
  summarise(n_cells = n(), median_nFeature = median(nFeature_RNA), median_nCount = median(nCount_RNA),
            median_pctMT = median(percent.mt), max_pctMT = max(percent.mt),
            min_nFeature = min(nFeature_RNA), max_nFeature = max(nFeature_RNA))
write.csv(qc_per_sample, file.path(d$results, "qc_per_sample.csv"), row.names = FALSE)
print(as.data.frame(qc_per_sample))

cat("\n== Clustering ==\n")
for (cl in c("clusters_cca", "clusters_harmony", "seurat_clusters")) {
  if (cl %in% colnames(md)) cat(cl, ":", length(unique(md[[cl]])), "clusters\n")
}
ct <- table(md$seurat_clusters, md$sample)
single_sample <- rownames(ct)[rowSums(ct > 0) == 1]

cat("\n== Replicate concordance (cells per cluster, repeat-biopsy donors) ==\n")
for (dn in reps$donor) print(table(md$sample[md$donor == dn], md$seurat_clusters[md$donor == dn]))

checks <- c(
  "counts layer joined"           = length(Layers(combined, search = "counts")) == 1,
  "no duplicated barcodes"        = !any(duplicated(colnames(combined))),
  "nFeature within QC bounds"     = all(md$nFeature_RNA > QC_MIN_FEATURES & md$nFeature_RNA < QC_MAX_FEATURES),
  "percent.mt below QC bound"     = all(md$percent.mt < QC_MAX_PERCENT_MT),
  "one condition per sample"      = !any(duplicated(distinct(md, sample, condition)$sample)),
  "integrated.harmony present"    = GLOBAL_REDUCTION %in% Reductions(combined),
  "umap.harmony present"          = GLOBAL_UMAP %in% Reductions(combined),
  "repeat-biopsy donors detected" = nrow(reps) > 0
)
cat("\n== Checks ==\n"); print(checks)
if (length(single_sample)) {
  cat("WARNING: clusters present in only one sample:", paste(single_sample, collapse = ", "), "\n")
}
sink()
write.csv(data.frame(check = names(checks), pass = unname(checks)),
          file.path(d$results, "qc_checks.csv"), row.names = FALSE)
save_session(STEP, d)
if (!all(checks)) stop("QC checks failed: ", paste(names(checks)[!checks], collapse = "; "))
log_step("All QC checks passed")
