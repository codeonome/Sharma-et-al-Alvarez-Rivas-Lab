# =====================================================================
# Step 05 -- Targeted validation of T-cell clusters 3, 8 and 9 (res 0.4)
#   cluster 3 -> Treg vs Th17 vs Th17-like Treg   (FOXP3 / IKZF2 / RORC)
#   cluster 8 -> Tph vs Tfh                        (CXCR5 / PDCD1 / BCL6)
#   cluster 9 -> population vs technical state     (QC metrics, ribosomal %)
# Input  : objects/04_tcell_reclustered.rds
# Output : results/05_tcell_validation/validation_cluster{3,8,9}_*.csv
#          figures/05_tcell_validation/*
# Read-only on the object.
# =====================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr); library(ggplot2) })
source("config.R"); source("R/helpers.R")
STEP <- "05_tcell_validation"
d <- step_dirs(STEP)

tcell_obj <- readRDS(OBJ$tcell_cluster)
stopifnot("tcell_clusters" %in% colnames(tcell_obj@meta.data))
Idents(tcell_obj) <- "tcell_clusters"
sink(file.path(d$results, "validation_report.txt"), split = TRUE)

# ---------------------------------------------------------------------
# Cluster 3: Treg / Th17 modules
# ---------------------------------------------------------------------
cluster3_panel <- c("FOXP3", "IKZF2", "RORC", "IL2RA", "CTLA4", "TNFRSF18", "TNFRSF4",
                    "CCR6", "IL1R1", "IL17A", "IL23R", "CCL20")
plot_violin_panel(tcell_obj, cluster3_panel, "tcell_clusters",
                  "Cluster 3 validation -- Treg vs Th17 markers",
                  file.path(d$figures, "cluster3_treg_th17_panel"), highlight = "3", ncol = 3)
write.csv(frac_expressing(tcell_obj, cluster3_panel, "tcell_clusters"),
          file.path(d$results, "validation_cluster3_frac_expressing.csv"), row.names = FALSE)

c3 <- FetchData(tcell_obj, vars = c("FOXP3", "RORC", "IKZF2"), cells = WhichCells(tcell_obj, idents = "3"))
cat("Cluster 3 co-expression (FOXP3 x RORC):\n"); print(table(FOXP3 = c3$FOXP3 > 0, RORC = c3$RORC > 0))
cat(sprintf("FOXP3+ cells also RORC+: %.1f%% | IKZF2+: %.1f%%\n",
            100 * mean(c3$RORC[c3$FOXP3 > 0] > 0), 100 * mean(c3$IKZF2[c3$FOXP3 > 0] > 0)))

# ---------------------------------------------------------------------
# Cluster 8: Tph (CXCR5-, PDCD1+, CXCL13+) vs Tfh (CXCR5+, BCL6+)
# ---------------------------------------------------------------------
cluster8_panel <- c("CXCR5", "PDCD1", "BCL6", "CXCL13", "TOX2", "PTPN13", "ICOS", "IL21", "MAF", "CD200")
plot_violin_panel(tcell_obj, cluster8_panel, "tcell_clusters",
                  "Cluster 8 validation -- Tph vs Tfh markers",
                  file.path(d$figures, "cluster8_tph_tfh_panel"), highlight = "8", ncol = 3)
write.csv(frac_expressing(tcell_obj, cluster8_panel, "tcell_clusters"),
          file.path(d$results, "validation_cluster8_frac_expressing.csv"), row.names = FALSE)

c8 <- FetchData(tcell_obj, vars = intersect(c("CXCR5", "PDCD1", "BCL6", "CXCL13"), rownames(tcell_obj)),
                cells = WhichCells(tcell_obj, idents = "8"))
cat("\nCluster 8 co-expression (CXCR5 x PDCD1):\n"); print(table(CXCR5 = c8$CXCR5 > 0, PDCD1 = c8$PDCD1 > 0))

# ---------------------------------------------------------------------
# Cluster 9: QC metrics vs all other clusters
# ---------------------------------------------------------------------
if (!"percent.ribo" %in% colnames(tcell_obj@meta.data)) {
  tcell_obj[["percent.ribo"]] <- PercentageFeatureSet(tcell_obj, pattern = "^RP[SL]")
}
qc_by_cluster <- tcell_obj@meta.data %>%
  group_by(tcell_clusters) %>%
  summarise(n_cells = n(), median_nFeature = median(nFeature_RNA), median_nCount = median(nCount_RNA),
            median_pct_mt = round(median(percent.mt), 2), median_pct_ribo = round(median(percent.ribo), 2),
            .groups = "drop") %>%
  arrange(median_nFeature)
write.csv(qc_by_cluster, file.path(d$results, "validation_cluster9_qc_metrics.csv"), row.names = FALSE)
cat("\nQC metrics per cluster:\n"); print(as.data.frame(qc_by_cluster))

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
qc_long <- tcell_obj@meta.data %>%
  select(cluster = tcell_clusters, all_of(qc_metrics)) %>%
  pivot_longer(cols = all_of(qc_metrics), names_to = "metric", values_to = "value") %>%
  mutate(focus = ifelse(cluster == "9", "yes", "no"))
p_qc <- ggplot(qc_long, aes(x = cluster, y = value, fill = focus)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = c(no = "grey70", yes = "firebrick"), guide = "none") +
  facet_wrap(~metric, ncol = 2, scales = "free_y") + theme_classic() +
  ggtitle("Cluster 9 validation -- QC metrics by cluster")
ggsave(file.path(d$figures, "cluster9_qc_metrics.svg"), p_qc, width = 11, height = 7)

c9_meta <- tcell_obj@meta.data %>% filter(tcell_clusters == "9")
cat("\nCluster 9 by sample:\n"); print(table(c9_meta$sample))
cat("Cluster 9 by condition:\n"); print(table(c9_meta$condition))
cat("\nCore T-cell marker detection by cluster:\n")
print(as.data.frame(frac_expressing(tcell_obj, c("CD3D", "CD3E", "TRAC", "PTPRC"), "tcell_clusters")))
sink()
save_session(STEP, d)
log_step("Done")
