# =====================================================================
# Step 02 -- Global clustering: resolution sweep on the Harmony embedding
# Input  : objects/00_integrated.rds
# Output : objects/02_global_clustered.rds  (global_res.<r> columns)
#          results/02_global_clustering/cluster_quality_flags.csv
#          figures/02_global_clustering/umap_res<r>.png
# Resolution used downstream: GLOBAL_RES_FINAL (0.1, 12 clusters).
# =====================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
source("config.R"); source("R/helpers.R")
STEP <- "02_global_clustering"
d <- step_dirs(STEP)

combined <- readRDS(OBJ$integrated)

# Recompute the graph explicitly on the Harmony embedding, then cluster at
# every resolution in one call (one RNA_snn_res.<r> column each).
combined <- FindNeighbors(combined, reduction = GLOBAL_REDUCTION, dims = GLOBAL_DIMS)
combined <- FindClusters(combined, resolution = GLOBAL_RES_SWEEP, random.seed = SEED_CLUSTERS)
combined <- rename_res_cols(combined, "global_res.")
cols <- res_cols(combined, "global_res.")
print(data.frame(resolution = sub("global_res.", "", cols, fixed = TRUE),
                 n_clusters = sapply(cols, function(cc) length(unique(combined[[cc]][, 1])))))

# UMAPs at the two lowest resolutions (UMAP coordinates are fixed; only colours change)
plot_res <- function(res) {
  cc <- paste0("global_res.", res)
  df <- as.data.frame(Embeddings(combined, GLOBAL_UMAP))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$cluster <- combined[[cc]][, 1]
  centroids <- df %>% group_by(cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.5) +
    geom_text(data = centroids, aes(label = cluster), color = "black", size = 4, fontface = "bold") +
    theme_void() + theme(legend.position = "none") +
    ggtitle(sprintf("res = %s (%d clusters)", res, length(unique(df$cluster))))
  ggsave(file.path(d$figures, sprintf("umap_res%s.png", res)), p, width = 8, height = 7, dpi = 150)
}
plot_res("0.1"); plot_res("0.2")

# Tiny (< 50 cells) and single-sample clusters at each resolution
flag_summary <- lapply(cols, function(cc) {
  ct <- table(combined[[cc]][, 1], combined$sample)
  single_sample <- rownames(ct)[rowSums(ct > 0) == 1]
  data.frame(resolution = sub("global_res.", "", cc, fixed = TRUE),
             n_clusters = nrow(ct),
             n_single_sample_clusters = length(single_sample),
             single_sample_clusters   = paste(single_sample, collapse = ";"),
             n_tiny_clusters          = sum(rowSums(ct) < 50))
}) %>% bind_rows()
write.csv(flag_summary, file.path(d$results, "cluster_quality_flags.csv"), row.names = FALSE)

saveRDS(combined, OBJ$global)
save_session(STEP, d)
log_step("Done: ", OBJ$global)
