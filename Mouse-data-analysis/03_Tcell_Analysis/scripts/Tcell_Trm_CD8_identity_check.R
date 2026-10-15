suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(svglite)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(cluster)
})

set.seed(123)

# ==============================================================================
# CONFIRM IDENTITY OF THE "Trm-like CD8" CLUSTER (cluster 6 @ RNA_snn_res.0.8)
# ==============================================================================
# A lineage-marker check flagged this 781-cell cluster as lineage-ambiguous (only ~50% CD8-
# marker+, ~26% Cd4+/CD8-). Three complementary lines of evidence:
#   Line 1  — marker profile of cluster 6 vs all clusters (heatmap + dotplot)
#   Line 1b — audit the top-100 FindAllMarkers genes the original label used
#   Line 2  — per-cell lineage composition (crude >0 cross-tab + module scores)
#   Line 3  — sub-cluster cluster 6 to test whether CD4 and CD8 separate
# Scope: DIAGNOSE + RECOMMEND. No labels changed, no downstream object written.

tcr_pat <- "(^Trav|^Trbv|^Trgv|^Trdv|^Trac$|^Trbc1$|^Trbc2$|^Traj|^Trbj)"
resolve_gene <- function(gn, rn) {
  hits <- grep(paste0("^", gn, "$"), rn, ignore.case = TRUE, value = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

pal_tcrab <- c(
  "Teff CD8-like"     = "#4E79A7", "Th2-like"          = "#F28E2B",
  "Tfh/Tph-like"      = "#E15759", "Tregs-like"        = "#76B7B2",
  "Trm type mixed CD8+ CD4" = "#59A14F", "Naive T-like"      = "#EDC948",
  "Th17-like"         = "#B07AA1", "Trm/Exh-like CD8"  = "#FF9DA7",
  "Teff/Tem-like CD8" = "#9C755F", "Prolif T"          = "#BAB0AC"
)

TARGET <- "Trm type mixed CD8+ CD4"

obj <- readRDS("../objects/tcell_TCRab_annotated_final04_clean.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
cat(sprintf("Loaded %d cells; target cluster '%s': %d cells\n",
            ncol(obj), TARGET, sum(obj$celltype_merged == TARGET, na.rm = TRUE)))

# ---- marker panel (grouped) --------------------------------------------------
panel <- list(
  Lineage      = c("Cd3e", "Cd4", "Cd8a", "Cd8b1"),
  Trm          = c("Itgae", "Itga1", "Cd69", "Cxcr6", "Runx3", "Zfp683"),
  Exhaustion   = c("Pdcd1", "Havcr2", "Lag3", "Tox", "Tigit", "Ctla4"),
  Cytotoxic    = c("Gzmb", "Gzmk", "Prf1", "Nkg7", "Ccl5", "Eomes", "Ifng"),
  CD4_helper   = c("Foxp3", "Il2ra", "Cd40lg", "Icos", "Bcl6", "Cxcr5"),
  Naive_Tcm    = c("Tcf7", "Sell", "Ccr7", "Lef1")
)
# resolve to actual symbols, drop absents
panel <- lapply(panel, function(gs) {
  r <- sapply(gs, resolve_gene, rownames(obj)); unname(r[!is.na(r)])
})
gene_order <- unlist(panel, use.names = FALSE)
gene_group <- rep(names(panel), lengths(panel))
cat("Marker panel resolved:", length(gene_order), "genes\n")

ct_levels <- names(pal_tcrab)[names(pal_tcrab) %in% unique(na.omit(obj$celltype_merged))]

# ==============================================================================
# LINE 1 — marker profile: cluster 6 vs all clusters
# ==============================================================================
avg <- AverageExpression(obj, assays = "RNA", layer = "data",
                         features = gene_order, group.by = "celltype_merged")[[1]]
avg <- avg[, ct_levels, drop = FALSE]

# per-cluster mean + percent-expressing table
expr_all <- FetchData(obj, vars = gene_order, layer = "data") %>%
  rownames_to_column("cell") %>%
  left_join(obj@meta.data %>% rownames_to_column("cell") %>%
              select(cell, celltype = celltype_merged), by = "cell") %>%
  filter(!is.na(celltype))

means_tbl <- expr_all %>%
  pivot_longer(all_of(gene_order), names_to = "gene", values_to = "expr") %>%
  group_by(celltype, gene) %>%
  summarise(mean_expr = round(mean(expr), 3),
            pct_expr  = round(100 * mean(expr > 0), 1), .groups = "drop")
write.csv(means_tbl, "../objects/Trm_CD8_marker_means_by_cluster.csv", row.names = FALSE)
cat("Saved: Trm_CD8_marker_means_by_cluster.csv\n")

# heatmap (z-scored per gene)
z <- t(scale(t(as.matrix(avg))))
z[is.na(z)] <- 0; z[z > 2] <- 2; z[z < -1.5] <- -1.5
row_ann <- data.frame(Module = gene_group, row.names = gene_order)
svglite("../figures/Trm_CD8_identity_heatmap.svg", width = 9, height = 11)
pheatmap(z[gene_order, ], cluster_rows = FALSE, cluster_cols = TRUE,
         clustering_distance_cols = "correlation", clustering_method = "ward.D2",
         annotation_row = row_ann, gaps_row = cumsum(lengths(panel))[-length(panel)],
         color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
         main = "Marker profile by cell type (z-scored) — is 'Trm-like CD8' intermediate?",
         fontsize_row = 8, fontsize_col = 9, border_color = NA)
dev.off()
cat("Saved: Trm_CD8_identity_heatmap.svg\n")

# dotplot
dot <- DotPlot(obj, features = gene_order, group.by = "celltype_merged") +
  scale_color_gradientn(colors = c("yellow", "orange", "red")) +
  coord_flip() +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Marker dot plot by cell type", x = NULL, y = NULL)
svglite("../figures/Trm_CD8_identity_dotplot.svg", width = 9, height = 11)
print(dot); dev.off()
cat("Saved: Trm_CD8_identity_dotplot.svg\n")

# ==============================================================================
# LINE 1b — audit the FindAllMarkers genes the original label rested on
# ==============================================================================
mk <- read.csv("../objects/markers_TCRab_res0.8_top500.csv", stringsAsFactors = FALSE)
c6 <- mk %>% filter(cluster == 6) %>% arrange(desc(avg_log2FC)) %>%
  mutate(rank = row_number())

gd_pat <- "^Tr[gd][cvj]"  # gamma/delta TCR gene transcripts
ieg    <- c("Fos","Fosb","Egr1","Egr2","Jun","Junb","Dusp1","Tnf","Ccl4","Nr4a1","Hspa1b")
trm_g  <- c("Itgae","Itga1","Cd69","Cxcr6","Id2","Runx3","Zfp683","Jaml","Adgrg1","Cdh1")
c6 <- c6 %>% mutate(tag = case_when(
  gene %in% c("Cd8a","Cd8b1")           ~ "CD8-specific",
  gene %in% c("Cd4","Cd40lg")           ~ "CD4-specific",
  grepl(gd_pat, gene)                    ~ "gamma-delta TCR gene",
  gene %in% ieg                          ~ "immediate-early/activation",
  gene %in% trm_g                        ~ "Trm/IEL (lineage-neutral)",
  TRUE                                   ~ "other"))
write.csv(c6, "../objects/Trm_CD8_top100_marker_audit.csv", row.names = FALSE)

cat("\n=== Cluster-6 marker audit ===\n")
cat("Top 15 markers by avg_log2FC:\n")
print(c6 %>% head(15) %>% select(rank, gene, avg_log2FC, pct.1, pct.2, tag))
cat("\nRank of canonical lineage / Trm genes:\n")
key <- c("Cd8a","Cd8b1","Cd4","Itga1","Id2","Itgae","Cd69","Cxcr6")
print(c6 %>% filter(gene %in% key) %>% select(rank, gene, avg_log2FC, pct.1, pct.2))
cat("\nTag composition of top 100:\n")
print(c6 %>% head(100) %>% count(tag, name = "n_genes") %>% arrange(desc(n_genes)))
cat("Saved: Trm_CD8_top100_marker_audit.csv\n")

# ==============================================================================
# LINE 2 — per-cell lineage composition of cluster 6
# ==============================================================================
trm_cells <- colnames(obj)[which(obj$celltype_merged == TARGET)]
raw <- FetchData(obj, vars = c("Cd4", "Cd8a", "Cd8b1", "Foxp3", "Cd3e", "Mki67"),
                 cells = trm_cells, layer = "counts")
raw$cd8 <- raw$Cd8a > 0 | raw$Cd8b1 > 0
raw$cd4 <- raw$Cd4 > 0
comp <- data.frame(
  category = c("CD8+/CD4-", "CD4+/CD8-", "double-positive", "double-negative"),
  n = c(sum(raw$cd8 & !raw$cd4), sum(raw$cd4 & !raw$cd8),
        sum(raw$cd8 & raw$cd4), sum(!raw$cd8 & !raw$cd4)))
comp$pct <- round(100 * comp$n / length(trm_cells), 1)
write.csv(comp, "../objects/Trm_CD8_lineage_composition.csv", row.names = FALSE)
cat("\n=== Cluster 6 crude (>0) lineage composition ===\n"); print(comp)
cat(sprintf("Foxp3+ %.0f%% | Mki67+ %.0f%% | Cd3e+ %.0f%%\n",
            100*mean(raw$Foxp3>0), 100*mean(raw$Mki67>0), 100*mean(raw$Cd3e>0)))

p_comp <- ggplot(comp, aes(x = "Trm type mixed CD8+ CD4", y = n, fill = category)) +
  geom_col(width = 0.6, color = "white") +
  scale_fill_manual(values = c("CD8+/CD4-"="#4E79A7","CD4+/CD8-"="#E15759",
                               "double-positive"="#B07AA1","double-negative"="grey75")) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face="bold", hjust=0.5)) +
  labs(title = "Cluster 6 lineage composition (raw counts >0)", x = NULL, y = "cells", fill = NULL)
svglite("../figures/Trm_CD8_lineage_composition.svg", width = 5, height = 5)
print(p_comp); dev.off()
cat("Saved: Trm_CD8_lineage_composition.svg\n")

# ---- gamma-delta / activation / TCR-status diagnostics (cluster 6 vs object) --
gd_present <- intersect(c("Trgc2","Trgv2","Trgc1","Trdc"), rownames(obj))
gd_c6  <- FetchData(obj, vars = gd_present, cells = trm_cells, layer = "counts")
gd_all <- FetchData(obj, vars = gd_present, layer = "counts")
ieg_present <- intersect(ieg, rownames(obj))
act_c6  <- rowMeans(FetchData(obj, vars = ieg_present, cells = trm_cells, layer = "data"))
act_all <- rowMeans(FetchData(obj, vars = ieg_present, layer = "data"))
diag_tbl <- data.frame(
  metric = c(paste0("pct_", gd_present, "_pos"), "pct_ANY_gd_pos",
             "mean_activation_score", "pct_single_chain"),
  cluster6 = c(round(100*colMeans(gd_c6 > 0), 1), round(100*mean(rowSums(gd_c6>0)>0), 1),
               round(mean(act_c6), 3),
               round(100*mean(obj$tcr_status[trm_cells] == "Single chain (TRA or TRB)"), 1)),
  whole_object = c(round(100*colMeans(gd_all > 0), 1), round(100*mean(rowSums(gd_all>0)>0), 1),
                   round(mean(act_all), 3),
                   round(100*mean(obj$tcr_status == "Single chain (TRA or TRB)", na.rm=TRUE), 1))
)
write.csv(diag_tbl, "../objects/Trm_CD8_gd_activation_diagnostics.csv", row.names = FALSE)
cat("\n=== Cluster 6: gamma-delta / activation / TCR-status vs whole object ===\n")
print(diag_tbl)
cat("Saved: Trm_CD8_gd_activation_diagnostics.csv\n")

# dropout-robust module scores (CD4 vs CD8), computed on whole object, read on cluster 6
# CD4 set is lineage-specific only (no pan-T Cd3e / no activation Icos, which would inflate it)
cd8_set <- list(intersect(c("Cd8a","Cd8b1"), rownames(obj)))
cd4_set <- list(intersect(c("Cd4","Cd40lg"), rownames(obj)))
obj <- AddModuleScore(obj, features = cd8_set, name = "CD8score", seed = 123)
obj <- AddModuleScore(obj, features = cd4_set, name = "CD4score", seed = 123)
ms <- obj@meta.data[trm_cells, c("CD8score1", "CD4score1")] %>%
  rownames_to_column("cell") %>%
  mutate(call = case_when(CD8score1 - CD4score1 >  0.05 ~ "CD8",
                          CD4score1 - CD8score1 >  0.05 ~ "CD4",
                          TRUE ~ "ambiguous"))
write.csv(ms, "../objects/Trm_CD8_module_scores.csv", row.names = FALSE)
cat("\n=== Cluster 6 module-score classification (dropout-robust) ===\n")
print(ms %>% count(call, name = "n_cells") %>% mutate(pct = round(100*n_cells/nrow(ms),1)))

p_ms <- ggplot(ms, aes(CD8score1, CD4score1, color = call)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_color_manual(values = c("CD8"="#4E79A7","CD4"="#E15759","ambiguous"="grey65")) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face="bold", hjust=0.5)) +
  labs(title = "Cluster 6: CD8 vs CD4 module score (per cell)",
       x = "CD8 module score", y = "CD4 module score", color = NULL)
svglite("../figures/Trm_CD8_module_scores.svg", width = 6, height = 5)
print(p_ms); dev.off()
cat("Saved: Trm_CD8_module_scores.svg\n")

# ==============================================================================
# LINE 3 — sub-cluster cluster 6 (cluster-2 template)
# ==============================================================================
cl6 <- subset(obj, cells = trm_cells)
cl6 <- JoinLayers(cl6)
DefaultAssay(cl6) <- "RNA"
cl6 <- NormalizeData(cl6, verbose = FALSE)
cl6 <- FindVariableFeatures(cl6, nfeatures = 2000, verbose = FALSE)
VariableFeatures(cl6) <- setdiff(VariableFeatures(cl6),
                                 grep(tcr_pat, VariableFeatures(cl6), value = TRUE, ignore.case = TRUE))
cl6 <- ScaleData(cl6, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
cl6 <- RunPCA(cl6, npcs = 30, seed.use = 123, verbose = FALSE)
cl6 <- RunUMAP(cl6, reduction = "pca", dims = 1:10, reduction.name = "umap.cl6", seed.use = 123, verbose = FALSE)
cl6 <- FindNeighbors(cl6, reduction = "pca", dims = 1:10, k.param = 15, verbose = FALSE)
for (res in c(0.2, 0.4, 0.6, 0.8)) {
  col_nm <- paste0("cl6_res", gsub("\\.", "", sprintf("%.1f", res)))
  cl6 <- FindClusters(cl6, resolution = res, algorithm = 1, cluster.name = col_nm, verbose = FALSE)
  cat(sprintf("  res %.1f: %d sub-clusters\n", res, length(unique(cl6@meta.data[[col_nm]]))))
}
pca_cl6 <- cl6@reductions$pca@cell.embeddings[, 1:10]
dist_cl6 <- dist(pca_cl6)
compute_sil <- function(labels, d) {
  if (length(unique(labels)) < 2) return(0)
  mean(silhouette(as.integer(factor(labels)), d)[, 3])
}
sil <- sapply(c("0.2","0.4","0.6","0.8"), function(r)
  compute_sil(cl6@meta.data[[paste0("cl6_res", gsub("\\.","",r))]], dist_cl6))
cat("  Silhouette:", paste(names(sil), round(sil,3), sep="=", collapse=", "), "\n")
best_res <- names(which.max(sil))
best_col <- paste0("cl6_res", gsub("\\.", "", best_res))
cl6$subcluster <- cl6@meta.data[[best_col]]
cat("  Best resolution:", best_res, "\n")

sub_markers <- FetchData(cl6, vars = c("Cd4","Cd8a","Cd8b1","Itgae","Cxcr6","Gzmk","Foxp3"),
                         layer = "data") %>%
  rownames_to_column("cell") %>%
  left_join(cl6@meta.data %>% rownames_to_column("cell") %>% select(cell, subcluster), by = "cell") %>%
  pivot_longer(-c(cell, subcluster), names_to = "gene", values_to = "expr") %>%
  group_by(subcluster, gene) %>%
  summarise(mean_expr = round(mean(expr),3), pct_expr = round(100*mean(expr>0),1), .groups = "drop")
n_by_sub <- cl6@meta.data %>% count(subcluster, name = "n_cells")
write.csv(sub_markers, "../objects/Trm_CD8_subcluster_markers.csv", row.names = FALSE)
cat("\n=== Sub-cluster sizes ===\n"); print(n_by_sub)
cat("\n=== Sub-cluster Cd4/Cd8 profile (pct expressing) — do lineages separate? ===\n")
print(sub_markers %>% filter(gene %in% c("Cd4","Cd8a","Cd8b1")) %>%
        select(subcluster, gene, pct_expr) %>%
        pivot_wider(names_from = gene, values_from = pct_expr))
# what DOES distinguish the sub-clusters? top DE genes between them
Idents(cl6) <- "subcluster"
if (length(unique(cl6$subcluster)) >= 2) {
  sub_de <- FindAllMarkers(cl6, only.pos = TRUE, logfc.threshold = 0.5,
                           min.pct = 0.25, verbose = FALSE)
  sub_de_top <- sub_de %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 15) %>% ungroup()
  write.csv(sub_de_top, "../objects/Trm_CD8_subcluster_top_DE.csv", row.names = FALSE)
  cat("\n=== Top genes distinguishing the sub-clusters (what axis separates them) ===\n")
  print(sub_de_top %>% group_by(cluster) %>%
          summarise(top_genes = paste(head(gene, 10), collapse = ", "), .groups = "drop"))
  cat("Saved: Trm_CD8_subcluster_top_DE.csv\n")
}

# sub-cluster UMAP + Cd4/Cd8 feature UMAPs
emb <- Embeddings(cl6, "umap.cl6") %>% as.data.frame() %>% rownames_to_column("cell")
xg <- colnames(emb)[2]; yg <- colnames(emb)[3]
fd <- FetchData(cl6, vars = c("Cd4","Cd8a","Cd8b1"), layer = "data") %>% rownames_to_column("cell")
df <- emb %>%
  left_join(cl6@meta.data %>% rownames_to_column("cell") %>% select(cell, subcluster), by="cell") %>%
  left_join(fd, by = "cell")
base_theme <- theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        plot.title = element_text(face="bold", hjust=0.5))
p_sub <- ggplot(df, aes(.data[[xg]], .data[[yg]], color = subcluster)) +
  geom_point(size=1) + base_theme + labs(title="Sub-clusters", x="UMAP1", y="UMAP2", color="sub")
feat_umap <- function(g) ggplot(df[order(df[[g]]),], aes(.data[[xg]], .data[[yg]], color=.data[[g]])) +
  geom_point(size=1) +
  scale_color_gradientn(colors=c("grey90","#FEE5D9","#FC9272","#DE2D26","#67000D")) +
  base_theme + labs(title=g, x="UMAP1", y="UMAP2", color="expr")
combined <- p_sub + feat_umap("Cd8a") + feat_umap("Cd8b1") + feat_umap("Cd4") +
  plot_layout(nrow = 2)
svglite("../figures/Trm_CD8_subcluster_umap.svg", width = 11, height = 9)
print(combined); dev.off()
cat("Saved: Trm_CD8_subcluster_umap.svg\n")

cat("\n=== Trm-like CD8 identity check complete ===\n")
