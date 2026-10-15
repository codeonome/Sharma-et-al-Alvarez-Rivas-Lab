#' ---
#' title: "Myeloid_cells_analysis_v2"
#' author: "Apoorva Sharma"
#' format: html
#' editor: visual
#' ---
#'
#' **Extended Figure 7 — Myeloid subcluster analysis**
#' v2 changes:
#'   - Curated named palette (myeloid_cols) used consistently across UMAP, pies, heatmap
#'   - Pie charts in 2x2 grid (BAS/PD top, REC/RCU bottom)
#'   - Combined final figure: UMAP (A) + pies (B) + heatmap (C)
#'   - SVG output only
#'
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(scRepertoire)
  library(assertthat)
  library(ggrepel)
  library(glmGamPoi)
  library(scGate)
  library(pheatmap)
  library(grid)
})

## -----------------------------------------------------------------------------
# 1. Load the globally annotated object from Module 01
rna_tcr <- readRDS("../../01_Global_Analysis/objects/global_clustered_annotated.rds")
cat("Major cell types in global object:", levels(rna_tcr@meta.data$major_celltype), "\n")

# 2. Subset to myeloid cells
myeloid_only <- subset(rna_tcr, subset = major_celltype == "Myeloid")
myeloid_only$major_celltype <- droplevels(myeloid_only$major_celltype)
saveRDS(myeloid_only, "../objects/myeloid_only_unclustered01.RDS")
rna_tcr_myeloid <- myeloid_only

## -----------------------------------------------------------------------------
# 3. Preprocessing — Normalize, Variable Features, Scale, PCA
set.seed(123)
DefaultAssay(rna_tcr_myeloid) <- "RNA"
rna_tcr_myeloid <- NormalizeData(rna_tcr_myeloid)
rna_tcr_myeloid <- FindVariableFeatures(rna_tcr_myeloid, nfeatures = 3000)
rna_tcr_myeloid <- ScaleData(rna_tcr_myeloid, vars.to.regress = "percent.mt", verbose = FALSE)
rna_tcr_myeloid <- RunPCA(rna_tcr_myeloid, npcs = 30, verbose = FALSE, seed.use = 123)

## -----------------------------------------------------------------------------
# 4. Clustering — 20 PCs (selected for analysis)
set.seed(123)
rna_tcr_myeloid <- FindNeighbors(rna_tcr_myeloid, reduction = "pca", dims = 1:20)
rna_tcr_myeloid <- FindClusters(rna_tcr_myeloid, resolution = 0.6, cluster.name = "pca_clusters_20PC")
rna_tcr_myeloid <- RunUMAP(rna_tcr_myeloid, reduction = "pca", dims = 1:20,
                            reduction.name = "umap.pca_myeloid")

set.seed(123)
rna_tcr_myeloid <- FindClusters(rna_tcr_myeloid, resolution = 0.3,
                                 cluster.name = "pca_clusters_myeloid_0.3")
Idents(rna_tcr_myeloid) <- "pca_clusters_myeloid_0.3"
reduction_used <- "umap.pca_myeloid"

saveRDS(rna_tcr_myeloid, "../objects/myeloid_only_clustered_unfilter02.rds")

## -----------------------------------------------------------------------------
# 5. Find All Markers (for heatmap; computed on unfiltered object)
myleoid_markers <- FindAllMarkers(
  rna_tcr_myeloid,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

top500 <- myleoid_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, desc(pct.1 - pct.2), .by_group = TRUE) %>%
  slice_head(n = 500) %>%
  ungroup()

write.csv(top500, "../objects/cca0.3_top500_markers_myeloid01.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
# 6. Heatmap — average expression z-scored, computed on UNFILTERED object
gene_text <- "
CD14, TLR4, TREM1, CLEC7A, NLRP3, IL1B, LST1, TYROBP, LILRB4A, MMP9,
Csf3r, Ly6g, Cxcr2, Ifit1, Ifit3, Isg15, Oasl1, Oasl2, Rsad2, Ddx60, Cxcl10, Stat1,
Mpeg1, Cd68, Lilrb4a, Clec4e, Clec4n, Slc11a1, Hif1a, Egln3, Txnip, Ctsb,
Adgre1, Mrc1, Trem2, Mertk, Apoe, Hmox1, Fn1, Igf1, Cd36, Maf,
Clec4a1, Clec4a3, Cd302, Pld4, Itga4, Itgb7, Ccr2, Ms4a4a, Runx3, Adap1,
Camp, Lcn2, Retnlg, Ngp, Mmp8, S100a8, S100a9, Ltf, Cd177, Padi4,
Zbtb46, Flt3, Clec10a, Ly75, Ccr7, Cd209a, Cd74, H2-Eb1, H2-Aa, Batf3, Ccl17,
Clec9a, Irf8, H2-Ab1, Cd86, Cd207, Epcam, Itgae, Cldn1, Cdh1, Ciita, Tjp1,
Mcpt8, Fcer1a, Ms4a2, Il4, Il13, Cpa3, Il3ra, Gata2, Hdc, Prss34,
P2ry12, Tmem119, Fcrls, Gpr34, Sall1, Pros1, Cx3cr1, C1qa
"
genes_input <- unlist(str_split(gene_text, "[,\\n]")) %>%
  str_trim() %>% .[. != ""] %>% unique()

avg <- AverageExpression(rna_tcr_myeloid, features = genes_input,
                          assays = "RNA", slot = "data")[[1]]
z <- t(scale(t(avg), center = TRUE, scale = TRUE))
z[is.na(z)] <- 0
z[z >  2]    <-  2
z[z < -1.5]  <- -1.5
colnames(z) <- stringr::str_remove(colnames(z), "^g")

## -----------------------------------------------------------------------------
# 7. Annotation labels and CONSISTENT COLOR PALETTE
label_short <- c(
  "0"  = "Act Mon/PMN/Mac",
  "1"  = "IFN-high PMN",
  "2"  = "Act-Hyp-Mon/Mac",
  "3"  = "M2-Mac",
  "4"  = "moDC",
  "5"  = "deg_PMN",
  "6"  = "cDC2",
  "7"  = "cDC1",
  "8"  = "Lang",
  "10" = "Bas",
  "12" = "Microglia-like"
)

## -----------------------------------------------------------------------------
# 8. Filter out contamination clusters 9 (B-cell/High-mito) and 11 (T-cell doublet)
drop_ids         <- c("9", "11")
rna_myeloid_filt <- subset(rna_tcr_myeloid,
                            idents = setdiff(levels(Idents(rna_tcr_myeloid)), drop_ids))

cl_present <- levels(Idents(rna_myeloid_filt))
rna_myeloid_filt$celltype_short <- unname(label_short[as.character(Idents(rna_myeloid_filt))])
rna_myeloid_filt$celltype_short <- factor(
  rna_myeloid_filt$celltype_short,
  levels = label_short[setdiff(names(label_short), setdiff(names(label_short), cl_present))]
)
Idents(rna_myeloid_filt) <- "celltype_short"

saveRDS(rna_myeloid_filt, "../objects/myeloid_only_clustered_filter_annotated03.rds")

## -----------------------------------------------------------------------------
# Color palette — original Set3/hcl.colors (generated from filtered object levels)
labs   <- levels(Idents(rna_myeloid_filt))
n_labs <- length(labs)
pal    <- RColorBrewer::brewer.pal(max(3, min(12, n_labs)), "Set3")
if (n_labs > length(pal)) pal <- grDevices::hcl.colors(n_labs, "Dynamic")
names(pal) <- labs

## -----------------------------------------------------------------------------
# PANEL A — Annotated UMAP (no legend; shared via pies)
total_cells <- ncol(rna_myeloid_filt)

p_annot <- DimPlot(
  rna_myeloid_filt,
  reduction = reduction_used,
  group.by  = "celltype_short",
  label     = TRUE, repel = TRUE,
  pt.size   = 0.1,
  cols      = pal
) +
  coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_minimal() +
  theme(
    legend.position = "none",         # legend handled by pie panel
    legend.title    = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_blank()
  )

## -----------------------------------------------------------------------------
# PANEL B — Pie charts 2×2 (BAS/PD top row, REC/RCU bottom row)
ct_tbl <- as.data.frame(table(
  condition = rna_myeloid_filt$condition,
  celltype  = rna_myeloid_filt$celltype_short
)) %>%
  dplyr::rename(n = Freq) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    prop  = n / sum(n),
    label = ifelse(prop >= 0.03, paste0(round(prop * 100, 1), "%"), "")
  ) %>%
  dplyr::ungroup()

ct_tbl$celltype  <- factor(ct_tbl$celltype,  levels = levels(Idents(rna_myeloid_filt)))
ct_tbl$condition <- factor(ct_tbl$condition, levels = c("BAS", "PD", "REC", "RCU"))

p_pies <- ggplot(ct_tbl, aes(x = "", y = prop, fill = celltype)) +
  geom_col(width = 1, color = NA) +
  coord_polar(theta = "y") +
  facet_wrap(~ condition, nrow = 2) +   # 2×2 grid
  scale_fill_manual(values = pal) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "black") +
  theme_void(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 12),
    legend.title     = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 8)
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE,
                              override.aes = list(size = 4)))

## -----------------------------------------------------------------------------
# PANEL C — Heatmap using pheatmap (original style from lab scripts)

# Build row annotation: n_cells + Celltype color (matches UMAP palette)
get_base <- function(x) sub("^[^0-9-]*", "", x)   # strips "g" prefix if present
base_ids <- get_base(colnames(z))
size_map <- setNames(as.integer(table(Idents(rna_tcr_myeloid))),
                     names(table(Idents(rna_tcr_myeloid))))
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Celltype: map cluster ID → celltype name (NA for dropped clusters 9/11)
celltype_ann <- label_short[base_ids]

col_ann <- data.frame(
  n_cells  = sapply(base_ids, function(id) size_map[[id]] %||% 0L),
  Celltype = factor(celltype_ann, levels = labs),
  row.names = colnames(z),
  check.names = FALSE
)

ann_colors <- list(
  n_cells  = colorRampPalette(c("white", "steelblue4"))(100),
  Celltype = pal   # same palette as UMAP so numbers ↔ colors match
)

bk <- seq(min(z), max(z), length.out = 101)

ph <- pheatmap::pheatmap(
  mat                      = t(z),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method        = "ward.D2",
  show_rownames            = TRUE,
  show_colnames            = TRUE,
  border_color             = NA,
  annotation_row           = col_ann,
  annotation_colors        = ann_colors,
  angle_col                = 90,
  fontsize_col             = 7,
  breaks                   = bk,
  main                     = "Myeloid subclusters (avg exp, z-scored by gene)",
  silent                   = TRUE
)

svglite::svglite("../figures/myeloid_clusters_complex_heatmap.svg", width = 15, height = 10)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()
cat("Saved: ../figures/myeloid_clusters_complex_heatmap.svg\n")

## -----------------------------------------------------------------------------
# COMBINED FINAL FIGURE — UMAP (A) + Pies 2×2 (B) + Heatmap (C)
# pheatmap gtable → patchwork wrap_elements

p_hm_wrap <- patchwork::wrap_elements(full = ph$gtable)

top_row <- p_annot + p_pies +
  plot_layout(widths = c(1, 1.3))

combined_fig <- top_row / plot_spacer() / p_hm_wrap +
  plot_layout(heights = c(1.3, 0.15, 1.5)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave(
  "../figures/Extended_Fig7_combined.svg",
  combined_fig,
  width  = 16,
  height = 20,
  dpi    = 600,
  units  = "in"
)
cat("Saved: ../figures/Extended_Fig7_combined.svg\n")

## -----------------------------------------------------------------------------
cat("Total myeloid cells displayed:", total_cells, "\n")
print(table(Idents(rna_myeloid_filt)))

if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
writeLines(capture.output(sessionInfo()), "../results/sessionInfo.txt")
