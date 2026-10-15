#' ---
#' title: "Global_analysis_postprocessing_v2"
#' author: "Apoorva Sharma"
#' format: html
#' editor: visual
#' ---
#'
#' **Fig 4B + Extended Data Fig 4a — v2**
#' Changes from v1:
#'   - Fig4_A_B_combined.svg: UMAP (A) + DotPlot (B) as a single SVG, panel labels A/B
#'   - DotPlot rebuilt with facet_grid so cell-type group labels appear on the right
#'   - SVG output only (no PNG)
#'
## -----------------------------------------------------------------------------

getwd()

# Load necessary libraries
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

## -----------------------------------------------------------------------------
set.seed(1234)

rna_tcr <- readRDS("../objects/rivas_3111_integrated_scrna_tcr.rds")

## -----------------------------------------------------------------------------
DefaultAssay(rna_tcr) <- "RNA"
Idents(rna_tcr) <- "cca_clusters_0.05"

## -----------------------------------------------------------------------------
# Ensure cluster factor is numeric-ordered
cl     <- levels(factor(rna_tcr$cca_clusters_0.05))
cl_num <- as.character(sort(as.numeric(as.character(cl))))
rna_tcr$cca_clusters_0.05 <- factor(rna_tcr$cca_clusters_0.05, levels = cl_num)

## -----------------------------------------------------------------------------
# PANEL A — UMAP coloured by cluster
# NRC-inspired palette: 8 visually distinct, publication-quality colors
cluster_cols <- c(
  "0" = "#4E79A7",
  "1" = "#F28E2B",
  "2" = "#E15759",
  "3" = "#76B7B2",
  "4" = "#59A14F",
  "5" = "#EDC948",
  "6" = "#B07AA1",
  "7" = "#FF9DA7"
)

# combine = FALSE returns a list; [[1]] extracts the single ggplot object
p_umap <- DimPlot(
  rna_tcr,
  reduction = "umap.cca",
  group.by  = "cca_clusters_0.05",
  label     = TRUE,
  cols      = cluster_cols,
  combine   = FALSE
)[[1]] +
  theme(legend.position = "right",
        aspect.ratio    = 1)          # force square plot panel

## -----------------------------------------------------------------------------
# Join layers before FindAllMarkers (required for Seurat v5 split layers)
if (inherits(rna_tcr, "Seurat")) rna_tcr <- JoinLayers(rna_tcr)
options(lifecycle_verbosity = "quiet")

immun_markers <- FindAllMarkers(
  rna_tcr,
  only.pos        = TRUE,
  logfc.threshold = 0.25,
  min.pct         = 0.1,
  test.use        = "wilcox"
)

top100 <- immun_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, desc(pct.1 - pct.2), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 100) %>%
  ungroup()

write.csv(top100, "../objects/cca0.05_top100_markers.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
# PANEL B — DotPlot with cell-type group labels on the right
#
# Canonical markers, ordered bottom-to-top within each group as they appear
# in the figure (Myeloid at bottom, Stromal at top).
marker_list <- list(
  "Myeloid"  = c("Cd14", "Fcgr1", "Csf3r", "C1qa", "C1qb", "Msr1", "Mafb", "Trem2", "Clec4e"),
  "T Cells"  = c("Cd3d", "Cd3e", "Cd3g", "Trac", "Trbc1", "Cd247", "Lck", "Zap70", "Tcf7", "Trbc2"),
  "B Cells"  = c("Cd19", "Cd79a", "Cd79b", "Ms4a1", "Pax5", "Ebf1", "Pou2af1", "Fcmr", "Cd22", "Bank1"),
  "Stromal"  = c("Col1a1", "Col3a1", "Dcn", "Lum", "Postn", "Pdgfra", "Fbn1", "Thbs2", "Dpt", "Col6a3")
)

# Keep only genes present in the dataset
marker_list_filt <- lapply(marker_list, function(g) intersect(g, rownames(rna_tcr)))
cat("Genes per group after filtering:\n")
print(lengths(marker_list_filt))

# Flat gene vector in display order (Myeloid bottom → Stromal top)
marker_genes <- unlist(marker_list_filt)

# Use Seurat's DotPlot to compute the scaled expression / pct values
p_dot_raw <- DotPlot(
  rna_tcr,
  features = marker_genes,
  group.by = "cca_clusters_0.05",
  cols     = c("#E6E1F4", "#3B4CC0"),
  dot.scale = 7
)

# Extract underlying data and attach group labels
dot_data <- p_dot_raw$data

gene_group_df <- do.call(rbind, lapply(names(marker_list_filt), function(grp) {
  genes <- marker_list_filt[[grp]]
  if (length(genes) == 0) return(NULL)
  data.frame(features.plot = as.character(genes), group = grp,
             stringsAsFactors = FALSE)
}))

dot_data <- merge(dot_data, gene_group_df, by = "features.plot", all.x = TRUE)

# facet_grid row order: first level = top panel, last = bottom panel
# Stromal top → B Cells → T Cells → Myeloid bottom
dot_data$group <- factor(dot_data$group,
                         levels = c("Stromal", "B Cells", "T Cells", "Myeloid"))

dot_data$features.plot <- factor(dot_data$features.plot,
                                  levels = marker_genes)  # first = bottom

# Rebuild as ggplot with facet_grid (row strip lands on the right by default)
p_dot <- ggplot(dot_data,
                aes(x = id, y = features.plot, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(
    low   = "#E6E1F4",
    high  = "#3B4CC0",
    name  = "Average\nExpression"
  ) +
  scale_size_continuous(
    name   = "Percent\nExpressed",
    range  = c(0, 7),
    breaks = c(0, 25, 50, 75)
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Marker Expression by Cluster") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 9),
    axis.text.y      = element_text(size  = 8, face = "bold"),
    strip.text.y     = element_text(face  = "bold", size = 10, angle = -90),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    strip.placement  = "outside",
    panel.spacing    = unit(0.3, "lines"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position  = "bottom",
    legend.box       = "horizontal"
  )

## -----------------------------------------------------------------------------
# COMBINE A + B with patchwork, add panel tags
combined <- p_umap + p_dot +
  plot_layout(widths = c(0.5, 0.45)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave(
  "../figures/Fig4_A_B_combined.svg",
  combined,
  width  = 12,
  height = 10,
  dpi    = 600,
  units  = "in"
)

cat("Saved: ../figures/Fig4_A_B_combined.svg\n")

## -----------------------------------------------------------------------------
# ANNOTATION — cluster → cell type mapping
map_types <- c(
  "0" = "Myeloid", "2" = "Myeloid", "5" = "Myeloid", "6" = "Myeloid",
  "1" = "T",       "4" = "T",
  "3" = "B",
  "7" = "Stromal"
)

Idents(rna_tcr) <- "cca_clusters_0.05"

rna_tcr$major_celltype <- plyr::mapvalues(
  x    = as.character(rna_tcr$cca_clusters_0.05),
  from = names(map_types),
  to   = unname(map_types)
)
rna_tcr$major_celltype <- factor(rna_tcr$major_celltype,
                                  levels = c("T", "B", "Myeloid", "Stromal"))

if (!"condition" %in% colnames(rna_tcr@meta.data)) {
  rna_tcr$condition <- factor(
    sub("_.*", "", rna_tcr$orig.ident),
    levels = c("BAS", "PD", "REC", "RCU")
  )
}

## -----------------------------------------------------------------------------
# UMAP split by condition (SVG only)
# Okabe-Ito palette — colorblind-safe, standard in biology publications
type_cols <- c(
  "T"       = "#C85A32",   # deep blue
  "B"       = "#E69F00",   # amber
  "Myeloid" =  "#0072B2",   # teal green
  "Stromal" =   "#2E9A44"  # mauve/pink
)

theme_umap_strong <- theme_minimal(base_size = 16) +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    axis.ticks      = element_blank(),
    strip.text      = element_text(face = "bold", size = 18, colour = "black"),
    plot.title      = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 15, face = "bold", colour = "black")
  )

p_split <- DimPlot(
  rna_tcr,
  reduction = "umap.cca",
  group.by  = "major_celltype",
  split.by  = "condition",
  cols      = type_cols,
  pt.size   = 1.9,
  raster    = TRUE
) +
  ggtitle("Major immune cell types by condition") +
  theme_umap_strong

ggsave("../figures/UMAP_major_types_by_condition.svg",
       p_split, width = 11, height = 5, dpi = 600, units = "in")
cat("Saved: ../figures/UMAP_major_types_by_condition.svg\n")

## -----------------------------------------------------------------------------
# PIE CHARTS (SVG only)
rna_tcr$condition     <- factor(rna_tcr$condition,     levels = c("BAS", "PD", "REC", "RCU"))
rna_tcr$major_celltype <- factor(rna_tcr$major_celltype, levels = c("T", "B", "Myeloid", "Stromal"))

pie_df <- rna_tcr@meta.data |>
  tibble::as_tibble() |>
  dplyr::group_by(condition, major_celltype) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  tidyr::complete(
    condition      = factor(levels(rna_tcr$condition),      levels = levels(rna_tcr$condition)),
    major_celltype = factor(levels(rna_tcr$major_celltype), levels = levels(rna_tcr$major_celltype)),
    fill = list(n = 0)
  ) |>
  dplyr::group_by(condition) |>
  dplyr::mutate(
    prop  = n / sum(n),
    label = ifelse(prop >= 0.03, scales::percent(prop, accuracy = 0.1), "")
  ) |>
  dplyr::ungroup()

p_pies <- ggplot(pie_df, aes(x = "", y = prop, fill = major_celltype)) +
  geom_col(width = 1, color = NA) +     # no borders between slices
  coord_polar(theta = "y") +
  facet_wrap(~ condition, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = type_cols) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 4.9, fontface = "bold", color = "black") +
  labs(title = "Cell type composition per condition") +
  theme_void(base_size = 16) +
  theme(
    strip.text       = element_text(face = "bold", size = 18),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 15, face = "bold"),
    plot.title       = element_text(face = "bold", size = 20, hjust = 0.5)
  )

ggsave("../figures/pie_celltypes_by_condition.svg",
       p_pies, width = 10, height = 4, dpi = 600, units = "in")
cat("Saved: ../figures/pie_celltypes_by_condition.svg\n")

## -----------------------------------------------------------------------------
# COMBINED: split UMAP (top) + pie charts (bottom)
p_split_noleg <- p_split +
  theme(legend.position = "none", plot.title = element_blank())

p_pies_notitle <- p_pies +
  theme(plot.title = element_blank())

combined_umap_pies <- p_split_noleg / p_pies_notitle +
  plot_layout(heights = c(1.8, 1))

ggsave("../figures/UMAP_pies_by_condition_combined.svg",
       combined_umap_pies, width = 14, height = 9, dpi = 600, units = "in")
cat("Saved: ../figures/UMAP_pies_by_condition_combined.svg\n")

## -----------------------------------------------------------------------------
# Sanity checks
assertthat::assert_that(
  sum(table(rna_tcr$cca_clusters_0.05)) == sum(table(rna_tcr$major_celltype))
)

## -----------------------------------------------------------------------------
saveRDS(rna_tcr, "../objects/global_clustered_annotated.rds")
cat("Saved annotated Seurat object.\n")

## -----------------------------------------------------------------------------
sessionInfo()
