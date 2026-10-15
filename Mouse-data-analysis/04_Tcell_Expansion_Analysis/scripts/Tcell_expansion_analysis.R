#' ---
#' title: "Module 04 — T Cell Clonal Expansion Analysis"
#' author: "Apoorva Sharma"
#' ---
#'
#' Clonal expansion analysis of αβ T cells.
#' Input:  tcell_TCRab_annotated_final04_clean.rds  (Module 03 output, MAIT/innate-like NA cells removed)
#' Output: expansion-annotated object + 7 publication figures + summary CSVs
#'
#' **Clonal Definitions (critical — see v1 rationale below):**
#'
#' - **Clonotype identity:** `CTstrict` — the most specific identifier, combining
#'   α/β V(D)J(C) genes and CDR3 nucleotide sequences of all detected chains. Two cells are
#'   considered the same clone only if their full CTstrict strings match. For
#'   single-chain cells the missing chain is an "NA;NA" placeholder, so they are
#'   matched on their one detected chain.
#'
#' - **Clone size:** counted within each biological sample (`orig.ident`) SEPARATELY.
#'   This prevents cells from different samples — even of the same condition — from
#'   being pooled as a single expanding clone. A T cell that proliferates in vivo
#'   will produce daughter cells in THE SAME sample, not across samples.
#'
#' - **Expansion threshold:** >= 2 cells sharing the same CTstrict within the same
#'   sample = "Expanded". A singleton (1 cell per sample) = "Non-Expanded".
#'
#' - **Sample-averaged aggregation (used for all condition-level plots):**
#'   Percentages are computed within each individual sample FIRST, then averaged
#'   across samples within a condition. This gives each biological replicate equal
#'   weight regardless of sequencing depth. A deeply sequenced sample cannot
#'   artificially dominate the condition's signature.
#'
#' - **Heatmap vs Alluvial denominator:**
#'   - Alluvial plots: denominator = total T cells per sample
#'     → answers "What % of all T cells are expanded CD8 Teff?"
#'   - % Expanded heatmap: denominator = total cells of THAT CELL TYPE per sample
#'     → answers "What % of Tregs are clonally expanded?" (subtype-intrinsic rate)
#'   These two views are complementary, not redundant.

## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(svglite)
  library(scales)
  library(tibble)
  library(ggalluvial)
  library(ggrastr)
})

## ============================================================================
## PALETTE — matches Module 03 exactly
## ============================================================================
pal_tcrab <- c(
  "Teff CD8-like"     = "#4E79A7",
  "Th2-like"          = "#F28E2B",
  "Tfh/Tph-like"      = "#E15759",
  "Tregs-like"        = "#76B7B2",
  "Trm type mixed CD8+ CD4" = "#59A14F",
  "Naive T-like"      = "#EDC948",
  "Th17-like"         = "#B07AA1",
  "Trm/Exh-like CD8"  = "#FF9DA7",
  "Teff/Tem-like CD8" = "#9C755F",
  "Prolif T"          = "#BAB0AC"
)

condition_levels <- c("BAS", "PD", "REC", "RCU")
pal_status       <- c("Non-Expanded" = "#BFBFBF", "Expanded" = "#E74C3C")

## ============================================================================
## STAGE 1 — LOAD + COMPUTE EXPANSION STATUS
## ============================================================================

t_cells_TCRab <- readRDS("../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds")
cat("Loaded:", ncol(t_cells_TCRab), "cells,",
    length(unique(na.omit(t_cells_TCRab$celltype_merged))), "cell types\n")

clone_col          <- "CTstrict"
expansion_threshold <- 2L

meta0 <- t_cells_TCRab@meta.data %>% rownames_to_column("cell")

# Per-sample clone size: count cells sharing the same CTstrict within orig.ident
per_sample_size <- meta0 %>%
  filter(!is.na(.data[[clone_col]])) %>%
  count(orig.ident, !!sym(clone_col), name = "clone_size_in_sample")

meta1 <- meta0 %>%
  left_join(per_sample_size, by = c("orig.ident", clone_col)) %>%
  mutate(
    clone_size_in_sample = ifelse(is.na(clone_size_in_sample), 0L, as.integer(clone_size_in_sample)),
    expansion_status     = ifelse(clone_size_in_sample >= expansion_threshold, "Expanded", "Non-Expanded"),
    expansion_status     = factor(expansion_status, levels = c("Non-Expanded", "Expanded")),
    condition            = factor(condition, levels = condition_levels)
  ) %>%
  column_to_rownames("cell")

t_cells_TCRab@meta.data <- meta1

cat("\nExpansion status × condition:\n")
print(table(t_cells_TCRab$expansion_status, t_cells_TCRab$condition))
cat("\n")

## ============================================================================
## STAGE 2 — SUMMARY TABLES
## ============================================================================

# Cell-level table
cells_tbl <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, orig.ident, condition, all_of(clone_col), clone_size_in_sample, expansion_status) %>%
  arrange(orig.ident, desc(clone_size_in_sample))

write_csv(cells_tbl, "../objects/cells_expansion_table.csv")
cat("Saved: ../objects/cells_expansion_table.csv\n")

# Per-sample × clonotype summary
clono_by_sample <- per_sample_size %>%
  arrange(orig.ident, desc(clone_size_in_sample))
write_csv(clono_by_sample, "../objects/clonotype_by_sample_summary.csv")
cat("Saved: ../objects/clonotype_by_sample_summary.csv\n")

# Global clonotype summary
clono_global <- cells_tbl %>%
  filter(!is.na(.data[[clone_col]])) %>%
  group_by(.data[[clone_col]]) %>%
  summarise(
    freq_global      = n(),
    n_samples        = n_distinct(orig.ident),
    n_conditions     = n_distinct(condition),
    .groups = "drop"
  ) %>%
  arrange(desc(freq_global))
write_csv(clono_global, "../objects/clonotype_global_summary.csv")
cat("Saved: ../objects/clonotype_global_summary.csv\n")

## ============================================================================
## STAGE 3 — FIGURES
## ============================================================================
# Cell-type levels for stratified plots (NA = removed MAIT cells, excluded)
ct_levels <- names(pal_tcrab)[names(pal_tcrab) %in%
               unique(na.omit(t_cells_TCRab$celltype_merged))]

## ---- 3.1  UMAP: Expanded on top (all conditions + faceted) ------------------
umap_emb <- Embeddings(t_cells_TCRab, "umap.pca_TCRab") %>% as.data.frame()
umap_df  <- umap_emb %>%
  rownames_to_column("cell") %>%
  left_join(
    t_cells_TCRab@meta.data %>% rownames_to_column("cell") %>%
      select(cell, condition, expansion_status),
    by = "cell"
  )
xcol <- colnames(umap_emb)[1]; ycol <- colnames(umap_emb)[2]
# Negate UMAP_2 to match Module 03 orientation (purely cosmetic — embedding unchanged)
umap_df <- umap_df %>% mutate(UMAP_2_flip = -.data[[ycol]])

# All conditions combined
p_umap_all <- ggplot(umap_df, aes(.data[[xcol]], UMAP_2_flip)) +
  geom_point_rast(data = filter(umap_df, expansion_status == "Non-Expanded"),
                  color = pal_status["Non-Expanded"], alpha = 0.6, size = 0.25,
                  raster.dpi = 150) +
  geom_point_rast(data = filter(umap_df, expansion_status == "Expanded"),
                  color = pal_status["Expanded"], alpha = 0.85, size = 0.4,
                  raster.dpi = 150) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title  = element_text(face = "bold", hjust = 0.5),
        axis.line   = element_blank()) +
  labs(title = "Clonal expansion — all conditions",
       x = "UMAP_1", y = "UMAP_2")

svglite("../figures/expansion_UMAP_all.svg", width = 5, height = 4.5)
print(p_umap_all)
dev.off()
cat("Saved: ../figures/expansion_UMAP_all.svg\n")

# Faceted by condition
p_umap_cond <- ggplot(umap_df, aes(.data[[xcol]], UMAP_2_flip)) +
  geom_point_rast(data = filter(umap_df, expansion_status == "Non-Expanded"),
                  color = pal_status["Non-Expanded"], alpha = 0.6, size = 0.2,
                  raster.dpi = 150) +
  geom_point_rast(data = filter(umap_df, expansion_status == "Expanded"),
                  color = pal_status["Expanded"], alpha = 0.85, size = 0.35,
                  raster.dpi = 150) +
  facet_wrap(~ condition, nrow = 1) +
  theme_classic(base_size = 12) +
  theme(panel.grid       = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 13),
        axis.line        = element_blank()) +
  labs(x = "UMAP_1", y = "UMAP_2")

svglite("../figures/expansion_UMAP_byCondition.svg", width = 16, height = 4.2)
print(p_umap_cond)
dev.off()
cat("Saved: ../figures/expansion_UMAP_byCondition.svg\n")

## ---- 3.2  Bar: % Expanded vs Non-Expanded per condition ---------------------
# Sample-averaged: compute % within each sample, then average across samples per condition
per_sample_status <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  count(orig.ident, condition, expansion_status, name = "n") %>%
  group_by(orig.ident, condition) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

avg_status <- per_sample_status %>%
  group_by(condition, expansion_status) %>%
  summarise(
    n_samples = n_distinct(orig.ident),
    mean_pct  = mean(pct),
    sd_pct    = sd(pct),
    sem_pct   = sd_pct / sqrt(n_samples),
    .groups   = "drop"
  ) %>%
  # Non-Expanded first → bottom of bar; Expanded second → top of bar
  mutate(expansion_status = factor(expansion_status, levels = c("Non-Expanded", "Expanded")))

p_bar_pct <- ggplot(avg_status, aes(x = condition, y = mean_pct, fill = expansion_status)) +
  geom_col(position = position_fill(reverse = TRUE), color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(expansion_status == "Expanded",
                       paste0(round(mean_pct, 1), "%"), "")),
    position = position_fill(reverse = TRUE, vjust = 0.5),
    color = "white", fontface = "bold", size = 3.8
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(values = pal_status, name = NULL) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Expanded vs Non-Expanded (sample-averaged)",
       x = NULL, y = "% T cells (mean per sample)")

svglite("../figures/expansion_bar_pct.svg", width = 6, height = 5)
print(p_bar_pct)
dev.off()
cat("Saved: ../figures/expansion_bar_pct.svg\n")

# Mean counts per sample
avg_counts_status <- per_sample_status %>%
  group_by(condition, expansion_status) %>%
  summarise(
    n_samples = n_distinct(orig.ident),
    mean_n    = mean(n),
    sd_n      = sd(n),
    sem_n     = sd_n / sqrt(n_samples),
    .groups   = "drop"
  ) %>%
  mutate(expansion_status = factor(expansion_status, levels = c("Non-Expanded", "Expanded")))

p_bar_counts <- ggplot(avg_counts_status, aes(x = condition, y = mean_n, fill = expansion_status)) +
  geom_col(position = position_stack(reverse = TRUE), color = "white", linewidth = 0.3) +
  scale_fill_manual(values = pal_status, name = NULL) +
  theme_classic(base_size = 13) +
  theme(panel.grid       = element_blank(),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "none") +
  labs(title = "Mean T cell counts per sample (stacked)",
       x = "Condition", y = "Mean cell count per sample")

svglite("../figures/expansion_bar_counts.svg", width = 6, height = 5)
print(p_bar_counts)
dev.off()
cat("Saved: ../figures/expansion_bar_counts.svg\n")

## ---- 3.3  Alluvial: cell type composition of Expanded / Non-Expanded --------
# Two y-axis variants per status group:
#   (a) % of all αβ TCR+ T cells (denominator = total T cells per sample)
#   (b) Mean cell count per sample

# Denominator: total T cells per sample
den_all <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  count(orig.ident, condition, name = "n_all")

# Per sample × status × celltype counts (used for both variants)
#
# The `complete()` is REQUIRED, not cosmetic. Without it, a cell type × status
# combination absent from one replicate yields no row, so the mean() calls below
# divide by 1 instead of treating the absent replicate as 0 — inflating that
# value exactly 2x. This affected BAS only (Non-Expanded Trm/Exh-like CD8;
# Expanded Teff CD8-like; Expanded Th17-like) and made BAS's two alluvial panels
# sum to 102.14% instead of 100%. Zero-filling the full grid fixes both
# avg_pct_all and avg_count_all, which derive from this table.
num_samp <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  filter(!is.na(celltype_merged)) %>%
  mutate(celltype = factor(celltype_merged, levels = ct_levels)) %>%
  filter(!is.na(celltype)) %>%
  count(orig.ident, condition, expansion_status, celltype, name = "n_samp") %>%
  complete(nesting(orig.ident, condition), expansion_status, celltype,
           fill = list(n_samp = 0))

# (a) % of all T cells
avg_pct_all <- num_samp %>%
  left_join(den_all, by = c("orig.ident", "condition")) %>%
  mutate(pct_all = 100 * n_samp / n_all) %>%
  group_by(condition, expansion_status, celltype) %>%
  summarise(mean_pct = mean(pct_all), .groups = "drop") %>%
  mutate(condition = factor(condition, levels = condition_levels))

# (b) Mean cell count per sample
avg_count_all <- num_samp %>%
  group_by(condition, expansion_status, celltype) %>%
  summarise(mean_n = mean(n_samp), .groups = "drop") %>%
  mutate(condition = factor(condition, levels = condition_levels))

# Generic alluvial builder — y_var: column name in df; pct: add "%" suffix to axis
plot_alluvial <- function(df, status_label, title_label, y_var, y_label, pct = FALSE) {
  dat <- df %>%
    filter(expansion_status == status_label) %>%
    mutate(alluvium_id = as.integer(celltype))

  p <- ggplot(dat,
         aes(x = condition, stratum = celltype, alluvium = alluvium_id,
             y = .data[[y_var]], fill = celltype)) +
    geom_flow(alpha = 0.75, color = "white") +
    geom_stratum(width = 0.28, color = NA, alpha = 0.95) +
    scale_fill_manual(values = pal_tcrab, drop = FALSE) +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) +
    labs(title = title_label, x = NULL, y = y_label, fill = "Cell type")

  if (pct) p <- p + scale_y_continuous(labels = function(x) paste0(x, "%"))
  p
}

# ---- % versions
p_alluv_exp    <- plot_alluvial(avg_pct_all, "Expanded",
                                "Expanded clones — % of T cells",
                                y_var = "mean_pct",
                                y_label = "% of All T Cells (mean per sample)",
                                pct = TRUE)
p_alluv_nonexp <- plot_alluvial(avg_pct_all, "Non-Expanded",
                                "Non-Expanded — % of T cells",
                                y_var = "mean_pct",
                                y_label = "% of All T Cells (mean per sample)",
                                pct = TRUE)

# ---- count versions
p_alluv_exp_cnt    <- plot_alluvial(avg_count_all, "Expanded",
                                    "Expanded clones — mean cell count",
                                    y_var = "mean_n",
                                    y_label = "Mean cells per sample")
p_alluv_nonexp_cnt <- plot_alluvial(avg_count_all, "Non-Expanded",
                                    "Non-Expanded — mean cell count",
                                    y_var = "mean_n",
                                    y_label = "Mean cells per sample")

svglite("../figures/expansion_alluvial_expanded.svg",        width = 7.5, height = 5)
print(p_alluv_exp);     dev.off()
cat("Saved: ../figures/expansion_alluvial_expanded.svg\n")

svglite("../figures/expansion_alluvial_nonexpanded.svg",     width = 7.5, height = 5)
print(p_alluv_nonexp);  dev.off()
cat("Saved: ../figures/expansion_alluvial_nonexpanded.svg\n")

svglite("../figures/expansion_alluvial_expanded_counts.svg", width = 7.5, height = 5)
print(p_alluv_exp_cnt); dev.off()
cat("Saved: ../figures/expansion_alluvial_expanded_counts.svg\n")

svglite("../figures/expansion_alluvial_nonexpanded_counts.svg", width = 7.5, height = 5)
print(p_alluv_nonexp_cnt); dev.off()
cat("Saved: ../figures/expansion_alluvial_nonexpanded_counts.svg\n")

## ---- 3.4  Heatmap: % Expanded within each cell type × condition --------------
# Denominator: total cells of THAT CELL TYPE per sample (not total T cells).
# Each sample's % is computed first, then averaged across samples per condition.
# This answers "What % of CD8 Teff cells are clonally expanded?" per condition.
per_sample_ct <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  filter(!is.na(celltype_merged)) %>%
  mutate(celltype = factor(celltype_merged, levels = ct_levels)) %>%
  filter(!is.na(celltype)) %>%
  count(orig.ident, condition, celltype, expansion_status, name = "n")

pct_exp_per_sample <- per_sample_ct %>%
  group_by(orig.ident, condition, celltype) %>%
  summarise(
    pct_expanded = 100 * sum(n[expansion_status == "Expanded"], na.rm = TRUE) / sum(n),
    .groups = "drop"
  )

avg_pct_exp <- pct_exp_per_sample %>%
  group_by(condition, celltype) %>%
  summarise(mean_pct_exp = mean(pct_expanded), .groups = "drop") %>%
  mutate(condition = factor(condition, levels = condition_levels))

# Cell types ordered lowest → highest mean expansion (highest at top of plot)
ct_ord <- avg_pct_exp %>%
  group_by(celltype) %>%
  summarise(m = mean(mean_pct_exp), .groups = "drop") %>%
  arrange(m) %>%
  pull(celltype)

# Integer grid positions for geom_rect
hm_df <- avg_pct_exp %>%
  mutate(
    cx   = as.integer(factor(condition, levels = condition_levels)),
    cy   = as.integer(factor(celltype,  levels = ct_ord)),
    frac = mean_pct_exp / 100  # proportion 0-1
  )

th <- 0.44  # half tile size (slight white gap between tiles)

p_exp_heatmap <- ggplot(hm_df) +
  # Grey tile = 100% reference
  geom_rect(aes(xmin = cx - th, xmax = cx + th,
                ymin = cy - th, ymax = cy + th),
            fill = "#D8D8D8", color = "white", linewidth = 0.6) +
  # Red fill rising from bottom to pct_exp proportion of tile height
  geom_rect(aes(xmin = cx - th, xmax = cx + th,
                ymin = cy - th, ymax = cy - th + 2 * th * frac),
            fill = "#E74C3C", color = NA) +
  # % label — white text when fill > 50%, black otherwise
  geom_text(aes(x = cx, y = cy,
                label = paste0(round(mean_pct_exp, 1), "%"),
                color = ifelse(frac > 0.5, "white", "black")),
            size = 3.2, fontface = "bold") +
  scale_color_identity() +
  scale_x_continuous(breaks = seq_along(condition_levels),
                     labels = condition_levels, expand = c(0.04, 0.04)) +
  scale_y_continuous(breaks = seq_along(ct_ord),
                     labels = ct_ord,   expand = c(0.04, 0.04)) +
  theme_classic(base_size = 12) +
  theme(panel.grid  = element_blank(),
        axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(hjust = 1,     size = 11),
        plot.title  = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "% Clonally expanded within each cell type (sample-averaged)",
       x = NULL, y = NULL)

svglite("../figures/expansion_pct_heatmap.svg", width = 6, height = 10)
print(p_exp_heatmap)
dev.off()
cat("Saved: ../figures/expansion_pct_heatmap.svg\n")

## ---- 3.4b  Donut: % Expanded within each cell type × condition (fixed) ------
# Non-Expanded is derived as 100 - Expanded so every donut sums to exactly 100%.
# Values come from avg_pct_exp (per-sample first, then averaged across samples).
donut_df <- avg_pct_exp %>%
  rename(Expanded = mean_pct_exp) %>%
  mutate(`Non-Expanded` = 100 - Expanded) %>%
  pivot_longer(c(Expanded, `Non-Expanded`),
               names_to = "expansion_status", values_to = "pct") %>%
  mutate(expansion_status = factor(expansion_status, levels = c("Non-Expanded", "Expanded")),
         condition = factor(condition, levels = condition_levels))

# Label sits at midpoint of Expanded arc.
# Non-Expanded fills y 0 → (100-pct_exp), so Expanded midpoint = 100 - pct_exp/2
donut_labels <- donut_df %>%
  filter(expansion_status == "Expanded") %>%
  mutate(y_mid = 100 - pct / 2,
         label = paste0(round(pct, 1), "%"))

p_donuts <- ggplot(donut_df, aes(x = 1, y = pct, fill = expansion_status)) +
  geom_col(color = NA, width = 0.45) +
  geom_text(data = donut_labels,
            aes(x = 1, y = y_mid, label = label),
            inherit.aes = FALSE,
            color = "black", size = 2.6, fontface = "bold") +
  coord_polar(theta = "y", clip = "off") +
  xlim(0.30, 1.5) +
  facet_grid(rows = vars(condition), cols = vars(celltype),
             labeller = labeller(celltype = label_wrap_gen(width = 10))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = pal_status, name = NULL) +
  theme_void(base_size = 11) +
  theme(
    strip.text.x    = element_text(face = "bold", size = 10,
                                    margin = margin(b = 5), lineheight = 1.1),
    strip.text.y    = element_text(face = "bold", size = 11, angle = 0,
                                    margin = margin(r = 6)),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(10, "pt"),
    legend.position = "bottom"
  )

svglite("../figures/expansion_donuts.svg", width = 14, height = 6)
print(p_donuts)
dev.off()
cat("Saved: ../figures/expansion_donuts.svg\n")

## ---- 3.5  Per-condition CSVs ------------------------------------------------
for (cond in condition_levels) {
  avg_pct_exp %>%
    filter(condition == cond) %>%
    select(Cell_Type        = celltype,
           Mean_pct_expanded = mean_pct_exp) %>%
    arrange(desc(Mean_pct_expanded)) %>%
    write_csv(paste0("../objects/expansion_", cond, ".csv"))
  cat("Saved: ../objects/expansion_", cond, ".csv\n", sep = "")
}

## ---- 3.6  Expanded-pool composition (% of all EXPANDED cells) ---------------
# Denominator: all EXPANDED T cells in that condition — NOT cells of that type
# (that is 3.4 / expansion_*.csv) and NOT all T cells (that is 3.3 / the
# alluvials). This answers "of the expanded T cells at REC, what fraction are
# Tfh/Tph-like?" and is the metric behind the manuscript's 42% figure.
#
# Two aggregations are written side by side because they diverge when the two
# replicates differ in size:
#   pct_pooled  — replicates pooled, then composition computed. PREFERRED for
#                 composition: weights every cell equally.
#   pct_sampavg — composition computed per sample, then averaged. Matches the
#                 sample-averaged convention of 3.3/3.4, but gives a small
#                 replicate the same weight as a large one.
#
# NOTE: BAS has very few expanded cells (see expanded_pool_denominators.csv).
# Treat its percentages as unstable; do not report them as composition.
expanded_cells <- t_cells_TCRab@meta.data %>%
  rownames_to_column("cell") %>%
  filter(expansion_status == "Expanded", !is.na(celltype_merged)) %>%
  mutate(celltype = factor(celltype_merged, levels = ct_levels)) %>%
  filter(!is.na(celltype))

# (a) replicates pooled within condition
comp_pooled <- expanded_cells %>%
  count(condition, celltype, name = "n_expanded", .drop = FALSE) %>%
  group_by(condition) %>%
  mutate(pct_pooled = 100 * n_expanded / sum(n_expanded)) %>%
  ungroup()

# (b) per sample, then averaged across samples within condition
# `condition` is a factor, so .drop = FALSE expands orig.ident x ALL conditions;
# the filter drops those impossible sample x condition pairs, which would
# otherwise be empty groups yielding 0/0 = NaN and poisoning every mean below.
comp_sampavg <- expanded_cells %>%
  count(orig.ident, condition, celltype, name = "n", .drop = FALSE) %>%
  group_by(orig.ident, condition) %>%
  filter(sum(n) > 0) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  group_by(condition, celltype) %>%
  summarise(pct_sampavg = mean(pct), .groups = "drop")

expanded_pool_comp <- comp_pooled %>%
  left_join(comp_sampavg, by = c("condition", "celltype")) %>%
  mutate(condition = factor(condition, levels = condition_levels)) %>%
  arrange(condition, desc(pct_pooled)) %>%
  select(Condition = condition, Cell_Type = celltype,
         n_expanded, pct_pooled, pct_sampavg)

write_csv(expanded_pool_comp, "../objects/expanded_pool_composition.csv")
cat("Saved: ../objects/expanded_pool_composition.csv\n")

# The denominators themselves, so the n behind each % is never guessed at
expanded_pool_n <- expanded_cells %>%
  count(condition, name = "n_expanded_total") %>%
  mutate(condition = factor(condition, levels = condition_levels)) %>%
  arrange(condition)

write_csv(expanded_pool_n, "../objects/expanded_pool_denominators.csv")
cat("Saved: ../objects/expanded_pool_denominators.csv\n")

## ---- 3.7  Combined publication figure ---------------------------------------
# 8-column grid, canvas 18 × 15 in
# Row heights (ratio 1.06 : 0.9 : 1.06, total 3.02):
#   Row 1 = 5.3 in  | Row 2 = 4.5 in  | Row 3 = 5.3 in
#
# 12-column grid, 1.5 in/col at 18 in canvas
# Row 1: bar_pct (cols 1-4, 6 in)  bar_counts (cols 5-6, 3 in)  UMAP all (cols 7-12, 9 in)
# Row 2: UMAP by condition full width (cols 1-12, 18 in → 4 facets × 4.5 in each)
# Row 3: alluvial expanded % (cols 1-4, 6 in)  alluvial expanded counts (cols 5-12, 12 in)
p_combined <- wrap_plots(
    p_bar_pct, p_bar_counts, p_umap_all,
    p_umap_cond,
    p_alluv_exp + theme(legend.position = "none"),  # E: no legend, narrower
    p_alluv_exp_cnt,                                 # F: expanded count alluvial
    design = c(
      area(1, 1,  1, 4),   # bar_pct — 4 cols (6 in)
      area(1, 5,  1, 6),   # bar_counts — 2 cols (3 in)
      area(1, 7,  1, 12),  # UMAP all — 6 cols (9 in)
      area(2, 1,  2, 12),  # UMAP by condition — full width
      area(3, 1,  3, 6),   # alluvial expanded % — 6 cols (9 in, no legend)
      area(3, 7,  3, 12)   # alluvial expanded counts — 6 cols (9 in)
    ),
    heights = c(1.06, 0.9, 1.06)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

svglite("../figures/expansion_combined.svg", width = 18, height = 15)
print(p_combined)
dev.off()
cat("Saved: ../figures/expansion_combined.svg\n")

## ============================================================================
## STAGE 4 — SAVE ANNOTATED OBJECT
## ============================================================================
saveRDS(t_cells_TCRab, "../objects/tcell_TCRab_expansion_annotated05.rds")
cat("Saved: ../objects/tcell_TCRab_expansion_annotated05.rds\n")

cat("\n=== Module 04 complete ===\n")
cat("Objects: cells_expansion_table.csv, clonotype_by_sample_summary.csv,\n")
cat("         clonotype_global_summary.csv, expansion_BAS/PD/REC/RCU.csv,\n")
cat("         expanded_pool_composition.csv, expanded_pool_denominators.csv,\n")
cat("         tcell_TCRab_expansion_annotated05.rds\n")
cat("Figures: expansion_UMAP_all.svg, expansion_UMAP_byCondition.svg,\n")
cat("         expansion_bar_pct.svg, expansion_bar_counts.svg,\n")
cat("         expansion_alluvial_expanded.svg, expansion_alluvial_expanded_counts.svg,\n")
cat("         expansion_alluvial_nonexpanded.svg, expansion_alluvial_nonexpanded_counts.svg,\n")
cat("         expansion_pct_heatmap.svg, expansion_combined.svg\n")
