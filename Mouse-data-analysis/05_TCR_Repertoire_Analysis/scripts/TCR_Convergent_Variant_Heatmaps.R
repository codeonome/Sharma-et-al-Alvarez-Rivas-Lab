#' ---
#' title: "Convergent-clonotype per-cell heatmaps (V1/V2/V3)"
#' author: "Apoorva Sharma"
#' ---
#'
#' Per-cell heatmaps of the 3 shared nucleotide variants (V1/V2/V3) of the dominant
#' convergent clonotype (protein CALGANSGTYQRF_CASSQEGQGTNERLFF), restricted to the
#' REC + RCU cells that carry them. Columns = individual cells, rows = a curated 15-gene
#' panel (Tph+, Treg+, negative markers), expression z-scored per variant (per panel),
#' blue-white-red.
#' Built with pheatmap: one gtable per variant (narrow columns via cellwidth; REC|RCU
#' separated by gaps_col; condition annotation strip), arranged in a row with gridExtra;
#' gene labels + legends shown once.
#'
#' Cell set per variant = TABLE14 rows with condition in {REC, RCU}. In REC/RCU the only
#' lineages present are the CD4 clusters + the cluster-6 "mixed" cells (all marker-confirmed
#' not-CD8 by a lineage-marker check) — no CD8/Prolif cells exist there, so the condition filter alone defines
#' the "CD4 + confirmed-not-CD8 mixed" population.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(gridExtra)
  library(grid)
  library(svglite)
})

set.seed(123)

# Case-insensitive gene-name resolver (same helper as Module 03)
resolve_gene <- function(gn, rn) {
  hits <- grep(paste0("^", gn, "$"), rn, ignore.case = TRUE, value = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

# ==============================================================================
# 1. LOAD OBJECT
# ==============================================================================
obj <- readRDS("../../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
cat("Loaded object:", ncol(obj), "cells\n")

# ==============================================================================
# 2. CELL -> VARIANT MAP (from TABLE14, REC/RCU only)
# ==============================================================================
t14 <- read.csv("../objects/Clonotype_TABLE14_convergent_variants_allcells.csv",
                stringsAsFactors = FALSE)
map <- t14 %>%
  filter(condition %in% c("REC", "RCU")) %>%
  select(cell_id, variant_id, condition) %>%
  filter(cell_id %in% colnames(obj))          # guard: keep only cells present in object

variant_levels <- c("V1", "V2", "V3")
cat("\nCells per variant (REC/RCU):\n")
print(table(factor(map$variant_id, levels = variant_levels), map$condition))

# ==============================================================================
# 3. GENE PANEL — grouped hierarchy: Tph (+), Treg (+), Negative markers
# ==============================================================================
gene_groups <- list(
  "Tph (+)"  = c("Tbx21", "Batf", "Pdcd1", "Cxcr6", "Cd4", "Maf", "Ifng", "Tnfsf11"),
  "Treg (+)" = c("Ctla4", "Foxp3"),
  "Negative" = c("Sell", "Ccr7", "Gzmb", "Cd8b1", "Cd8a")
)
gene_df <- data.frame(
  gene  = unlist(gene_groups, use.names = FALSE),
  group = rep(names(gene_groups), lengths(gene_groups)),
  stringsAsFactors = FALSE
)
gene_df$resolved <- vapply(gene_df$gene, resolve_gene, character(1), rn = rownames(obj))
missing <- gene_df$gene[is.na(gene_df$resolved)]
if (length(missing)) cat("\nWARNING genes not found (dropped):", paste(missing, collapse = ", "), "\n")
gene_df    <- gene_df[!is.na(gene_df$resolved), ]
gene_order <- gene_df$gene                                   # display names, grouped order
resolved   <- setNames(gene_df$resolved, gene_df$gene)
cat("\nResolved genes:", paste(unname(resolved), collapse = ", "), "\n")

# Row-group structure for the hierarchy (gaps_row + left Group annotation strip)
grp_lv   <- names(gene_groups)[names(gene_groups) %in% gene_df$group]
grp_sz   <- as.integer(table(factor(gene_df$group, levels = grp_lv)))
gaps_row <- head(cumsum(grp_sz), -1)
row_ann  <- data.frame(Group = factor(gene_df$group, levels = grp_lv), row.names = gene_order)
grp_cols <- c("Tph (+)" = "#4E79A7", "Treg (+)" = "#59A14F", "Negative" = "grey65")[grp_lv]

# ==============================================================================
# 4. FETCH EXPRESSION (log-normalised data layer)
# ==============================================================================
expr_all <- FetchData(obj, vars = unname(resolved), cells = map$cell_id, layer = "data")
colnames(expr_all) <- gene_order

# ==============================================================================
# 5. PER-VARIANT MATRICES (genes x cells, z-scored per variant, clamped +/-2)
# ==============================================================================
build_variant <- function(v) {
  cells_v <- map %>%
    filter(variant_id == v) %>%
    arrange(factor(condition, levels = c("REC", "RCU")), cell_id)   # REC block then RCU block
  ids <- cells_v$cell_id

  mat <- t(as.matrix(expr_all[ids, gene_order, drop = FALSE]))        # genes x cells
  z   <- t(scale(t(mat)))                                            # z-score each gene across cells
  z[is.na(z)] <- 0                                                  # constant genes -> 0 (white)
  z   <- pmax(pmin(z, 2), -2)                                       # clamp to reference's +/-2

  list(z = z, raw = mat, cond = setNames(cells_v$condition, ids), ids = ids)
}
vdat <- setNames(lapply(variant_levels, build_variant), variant_levels)

# ---- provenance CSV (matrix actually plotted) --------------------------------
prov <- bind_rows(lapply(variant_levels, function(v) {
  z <- vdat[[v]]$z; raw <- vdat[[v]]$raw
  data.frame(
    variant_id = v,
    cell_id    = rep(colnames(z), each = nrow(z)),
    condition  = rep(vdat[[v]]$cond[colnames(z)], each = nrow(z)),
    gene       = rep(rownames(z), times = ncol(z)),
    raw_expr   = round(as.vector(raw), 4),
    z          = round(as.vector(z), 4),
    row.names  = NULL
  )
}))
write.csv(prov, "../objects/Clonotype_TABLE16_convergent_variant_heatmap_matrix.csv",
          row.names = FALSE)
cat("\nSaved: ../objects/Clonotype_TABLE16_convergent_variant_heatmap_matrix.csv\n")

# ==============================================================================
# 6. FIGURE — pheatmap, one gtable per variant, arranged in a row (gridExtra)
# ==============================================================================
bwr       <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)  # blue-white-red
breaks    <- seq(-2, 2, length.out = 101)                            # symmetric, matches +/-2 clamp
cond_cols <- c("REC" = "black", "RCU" = "#CB181D")                   # REC = black, RCU = red
CELL_W    <- 9    # per-cell column width (points) — narrow, thin stripes
CELL_H    <- 14   # per-gene row height (points)

# One pheatmap gtable per variant. Narrow columns (cellwidth); REC|RCU split by a white
# gaps_col (condition strip on top). Genes shown as a grouped hierarchy: gaps_row + a left
# Group annotation strip (Tph/Treg/Negative). pheatmap forces row labels to the RIGHT, so
# gene labels + all legends are shown once on the rightmost panel (V3).
make_ph <- function(v, last) {
  z     <- vdat[[v]]$z
  cond  <- vdat[[v]]$cond[colnames(z)]
  ann   <- data.frame(Condition = factor(cond, levels = c("REC", "RCU")),
                      row.names = colnames(z))
  n_rec <- sum(cond == "REC")

  pheatmap(
    z,
    color                = bwr,
    breaks               = breaks,
    cluster_rows         = FALSE,
    cluster_cols         = FALSE,
    scale                = "none",                 # z-scoring already done per variant
    cellwidth            = CELL_W,
    cellheight           = CELL_H,
    gaps_col             = n_rec,                   # REC | RCU separator
    gaps_row             = gaps_row,                # Tph | Treg | Negative separators
    annotation_col       = ann,
    annotation_row       = row_ann,                 # left Group strip (the hierarchy)
    annotation_colors    = list(Condition = cond_cols, Group = grp_cols),
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    show_colnames        = FALSE,
    show_rownames        = last,                    # gene labels once, on the rightmost panel (V3)
    main                 = sprintf("%s (n=%d)", v, ncol(z)),
    fontsize             = 10,
    fontsize_row         = 9,
    legend               = last,                    # single Expression legend (rightmost)
    annotation_legend    = last,                    # single Group + Condition legend (rightmost)
    border_color         = NA,                       # continuous tiles, like the reference
    silent               = TRUE
  )$gtable
}

g1 <- make_ph("V1", last = FALSE)
g2 <- make_ph("V2", last = FALSE)
g3 <- make_ph("V3", last = TRUE)

svglite("../figures/convergent_variant_heatmaps.svg", width = 9.5, height = 5.5)
gridExtra::grid.arrange(
  g1, g2, g3, nrow = 1, widths = c(1.0, 0.85, 2.3),
  top = grid::textGrob(
    "Convergent clonotype (CALGANSGTYQRF / CASSQEGQGTNERLFF) — shared nucleotide variants, REC + RCU",
    gp = grid::gpar(fontface = "bold", fontsize = 11)))
dev.off()
cat("Saved: ../figures/convergent_variant_heatmaps.svg\n")

# ==============================================================================
# 7. VARIANT V1 FOCUS (highly-expanded) — standalone heatmap + TCR sequence table
# ==============================================================================
# V1 is one exact CTstrict nucleotide clone, so its TCR sequence is identical across
# all V1 cells. Parse the defining sequence and tabulate alpha / beta genes + CDR3 aa + nt.
v1_ct <- unique(t14[t14$variant_id == "V1", c("CTstrict", "CTaa", "alpha_cdr3_nt", "beta_cdr3_nt")])
stopifnot(nrow(v1_ct) == 1)                                  # V1 must be a single exact clone
a_gene <- sub(";.*", "", sub("_.*", "", v1_ct$CTstrict))     # alpha genes  (before first ';', before '_')
b_gene <- sub(";.*", "", sub(".*_", "", v1_ct$CTstrict))     # beta  genes  (after '_', before ';')
aa     <- strsplit(v1_ct$CTaa, "_", fixed = TRUE)[[1]]

seq_tbl <- data.frame(
  Chain   = c("alpha (TRA)", "beta (TRB)"),
  Genes   = c(a_gene, b_gene),
  CDR3_aa = aa,
  CDR3_nt = c(v1_ct$alpha_cdr3_nt, v1_ct$beta_cdr3_nt),
  stringsAsFactors = FALSE
)
write.csv(seq_tbl, "../objects/Clonotype_TABLE17_V1_TCR_sequence.csv", row.names = FALSE)
cat("\nSaved: ../objects/Clonotype_TABLE17_V1_TCR_sequence.csv\n")
print(seq_tbl)

# Standalone V1 heatmap (same tile size / style) with the sequence table beneath
gv1 <- make_ph("V1", last = TRUE)
tt  <- gridExtra::ttheme_minimal(
  base_size = 8,
  core    = list(fg_params = list(fontfamily = "mono", hjust = 0, x = 0.02)),
  colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02)))
tg  <- gridExtra::tableGrob(seq_tbl, rows = NULL, theme = tt)

svglite("../figures/convergent_V1_heatmap.svg", width = 9, height = 5)
gridExtra::grid.arrange(
  gv1, tg, ncol = 1, heights = c(2.2, 1),
  top = grid::textGrob("Convergent clonotype variant V1 (highly expanded) — REC + RCU",
                       gp = grid::gpar(fontface = "bold", fontsize = 11)))
dev.off()
cat("Saved: ../figures/convergent_V1_heatmap.svg\n")

cat("\n=== Convergent variant heatmaps complete ===\n")
