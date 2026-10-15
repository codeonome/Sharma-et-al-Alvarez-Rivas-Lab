suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(circlize)
})


obj <- readRDS("../../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds")

# CD_broad is not stored in the object — derive from celltype_merged.
# 3-class: cluster 6 ("Trm type mixed CD8+ CD4") is a confirmed CD4/CD8-mixed
# population (see Module 03 §4.8c) — given its own "mixed" class so it is
# excluded from CD4/CD8-specific analyses but retained (labeled) in all-cell ones.
cd8_types_all <- c("Teff CD8-like", "Trm/Exh-like CD8", "Teff/Tem-like CD8")
cd4_types_all <- c("Tfh/Tph-like", "Th17-like", "Th2-like", "Tregs-like")
obj$CD_broad <- dplyr::case_when(
  obj$celltype_merged %in% cd8_types_all              ~ "CD8",
  obj$celltype_merged %in% cd4_types_all              ~ "CD4",
  obj$celltype_merged == "Trm type mixed CD8+ CD4"    ~ "mixed",
  TRUE                                                ~ NA_character_
)

conditions_all <- c("BAS", "PD", "REC", "RCU")

# ==============================================================================
# COLOR DESIGN NOTES
# ==============================================================================
# V gene colors are GLOBALLY CONSISTENT across all panels: the same TRBV gene
# gets the same color in every chord diagram. This is achieved by:
#   1. Collecting the union of top-12 TRBV genes across ALL subsets/panels
#   2. Sorting alphabetically for reproducibility
#   3. Assigning colors once from DISTINCT_V_PALETTE (30 hand-picked hex values
#      chosen for maximal hue + lightness contrast — Kelly's sequence + supplemental)
#
# GENE SELECTION (per panel):
#   Top 12 TRBV + top 12 TRBJ genes by raw cell count, computed independently
#   within each panel's subset. Only V–J pairs where BOTH genes are in their
#   respective top-12 lists are drawn.
#   Rationale: >12 sectors per side makes chord diagrams unreadable (overlapping
#   ribbons, color exhaustion, label collision). Top-12 typically captures ~80-90%
#   of V–J usage, preserving dominant repertoire structure while suppressing noise
#   from 1-2 cell pairings.
#
# LIMITATION: Per-panel top-12 is a display convenience, not a consistent
#   comparative framework. A gene used at ~8% in two panels may appear in one
#   but not the other if it ranks 12th in one and 13th in the other. This is a
#   display artifact. Do not use chord diagrams to conclude a specific V gene is
#   absent from a cell type — use the Quantification CSVs for that.
# ==============================================================================

# 30 SOLID, saturated colors arranged so adjacent entries are always different
# hues (red → blue → green → orange → purple → cyan …). Full-strength (not
# pastel) so ribbons read as distinct blocks of color. Hue alternates widely so
# whatever subset of genes appears together in a panel, neighboring sectors stay
# distinguishable. NOTE: ~20 V genes are colored in practice; 30 solid colors
# cannot all be mutually colorblind-safe (an inherent limit of many categories) —
# the Quantification CSVs, not color alone, are the source of truth for gene ID.
DISTINCT_V_PALETTE <- c(
  "#E6194B",  # red
  "#4363D8",  # blue
  "#3CB44B",  # green
  "#F58231",  # orange
  "#911EB4",  # purple
  "#2BB6D9",  # cyan
  "#F032E6",  # magenta
  "#9A6324",  # brown
  "#E6C200",  # gold
  "#469990",  # teal
  "#8CC63F",  # lime
  "#9B6EDC",  # violet
  "#1F6FB2",  # steel blue
  "#E8590C",  # burnt orange
  "#2E9E5B",  # emerald
  "#C2185B",  # rose
  "#0E9BB5",  # turquoise
  "#B5651D",  # ochre
  "#6B8E23",  # olive green
  "#7048C4",  # blue-violet
  "#D81B60",  # pink-red
  "#00897B",  # dark teal
  "#5D4037",  # dark brown
  "#3949AB",  # indigo
  "#8E24AA",  # amethyst
  "#C0CA33",  # yellow-green
  "#00838F",  # dark cyan
  "#AD1457",  # dark rose
  "#EF6C00",  # dark amber
  "#546E7A"   # blue grey
)

# ==============================================================================
# HELPERS
# ==============================================================================

clean_gene <- function(x) gsub("[^A-Za-z0-9+\\-]", "", trimws(x))

# Parse beta chains from a metadata data.frame, filter to a condition
parse_beta_for_cond <- function(meta, cond_col = "sample_category", cond_val) {
  meta %>%
    filter(.data[[cond_col]] == cond_val) %>%
    mutate(CTgene = na_if(CTgene, "NA")) %>%
    separate(CTgene, into = c("alpha_raw","beta_raw"), sep = "_", fill = "right") %>%
    separate(beta_raw, into = c("TRBV","TRBD","TRBJ","TRBC"),
             sep = "\\.", fill = "right", extra = "merge") %>%
    mutate(TRBV = clean_gene(TRBV), TRBJ = clean_gene(TRBJ)) %>%
    filter(!is.na(TRBV), TRBV != "", !is.na(TRBJ), TRBJ != "")
}

# Return the top-N TRBV gene names from any metadata slice (no condition filter).
# Used to build the global V gene → color mapping.
get_top_v_genes <- function(meta, topN_V = 12) {
  meta %>%
    mutate(CTgene = na_if(CTgene, "NA")) %>%
    separate(CTgene, into = c("alpha_raw","beta_raw"), sep = "_", fill = "right", extra = "merge") %>%
    separate(beta_raw, into = c("TRBV","TRBD","TRBJ","TRBC"),
             sep = "\\.", fill = "right", extra = "merge") %>%
    mutate(TRBV = clean_gene(TRBV)) %>%
    filter(!is.na(TRBV), TRBV != "") %>%
    count(TRBV, sort = TRUE) %>%
    slice_head(n = topN_V) %>%
    pull(TRBV)
}

# Build GLOBAL_V_COLORS: union of top-12 V genes across every subset that any
# combined figure will draw. Covers all CD4 and CD8 cell types so no panel
# ever encounters a gene missing from the map (which would render transparent).
.base_slices <- list(
  obj@meta.data,                                                              # all T cells
  obj@meta.data %>% filter(CD_broad == "CD4"),                               # all CD4
  obj@meta.data %>% filter(CD_broad == "CD4", expansion_status == "Expanded"),
  obj@meta.data %>% filter(CD_broad == "CD8"),                               # all CD8
  obj@meta.data %>% filter(CD_broad == "CD8", expansion_status == "Expanded")
)
.celltype_slices <- lapply(                                                   # each cell type expanded
  c(cd4_types_all, cd8_types_all),
  function(ct) obj@meta.data %>% filter(celltype_merged == ct, expansion_status == "Expanded")
)
.all_v_genes <- sort(unique(unlist(lapply(
  c(.base_slices, .celltype_slices), get_top_v_genes, topN_V = 12
))))
stopifnot(length(.all_v_genes) <= length(DISTINCT_V_PALETTE))
GLOBAL_V_COLORS <- setNames(DISTINCT_V_PALETTE[seq_along(.all_v_genes)], .all_v_genes)
message(sprintf("Global V gene palette: %d unique genes mapped.", length(.all_v_genes)))

# Draw ONE chord diagram into the current plot region.
# J sectors = white + dark border; ribbons colored by V gene.
# v_color_map: optional named vector (TRBV gene → color). When provided, colors
# are globally consistent across panels. Falls back to local sequential assignment.
draw_chord_panel <- function(dfc, title_label = "", topN_V = 12, topN_J = 12,
                              v_color_map = NULL) {

  if (is.null(dfc) || nrow(dfc) == 0) {
    plot.new(); title(paste0(title_label, "\n(no cells)"), cex.main = 0.9); return(invisible(NULL))
  }

  # Proportion-normalised weights
  dfc <- dfc %>% mutate(weight = 1 / nrow(dfc))
  pairs <- dfc %>%
    group_by(TRBV, TRBJ) %>%
    summarise(n = sum(weight), .groups = "drop") %>%
    arrange(desc(n))

  top_v <- dfc %>% count(TRBV, sort = TRUE) %>% slice_head(n = topN_V) %>% pull(TRBV)
  top_j <- dfc %>% count(TRBJ, sort = TRUE) %>% slice_head(n = topN_J) %>% pull(TRBJ)
  filt  <- pairs %>% filter(TRBV %in% top_v, TRBJ %in% top_j)

  if (nrow(filt) == 0) {
    plot.new(); title(paste0(title_label, "\n(no pairs)"), cex.main = 0.9); return(invisible(NULL))
  }

  beta_mat <- xtabs(n ~ TRBV + TRBJ, filt)

  # Colors ── V sectors colored, J sectors white
  n_v <- nrow(beta_mat); n_j <- ncol(beta_mat)
  if (!is.null(v_color_map)) {
    v_cols <- v_color_map[rownames(beta_mat)]
    missing <- which(is.na(v_cols))
    if (length(missing) > 0) {
      # Safety fallback: assign from palette tail for any unmapped gene
      v_cols[missing] <- DISTINCT_V_PALETTE[
        seq(length(DISTINCT_V_PALETTE), by = -1, length.out = length(missing))
      ]
    }
    names(v_cols) <- rownames(beta_mat)
  } else {
    v_cols <- setNames(DISTINCT_V_PALETTE[seq_len(n_v)], rownames(beta_mat))
  }
  grid_cols <- c(v_cols, setNames(rep("white", n_j), colnames(beta_mat)))

  # Grid borders ── V sectors: same as fill (seamless); J sectors: grey30
  grid_borders <- c(setNames(v_cols, names(v_cols)),
                    setNames(rep("grey30", n_j), colnames(beta_mat)))

  # Ribbon colors ── matrix colored by TRBV row (add slight transparency)
  ribbon_cols <- matrix(
    adjustcolor(v_cols[rownames(beta_mat)], alpha.f = 0.70),
    nrow = n_v, ncol = n_j,
    dimnames = dimnames(beta_mat)
  )

  circos.clear()
  circos.par(gap.degree = 2, track.margin = c(0.01, 0.03),
             start.degree = 90, points.overflow.warning = FALSE)

  chordDiagram(
    beta_mat,
    grid.col        = grid_cols,
    grid.border     = grid_borders,
    col             = ribbon_cols,
    order           = c(rownames(beta_mat), colnames(beta_mat)),
    annotationTrack = "grid",
    transparency    = 0,
    preAllocateTracks = list(track.height = 0.05)
  )

  circos.trackPlotRegion(track.index = 1,
    panel.fun = function(x, y) {
      s <- get.cell.meta.data("sector.index")
      circos.text(mean(get.cell.meta.data("xlim")),
                  get.cell.meta.data("ylim")[2] + mm_y(3),
                  labels = s, facing = "clockwise",
                  niceFacing = TRUE, cex = 1.0, font = 2)
    }, bg.border = NA)

  if (nchar(title_label) > 0)
    title(title_label, cex.main = 1.1, font.main = 2)

  circos.clear()
}

# ==============================================================================
# INDIVIDUAL CHORD DIAGRAMS (one SVG per condition)
# ==============================================================================

make_TCRB_chord <- function(
    obj,
    cell_filter      = NULL,
    expansion_filter = NULL,
    conditions       = conditions_all,
    topN_V = 12, topN_J = 12,
    output_prefix    = "TCRB_VJ"
) {
  obj_sub <- obj
  if (!is.null(cell_filter))
    obj_sub <- subset(obj_sub,
      cells = rownames(obj_sub@meta.data)[
        obj_sub@meta.data[[cell_filter$column]] == cell_filter$value])
  if (!is.null(expansion_filter))
    obj_sub <- subset(obj_sub,
      cells = rownames(obj_sub@meta.data)[
        obj_sub@meta.data[[expansion_filter$column]] == expansion_filter$value])

  meta_all_conds <- obj_sub@meta.data

  # Save quantification table
  quant_list <- lapply(conditions, function(cond) {
    dfc <- parse_beta_for_cond(meta_all_conds, cond_val = cond)
    if (nrow(dfc) == 0) return(NULL)
    dfc %>% count(TRBV, TRBJ, sort = TRUE) %>%
      mutate(Condition = cond, Pct_in_condition = 100 * n / sum(n))
  })
  quant_table <- bind_rows(quant_list)
  write.csv(quant_table,
            paste0("../objects/", output_prefix, "_Quantification.csv"),
            row.names = FALSE)

  for (cond in conditions) {
    dfc <- parse_beta_for_cond(meta_all_conds, cond_val = cond)
    if (nrow(dfc) == 0) { message("Skipping ", cond, " (no cells)"); next }

    svg_file <- paste0("../figures/", output_prefix, "_", cond, ".svg")
    svg(svg_file, width = 9, height = 9)
    draw_chord_panel(dfc, title_label = paste0("TCRβ V–J  |  ", cond), topN_V, topN_J)
    dev.off()
    message("Saved: ", svg_file)
  }
  invisible(quant_table)
}

# ==============================================================================
# COMBINED FIGURE: All CD4 T  ×  Expanded CD4 T  ×  4 conditions  (2 × 4 grid)
# ==============================================================================

make_combined_cd4_figure <- function(obj,
                                     conditions = conditions_all,
                                     outfile = "../figures/CD4_combined_chord.svg",
                                     topN_V = 12, topN_J = 12) {

  cd4_meta    <- obj@meta.data %>% filter(CD_broad == "CD4")
  cd4_exp_meta <- obj@meta.data %>% filter(CD_broad == "CD4",
                                            expansion_status == "Expanded")

  svg(outfile, width = 30, height = 16)
  par(mfrow = c(2, 4),
      mar   = c(3.5, 3.5, 4.0, 3.5),
      oma   = c(1, 8, 4, 1))

  row_configs <- list(
    list(meta = cd4_meta,     label = "All CD4 T"),
    list(meta = cd4_exp_meta, label = "Expanded CD4 T")
  )

  for (row_idx in seq_along(row_configs)) {
    meta_use  <- row_configs[[row_idx]]$meta
    row_label <- row_configs[[row_idx]]$label

    for (col_idx in seq_along(conditions)) {
      cond <- conditions[col_idx]
      dfc  <- parse_beta_for_cond(meta_use, cond_val = cond)
      draw_chord_panel(dfc, title_label = "", topN_V, topN_J,
                       v_color_map = GLOBAL_V_COLORS)

      # Column headers (top row only)
      if (row_idx == 1)
        mtext(cond, side = 3, line = 0.4, font = 2, cex = 1.4)
    }

    # Row label in outer left margin, centred on each row
    mtext(row_label, side = 2, outer = TRUE, line = 5,
          at = if (row_idx == 1) 0.75 else 0.25,
          font = 2, cex = 1.2, las = 0)
  }

  dev.off()
  message("Saved: ", outfile)
}

# ==============================================================================
# COMBINED FIGURE: All T cells  ×  4 conditions  (1 × 4 grid)
# ==============================================================================

make_combined_allcells_figure <- function(obj,
                                          conditions = conditions_all,
                                          outfile    = "../figures/ALLcells_combined.svg",
                                          topN_V = 12, topN_J = 12) {
  meta_all <- obj@meta.data

  svg(outfile, width = 36, height = 11)
  par(mfrow = c(1, 4),
      mar   = c(3.5, 3.5, 4.0, 3.5),
      oma   = c(1, 2, 4, 1))

  for (cond in conditions) {
    dfc <- parse_beta_for_cond(meta_all, cond_val = cond)
    draw_chord_panel(dfc, title_label = "", topN_V, topN_J,
                     v_color_map = GLOBAL_V_COLORS)
    mtext(cond, side = 3, line = 0.4, font = 2, cex = 1.4)
  }

  mtext("All T Cells  |  TCRβ V–J Usage",
        side = 3, outer = TRUE, line = 1.5, font = 2, cex = 1.6)

  dev.off()
  message("Saved: ", outfile)
}

# ==============================================================================
# COMBINED FIGURE: Expanded cells per cell type  ×  REC + RCU  (8 × 2 grid)
# ==============================================================================

make_combined_celltype_figure <- function(obj,
                                          conditions = c("REC", "RCU"),
                                          outfile    = "../figures/Expanded_celltype_combined.svg",
                                          topN_V = 12, topN_J = 12) {
  all_types <- c(cd4_types_all, cd8_types_all)
  n_types   <- length(all_types)    # 8
  n_conds   <- length(conditions)   # 2

  svg(outfile, width = 16, height = 56)
  par(mfrow = c(n_types, n_conds),
      mar   = c(3.5, 3.5, 4.0, 3.5),
      oma   = c(1, 14, 4, 1))

  for (ct in all_types) {
    meta_ct <- obj@meta.data %>%
      filter(celltype_merged == ct, expansion_status == "Expanded")

    for (cond in conditions) {
      dfc <- parse_beta_for_cond(meta_ct, cond_val = cond)
      draw_chord_panel(dfc, title_label = "", topN_V, topN_J,
                       v_color_map = GLOBAL_V_COLORS)

      if (ct == all_types[1])
        mtext(cond, side = 3, line = 0.4, font = 2, cex = 1.4)
    }
  }

  # Row labels — one per cell type, centered on each row
  at_y <- rev(seq(1 / (2 * n_types), 1 - 1 / (2 * n_types), length.out = n_types))
  mtext(all_types, side = 2, outer = TRUE, line = 6,
        at = at_y, font = 2, cex = 0.95, las = 0)

  mtext("Expanded T Cells  |  TCRβ V–J Usage  |  REC vs RCU",
        side = 3, outer = TRUE, line = 1.5, font = 2, cex = 1.6)

  dev.off()
  message("Saved: ", outfile)
}

# ==============================================================================
# RUN
# ==============================================================================

message("\n--- Combined CD4 figure (2×4) ---")
make_combined_cd4_figure(obj)

message("\n--- All T cells combined (1×4) ---")
make_combined_allcells_figure(obj)

message("\n--- Expanded cell type combined (8×2) ---")
make_combined_celltype_figure(obj)

cat("\nDone.\n")
