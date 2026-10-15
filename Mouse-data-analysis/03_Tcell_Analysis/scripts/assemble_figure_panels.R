#' Assemble multi-panel T cell figure — true SVG inlining
#'
#' Each panel's SVG content is parsed and inlined directly into the combined
#' SVG. Every element (text, paths, shapes) is a real SVG node — fully
#' editable in Inkscape or Illustrator without unlinking.
#'
#' IDs are prefixed per panel to prevent conflicts.
#' Gap fix: row 2 bottom-aligned, row 3 top-aligned → UMAP and pie rows meet.
#'
#' Output: ../figures/Figure_Tcell_combined.svg

## ---- file paths (relative to THIS script location) --------------------------
paths <- list(
  global = "../../01_Global_Analysis/figures/UMAP_pies_by_condition_combined.svg",
  tcr    = "../figures/TCR_status_by_condition.svg",
  annot  = "../figures/UMAP_TCRab_annotated.svg",
  cond   = "../figures/UMAP_TCRab_byCondition.svg",
  pies   = "../figures/Pie_TCRab_byCondition.svg"
)

missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0)
  stop("Missing files: ", paste(missing, collapse = ", "))
cat("All input files found.\n")

## ---- SVG inlining helpers ---------------------------------------------------

parse_svg <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  # Extract viewBox (handles single or double quotes)
  m <- regexpr("viewBox=['\"]([^'\"]+)['\"]", txt)
  vb_str <- regmatches(txt, m)
  nums <- as.numeric(
    strsplit(trimws(gsub("viewBox=['\"]|['\"]", "", vb_str)), "\\s+")[[1]]
  )

  # Strip XML declaration and outer <svg ...> wrapper
  body <- gsub("(?s)^.*?<svg[^>]*>", "", txt, perl = TRUE)
  body <- gsub("(?s)</svg>\\s*$",      "", body, perl = TRUE)

  list(body = trimws(body), w = nums[3], h = nums[4])
}

prefix_ids <- function(body, pfx) {
  # Prefix id= definitions (single and double quotes)
  body <- gsub("id='([^']+)'", sprintf("id='%s_\\1'", pfx), body, perl = TRUE)
  body <- gsub('id="([^"]+)"', sprintf('id="%s_\\1"', pfx), body, perl = TRUE)
  # Prefix url(#...) references (clip-path, fill, mask, etc.)
  body <- gsub("url\\(#([^)]+)\\)", sprintf("url(#%s_\\1)", pfx), body, perl = TRUE)
  body
}

panel_group <- function(path, id, dest_x, dest_y, dest_w, dest_h,
                         align = "mid") {
  p  <- parse_svg(path)
  s  <- min(dest_w / p$w, dest_h / p$h)         # "meet" scale
  rw <- p$w * s;  rh <- p$h * s
  tx <- dest_x + (dest_w - rw) / 2              # center horizontally
  ty <- switch(align,
    "min" = dest_y,
    "max" = dest_y + dest_h - rh,
    dest_y + (dest_h - rh) / 2
  )
  body <- prefix_ids(p$body, id)
  c(
    sprintf('  <g id="%s" transform="translate(%.4f,%.4f) scale(%.8f)">', id, tx, ty, s),
    body,
    '  </g>'
  )
}

panel_label <- function(x, y, txt, size_pt = 40) {
  sprintf(
    '  <text x="%.4f" y="%.4f" font-family="Arial,Helvetica,sans-serif" font-size="%s" font-weight="bold" fill="black">%s</text>',
    x + size_pt * 0.3, y + size_pt * 1.1, size_pt, txt
  )
}

## ---- layout in points (1pt = 1/72 in) --------------------------------------
# Using pt as the coordinate system — matches svglite's internal viewBox units

W_in <- 24;   H_in <- 22.5
W <- W_in * 72;  H <- H_in * 72   # 1728 × 1620 pt

# Row heights (in)   →   pt
row_h_in <- c(10.5, 5.5, 6.5)
row_h    <- row_h_in * 72          # c(756, 396, 468)
row_y    <- c(0, cumsum(row_h[-length(row_h)]))

# Column widths (pt)
B_w <- 0.55 * W;  B_x <- 0
C_w <- 0.45 * W;  C_x <- B_w

D1_w <- 0.35 * W;  D1_x <- 0
D2_w <- 0.65 * W;  D2_x <- D1_w

E_w <- W;  E_x <- 0

## ---- build combined SVG -----------------------------------------------------
cat("Parsing and inlining SVGs...\n")

svg_lines <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg"',
  '     xmlns:xlink="http://www.w3.org/1999/xlink"',
  sprintf('     width="%sin" height="%sin"', W_in, H_in),
  sprintf('     viewBox="0 0 %s %s">', W, H),
  '',
  '  <!-- ═══ Panel B: Global UMAP + pies by condition ═══ -->',
  panel_group(paths$global, "panelB", B_x, row_y[1], B_w, row_h[1], "mid"),
  panel_label(B_x, row_y[1], "B"),
  '',
  '  <!-- ═══ Panel C: TCR chain pairing ═══ -->',
  panel_group(paths$tcr, "panelC", C_x, row_y[1], C_w, row_h[1], "mid"),
  panel_label(C_x, row_y[1], "C"),
  '',
  '  <!-- ═══ Panel D left: Annotated T cell UMAP ═══ -->',
  panel_group(paths$annot, "panelD1", D1_x, row_y[2], D1_w, row_h[2], "max"),
  panel_label(D1_x, row_y[2], "D"),
  '',
  '  <!-- ═══ Panel D right: T cell UMAPs by condition ═══ -->',
  panel_group(paths$cond, "panelD2", D2_x, row_y[2], D2_w, row_h[2], "max"),
  '',
  '  <!-- ═══ Panel E: Pie charts by condition ═══ -->',
  panel_group(paths$pies, "panelE", E_x, row_y[3], E_w, row_h[3], "min"),
  panel_label(E_x, row_y[3], "E"),
  '',
  '</svg>'
)

out_path <- "../figures/Figure_Tcell_combined.svg"
writeLines(svg_lines, out_path)

sz <- file.size(out_path)
cat(sprintf("Saved: %s  (%.1f MB)\n", out_path, sz / 1e6))
cat(sprintf("Canvas: %s × %.1f inches  (%s × %s pt)\n", W_in, H_in, W, H))
cat("All panel elements are inlined — fully editable in Inkscape / Illustrator.\n")
