#' Panel D: Annotated T cell UMAP (left) + by-condition UMAPs (right)
#'
#' Combines two existing SVGs into one editable panel SVG using true inlining.
#' Output: ../figures/panel_D.svg

## ---- paths -------------------------------------------------------------------
paths <- list(
  annot = "../figures/UMAP_TCRab_annotated.svg",
  cond  = "../figures/UMAP_TCRab_byCondition.svg"
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = ", "))

## ---- SVG inlining helpers (same as assemble_figure_panels.R) ----------------
parse_svg <- function(path) {
  txt  <- paste(readLines(path, warn = FALSE), collapse = "\n")
  m    <- regexpr("viewBox=['\"]([^'\"]+)['\"]", txt)
  nums <- as.numeric(strsplit(
    trimws(gsub("viewBox=['\"]|['\"]", "", regmatches(txt, m))), "\\s+")[[1]])
  body <- gsub("(?s)^.*?<svg[^>]*>", "", txt, perl = TRUE)
  body <- gsub("(?s)</svg>\\s*$",      "", body, perl = TRUE)
  list(body = trimws(body), w = nums[3], h = nums[4])
}

prefix_ids <- function(body, pfx) {
  body <- gsub("id='([^']+)'", sprintf("id='%s_\\1'", pfx), body, perl = TRUE)
  body <- gsub('id="([^"]+)"', sprintf('id="%s_\\1"', pfx), body, perl = TRUE)
  body <- gsub("url\\(#([^)]+)\\)", sprintf("url(#%s_\\1)", pfx), body, perl = TRUE)
  body
}

panel_group <- function(path, id, dest_x, dest_y, dest_w, dest_h, align = "mid") {
  p  <- parse_svg(path)
  s  <- min(dest_w / p$w, dest_h / p$h)
  rw <- p$w * s;  rh <- p$h * s
  tx <- dest_x + (dest_w - rw) / 2
  ty <- switch(align,
    "min" = dest_y,
    "max" = dest_y + dest_h - rh,
    dest_y + (dest_h - rh) / 2
  )
  body <- prefix_ids(p$body, id)
  c(sprintf('  <g id="%s" transform="translate(%.4f,%.4f) scale(%.8f)">', id, tx, ty, s),
    body, '  </g>')
}

panel_label <- function(x, y, txt, size_pt = 40) {
  sprintf('  <text x="%.4f" y="%.4f" font-family="Arial,Helvetica,sans-serif" font-size="%s" font-weight="bold" fill="black">%s</text>',
    x + size_pt * 0.3, y + size_pt * 1.1, size_pt, txt)
}

## ---- layout (pt) ------------------------------------------------------------
# Annotated UMAP: 7×5 in  →  504×360 pt
# By-condition:  18×5 in  → 1296×360 pt
# Canvas: 22×5.5 in = 1584×396 pt

W_in <- 22;  H_in <- 5.5
W <- W_in * 72;  H <- H_in * 72   # 1584 × 396 pt

D1_w <- round(0.30 * W);  D1_x <- 0
D2_w <- W - D1_w;          D2_x <- D1_w

## ---- build ------------------------------------------------------------------
svg_lines <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"',
  sprintf('     width="%sin" height="%sin" viewBox="0 0 %s %s">', W_in, H_in, W, H),
  '',
  '  <!-- Annotated T cell UMAP -->',
  panel_group(paths$annot, "annot", D1_x, 0, D1_w, H, "mid"),
  panel_label(D1_x, 0, "D"),
  '',
  '  <!-- T cell UMAPs by condition -->',
  panel_group(paths$cond, "cond", D2_x, 0, D2_w, H, "mid"),
  '',
  '</svg>'
)

out <- "../figures/panel_D.svg"
writeLines(svg_lines, out)
cat("Saved:", out, sprintf("(%.1f × %.1f in)\n", W_in, H_in))
