#' Extended Figure: Heatmap (left) + Violin (right)
#'
#' True SVG inlining — all elements editable.
#' Output: ../figures/Figure_extended_heatmap_violin.svg

## ---- paths -------------------------------------------------------------------
paths <- list(
  heat = "../figures/heatmap_TCRab.svg",
  viol = "../figures/violin_markers_stacked.svg"
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing: ", paste(missing, collapse = ", "))

## ---- SVG inlining helpers ----------------------------------------------------
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

# halign: "left" = flush left, "right" = flush right, "mid" = centred
panel_group <- function(path, id, dest_x, dest_y, dest_w, dest_h,
                         valign = "mid", halign = "mid") {
  p  <- parse_svg(path)
  s  <- min(dest_w / p$w, dest_h / p$h)
  rw <- p$w * s;  rh <- p$h * s
  tx <- switch(halign,
    "left"  = dest_x,
    "right" = dest_x + dest_w - rw,
    dest_x + (dest_w - rw) / 2
  )
  ty <- switch(valign,
    "min" = dest_y,
    "max" = dest_y + dest_h - rh,
    dest_y + (dest_h - rh) / 2
  )
  body <- prefix_ids(p$body, id)
  c(sprintf('  <g id="%s" transform="translate(%.4f,%.4f) scale(%.8f)">', id, tx, ty, s),
    body, '  </g>')
}

panel_label <- function(x, y, txt, size_pt = 40) {
  sprintf(
    '  <text x="%.4f" y="%.4f" font-family="Arial,Helvetica,sans-serif" font-size="%s" font-weight="bold" fill="black">%s</text>',
    x + size_pt * 0.3, y + size_pt * 1.1, size_pt, txt
  )
}

## ---- layout -----------------------------------------------------------------
# Heatmap native size:        9.5 × 22 in
# Stacked violin native size: 8   × 14 in  (tall/narrow — scale to canvas height)
#
# Canvas height = heatmap height = 22 in (heatmap drives height)
# Stacked violin scaled to fit 22 in height → rendered width = 8 * (22/14) ≈ 12.6 in
# Canvas width = heatmap_w + gap + scaled_violin_w = 9.5 + 0.4 + 12.6 = 22.5 → 24 in

gap_in  <- 0.4
W_in    <- 24;  H_in <- 22
W       <- W_in * 72;  H <- H_in * 72

heat_w  <- 9.5 * 72   # heatmap native width in pt
viol_w  <- 12.6 * 72  # stacked violin width scaled to canvas height (8 * 22/14)
margin  <- (W - heat_w - gap_in * 72 - viol_w) / 2

A_x <- margin
A_w <- heat_w + gap_in * 72 / 2     # slot includes half the gap so heatmap right-aligns into it

B_x <- A_x + A_w
B_w <- W - B_x - margin

## ---- build ------------------------------------------------------------------
svg_lines <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"',
  sprintf('     width="%gin" height="%gin" viewBox="0 0 %g %g">', W_in, H_in, W, H),
  '',
  '  <!-- Marker heatmap (A) — right-aligned in its slot so it sits flush with violin -->',
  panel_group(paths$heat, "heatmap", A_x, 0, A_w, H, valign = "mid", halign = "right"),
  panel_label(A_x, 0, "A"),
  '',
  '  <!-- Violin marker plots (B) — left-aligned to minimise gap -->',
  panel_group(paths$viol, "violin", B_x, 0, B_w, H, valign = "mid", halign = "left"),
  panel_label(B_x, 0, "B"),
  '',
  '</svg>'
)

out <- "../figures/Figure_extended_heatmap_violin.svg"
writeLines(svg_lines, out)
cat("Saved:", out, sprintf("(%.1f × %.1f in)\n", W_in, H_in))
