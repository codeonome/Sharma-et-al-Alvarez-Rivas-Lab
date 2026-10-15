#' ---------------------------------------------------------------------------
#' Module 06 — B-helper module score: feature UMAPs, human and mouse
#'
#' One feature UMAP per species showing a curated 15-gene B-helper module score.
#' Cluster identity is carried by the companion violin and gene-heatmap figures,
#' which share the same cluster ordering.
#'
#' IMPORTANT — this is a WITHIN-species question, not a cross-species test.
#' The two lists share 11 genes but differ in 4 each (human adds CXCL13, IL6ST,
#' PASK, AFF3; mouse adds Bcl6, Il21r, Cd84, Ebi3). Because the inputs differ the
#' two scores are NOT numerically comparable between species, and no cross-species
#' claim should be drawn from them. The divergence is forced by biology: CXCL13 is
#' expressed in 0.0% of mouse T cells. For the cross-species question see scripts
#' 01-03, which answer it negatively.
#'
#' CXCR3 was deliberately excluded from both lists. It was the only gene running
#' opposite to the target cluster in either species (human 2.6% of C9 vs 6.9%
#' overall; mouse 10.5% of Tfh/Tph-like vs 17.7% overall) -- a Th1-axis marker
#' rather than a B-helper one, dragging both scores down.
#'
#' Input:  ../../Human/tcell_final_seurat_object.rds
#'         ../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds
#' Output: ../figures/bhelper_featureumap_{human,mouse}.svg
#'         ../objects/bhelper_score_by_cluster.csv
#'         ../objects/bhelper_gene_detection.csv
#'         ../objects/bhelper_scores_percell_{human,mouse}.csv.gz
#'         ../results/sessionInfo_04.txt
#'
#' Run from 06_CrossSpecies_Tph/scripts/:  Rscript 04_bhelper_featureumap.R
#' ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
  library(patchwork); library(svglite); library(scales); library(FSA); library(tidyr)
})

set.seed(1234)

HUMAN_RDS <- "../../Human/tcell_final_seurat_object.rds"
MOUSE_RDS <- "../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds"
OBJ <- "../objects"; FIG <- "../figures"; RES <- "../results"

# --- the two curated lists, 15 genes each. CXCR3 intentionally absent. --------
BHELPER_HUMAN <- c("PDCD1","ICOS","MAF","BATF","IL21","CD200","SH2D1A","TOX2",
                   "TNFRSF4","CD4","CXCL13","IL6ST","PASK","AFF3","CD40LG")
BHELPER_MOUSE <- c("Pdcd1","Icos","Maf","Batf","Il21","Cd200","Sh2d1a","Tox2",
                   "Tnfrsf4","Cd4","Cd40lg","Bcl6","Il21r","Cd84","Ebi3")

HUMAN_TARGET <- "C9";            HUMAN_ACCENT <- "#6FA8B8"
MOUSE_TARGET <- "Tfh/Tph-like";  MOUSE_ACCENT <- "#E15759"
CTRL_N <- 100; NBIN <- 24

# Cliff's delta from the Mann-Whitney U via ranks: 2*AUC - 1.
# |d| < 0.147 negligible, < 0.33 small, < 0.474 medium, else large.
cliffs_delta <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(x, y))
  2 * ((sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)) - 1
}

cat("=====================================================================\n")
cat(" Module 06 / script 04 — B-helper module score feature UMAPs\n")
cat("=====================================================================\n\n")

# --- gene gate: CXCR3 must not creep back in --------------------------------
stopifnot(length(BHELPER_HUMAN) == 15, length(BHELPER_MOUSE) == 15)
if (any(toupper(c(BHELPER_HUMAN, BHELPER_MOUSE)) == "CXCR3"))
  stop("CXCR3 is present in a gene list but was deliberately excluded")
cat("[0] gene lists: 15 human, 15 mouse; CXCR3 confirmed absent\n")
cat("    shared (11):", paste(sort(intersect(BHELPER_HUMAN, toupper(BHELPER_MOUSE))),
                              collapse = ", "), "\n")
cat("    human-only :", paste(setdiff(BHELPER_HUMAN, toupper(BHELPER_MOUSE)),
                              collapse = ", "), "\n")
cat("    mouse-only :", paste(setdiff(toupper(BHELPER_MOUSE), BHELPER_HUMAN),
                              collapse = ", "), "\n\n")

# ---------------------------------------------------------------------------
# Score one object and build its two-panel figure
# ---------------------------------------------------------------------------
run_species <- function(obj, species, genes, group_col, reduction, target,
                        accent, n_expected, outfile, fig_w, fig_h,
                        violin_file, violin_w, violin_h,
                        heatmap_file, heatmap_file_plain, heat_w, heat_h,
                        replicate_col, pt_size, flip_y = FALSE) {

  cat(sprintf("[%s] %d cells, %d genes\n", species, ncol(obj), nrow(obj)))

  missing <- setdiff(genes, rownames(obj))
  if (length(missing))
    stop(sprintf("GENE GATE (%s): not found in object: %s",
                 species, paste(missing, collapse = ", ")))
  cat(sprintf("    gene gate passed: all %d genes present\n", length(genes)))

  n_tgt <- sum(obj@meta.data[[group_col]] == target)
  if (n_tgt != n_expected)
    stop(sprintf("IDENTITY GATE (%s): %s should be %d cells, got %d",
                 species, target, n_expected, n_tgt))
  cat(sprintf("    identity gate passed: %s = %d cells\n", target, n_tgt))

  # AddModuleScore: mean of the set minus the mean of control genes drawn from
  # matched expression bins, so highly expressed genes cannot dominate.
  obj <- AddModuleScore(obj, features = list(genes), name = "BHelper",
                        ctrl = CTRL_N, nbin = NBIN, seed = 1234)
  score_col <- "BHelper1"
  stopifnot(score_col %in% colnames(obj@meta.data))

  md <- data.frame(cell = colnames(obj),
                   cluster = as.character(obj@meta.data[[group_col]]),
                   replicate = as.character(obj@meta.data[[replicate_col]]),
                   bhelper_score = obj@meta.data[[score_col]],
                   stringsAsFactors = FALSE)

  by_cluster <- md %>%
    group_by(cluster) %>%
    summarise(n_cells = n(),
              mean_score = mean(bhelper_score),
              median_score = median(bhelper_score), .groups = "drop") %>%
    arrange(desc(mean_score)) %>%
    mutate(species = species, is_target = cluster == target, rank = row_number())
  cat("\n    score by cluster:\n")
  print(as.data.frame(by_cluster %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          select(cluster, n_cells, mean_score, median_score, rank, is_target)),
        row.names = FALSE)

  peak <- by_cluster$cluster[1]
  cat(sprintf("\n    >>> score peaks in: %s   (target = %s) -> %s\n",
              peak, target, ifelse(peak == target, "as expected", "NOT the target")))

  # per-gene detection, overall and in the target cluster
  cnt <- GetAssayData(obj, layer = "counts")[genes, , drop = FALSE]
  in_tgt <- obj@meta.data[[group_col]] == target
  detection <- data.frame(
    species = species, gene = genes,
    pct_overall = round(Matrix::rowMeans(cnt > 0) * 100, 1),
    pct_target  = round(Matrix::rowMeans(cnt[, in_tgt, drop = FALSE] > 0) * 100, 1),
    row.names = NULL) %>%
    mutate(enriched_in_target = pct_target > pct_overall) %>%
    arrange(desc(pct_target))

  # ---- figure: score UMAP | cluster reference UMAP -------------------------
  emb <- Embeddings(obj, reduction = reduction)[, 1:2]
  # Display-only vertical flip. Negating UMAP_2 rotates the embedding 180 deg
  # about the x-axis; it changes nothing about the data or the scores, and is
  # applied identically to both panels so they stay registered.
  if (flip_y) emb[, 2] <- -emb[, 2]
  df <- data.frame(x = emb[, 1], y = emb[, 2],
                   score = obj@meta.data[[score_col]],
                   cluster = factor(obj@meta.data[[group_col]]))
  ax <- colnames(emb)

  base <- theme_classic(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          legend.key.size = unit(0.35, "cm"))

  # Colour scale: light grey at the low end so background cells recede and the
  # hotspot reads immediately. Limits clipped to the 1st-99th percentile with
  # oob = squish, otherwise a handful of extreme cells compress the whole ramp
  # into its bottom fifth and nothing is visible.
  lims <- unname(quantile(df$score, c(0.01, 0.99)))
  p_score <- ggplot(df[order(df$score), ], aes(x, y, colour = score)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradientn(
      colours = c("grey90", "grey84", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000"),
      # weight the ramp so the lower ~40% stays neutral grey and only genuinely
      # elevated cells take colour -- otherwise the whole plot reads as orange
      values = scales::rescale(c(0, 0.40, 0.58, 0.74, 0.88, 1)),
      limits = lims, oob = scales::squish, name = "B-helper\nscore") +
    labs(title = paste0(species, " \u2014 B-helper module score"),
         subtitle = sprintf("%d genes, AddModuleScore (ctrl=%d, nbin=%d); colour clipped to 1st-99th pct",
                            length(genes), CTRL_N, NBIN),
         x = ax[1], y = ax[2]) + base

  svglite(file.path(FIG, outfile), width = fig_w, height = fig_h)
  print(p_score)
  invisible(dev.off())
  cat(sprintf("    saved: ../figures/%s\n", outfile))

  # ---- per-cluster statistics ---------------------------------------------
  # Dunn post-hoc, target vs every other cluster, BH-adjusted across the full
  # matrix of pairwise comparisons. Reported WITH Cliff's delta: at these cell
  # numbers p-values are near-guaranteed and the effect size is what matters.
  dn <- dunnTest(bhelper_score ~ factor(cluster), data = md, method = "bh")$res
  tgt_scores <- md$bhelper_score[md$cluster == target]
  per_cluster <- lapply(setdiff(unique(md$cluster), target), function(other) {
    oth <- md$bhelper_score[md$cluster == other]
    row <- dn[dn$Comparison %in% c(paste(target, other, sep = " - "),
                                   paste(other, target, sep = " - ")), ][1, ]
    data.frame(species = species, target = target, other_cluster = other,
               median_target = median(tgt_scores), median_other = median(oth),
               median_diff = median(tgt_scores) - median(oth),
               cliffs_delta = cliffs_delta(tgt_scores, oth),
               Z = row$Z, P_unadj = row$P.unadj, P_adj_BH = row$P.adj)
  }) %>% bind_rows() %>% arrange(desc(cliffs_delta)) %>%
    mutate(effect = cut(abs(cliffs_delta), c(-Inf, .147, .33, .474, Inf),
                        labels = c("negligible", "small", "medium", "large")),
           sig = cut(P_adj_BH, c(-Inf, .001, .01, .05, Inf),
                     labels = c("***", "**", "*", "ns")))
  cat("\n    target vs each other cluster (Dunn BH + Cliff's delta):\n")
  print(as.data.frame(per_cluster %>%
          transmute(other_cluster, median_diff = round(median_diff, 3),
                    cliffs_delta = round(cliffs_delta, 3), effect,
                    P_adj_BH = signif(P_adj_BH, 3), sig)), row.names = FALSE)

  # ---- replicate-level check ----------------------------------------------
  # Cells within a donor/sample are not independent, so the cell-level p-values
  # above are pseudo-replicated. Collapse to one mean per replicate and repeat
  # the target-vs-rest contrast as a paired test -- this is the defensible one.
  rep_means <- md %>% mutate(grp = ifelse(cluster == target, "target", "rest")) %>%
    group_by(replicate, grp) %>%
    summarise(mean_score = mean(bhelper_score), n = n(), .groups = "drop")
  rep_wide <- rep_means %>% select(-n) %>%
    pivot_wider(names_from = grp, values_from = mean_score) %>%
    filter(!is.na(target), !is.na(rest))
  rep_test <- data.frame(
    species = species, n_replicates = nrow(rep_wide),
    median_target = median(rep_wide$target), median_rest = median(rep_wide$rest),
    median_paired_diff = median(rep_wide$target - rep_wide$rest),
    p_paired_wilcoxon = if (nrow(rep_wide) >= 3)
      wilcox.test(rep_wide$target, rep_wide$rest, paired = TRUE)$p.value else NA_real_)
  cat(sprintf("\n    replicate-level (%s, n=%d): paired diff %+.3f, Wilcoxon p = %s\n",
              replicate_col, rep_test$n_replicates, rep_test$median_paired_diff,
              ifelse(is.na(rep_test$p_paired_wilcoxon), "n/a",
                     signif(rep_test$p_paired_wilcoxon, 3))))

  # ---- per-cluster distribution -------------------------------------------
  # Ordered by median so the ranking is readable at a glance. Violin + box
  # rather than a bar of means: with cluster sizes from 366 to 6,734 cells the
  # spread matters as much as the centre.
  kw <- kruskal.test(bhelper_score ~ factor(cluster), data = md)
  # kruskal.test underflows to exactly 0 at these cell numbers; report the
  # double-precision floor rather than printing a literal "p = 0"
  kw_lab <- if (kw$p.value < .Machine$double.xmin) "p < 2.2e-308"
            else sprintf("p = %.3g", kw$p.value)
  ord <- by_cluster %>% arrange(desc(median_score)) %>% pull(cluster)
  vd <- md %>% mutate(cluster = factor(cluster, levels = ord),
                      is_target = cluster == target)
  # park the n= labels in reserved space below the data, so long violin tails
  # cannot run through the text
  yr <- range(vd$bhelper_score); pad <- 0.13 * diff(yr)
  # Dunn BH significance of each cluster vs the target, placed above each violin
  siglab <- per_cluster %>%
    transmute(cluster = factor(other_cluster, levels = ord),
              sig = as.character(sig), ytop = yr[2] + pad * 0.30) %>%
    bind_rows(data.frame(cluster = factor(target, levels = ord),
                         sig = "ref", ytop = yr[2] + pad * 0.30))
  nlab <- by_cluster %>% mutate(cluster = factor(cluster, levels = ord)) %>%
    arrange(cluster) %>% mutate(lab = paste0("n=", n_cells), ypos = yr[1] - pad * 0.55)

  p_v <- ggplot(vd, aes(cluster, bhelper_score, fill = is_target)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey60") +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey35") +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.25, fill = "white") +
    geom_text(data = nlab, aes(x = cluster, y = ypos, label = lab),
              inherit.aes = FALSE, size = 2.1, colour = "grey40") +
    geom_text(data = siglab, aes(x = cluster, y = ytop, label = sig),
              inherit.aes = FALSE, size = 2.4, colour = "grey20") +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.03))) +
    scale_fill_manual(values = c(`TRUE` = accent, `FALSE` = "grey85"), guide = "none") +
    labs(title = paste0(species, " \u2014 B-helper score by cluster"),
         subtitle = sprintf("ordered by median; %s highlighted. Kruskal-Wallis %s\nmarks = Dunn BH vs %s (*** <0.001, ** <0.01, * <0.05)",
                            target, kw_lab, target),
         x = NULL, y = "B-helper module score") +
    theme_classic(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1))

  svglite(file.path(FIG, violin_file), width = violin_w, height = violin_h)
  print(p_v); invisible(dev.off())
  cat(sprintf("    saved: ../figures/%s\n", violin_file))
  cat(sprintf("    Kruskal-Wallis across clusters: %s\n\n", kw_lab))

  # ---- gene x cluster heatmap ---------------------------------------------
  # Mean log-normalised expression per cluster, z-scored WITHIN each gene, so
  # every gene's distribution across clusters is on a comparable scale and a
  # highly expressed gene (ICOS) cannot dominate a lowly expressed one (IL21).
  # This is what shows which genes actually drive the module score in the target.
  dat <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
  cl <- md$cluster
  gm <- sapply(ord, function(k) Matrix::rowMeans(dat[, cl == k, drop = FALSE]))
  gp <- sapply(ord, function(k)
          Matrix::rowMeans(GetAssayData(obj, layer = "counts")[genes, cl == k, drop = FALSE] > 0) * 100)
  gz <- t(scale(t(gm)))
  gene_ord <- names(sort(gz[, target], decreasing = FALSE))   # bottom-up in plot

  hd <- as.data.frame(as.table(gz)) %>% setNames(c("gene", "cluster", "z")) %>%
    left_join(as.data.frame(as.table(gm)) %>% setNames(c("gene", "cluster", "mean_expr")),
              by = c("gene", "cluster")) %>%
    left_join(as.data.frame(as.table(gp)) %>% setNames(c("gene", "cluster", "pct_expr")),
              by = c("gene", "cluster")) %>%
    mutate(species = species,
           gene = factor(gene, levels = gene_ord),
           cluster = factor(cluster, levels = ord),
           is_target = cluster == target)

  # Base heatmap, shared by both variants.
  heat_base <- ggplot(hd, aes(cluster, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = subset(hd, is_target), colour = "black", linewidth = 0.9, fill = NA) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "z-scored\nmean expr") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1))

  # (a) annotated: % of cells expressing printed in each tile
  p_h <- heat_base +
    geom_text(aes(label = sprintf("%.0f", pct_expr)), size = 1.9, colour = "grey25") +
    labs(title = paste0(species, " \u2014 B-helper genes by cluster"),
         subtitle = paste0("fill = mean expression z-scored within each gene; number = % of cells expressing.\n",
                           "Clusters ordered by module score (left = highest); ",
                           target, " outlined."))
  svglite(file.path(FIG, heatmap_file), width = heat_w, height = heat_h)
  print(p_h); invisible(dev.off())
  cat(sprintf("    saved: ../figures/%s\n", heatmap_file))

  # (b) clean: fill only, no in-tile numbers. Same data and same ordering --
  # the % values remain in bhelper_gene_by_cluster_expression.csv.
  p_h2 <- heat_base +
    labs(title = paste0(species, " \u2014 B-helper genes by cluster"),
         subtitle = paste0("fill = mean expression z-scored within each gene. ",
                           "Clusters ordered by module score (left = highest); ",
                           target, " outlined."))
  svglite(file.path(FIG, heatmap_file_plain), width = heat_w, height = heat_h)
  print(p_h2); invisible(dev.off())
  cat(sprintf("    saved: ../figures/%s\n", heatmap_file_plain))

  write.csv(md, gzfile(file.path(OBJ, sprintf("bhelper_scores_percell_%s.csv.gz",
                                              tolower(species)))), row.names = FALSE)
  list(by_cluster = by_cluster %>% mutate(kruskal_p = kw$p.value),
       detection = detection, per_cluster = per_cluster,
       rep_test = rep_test, heat = hd)
}

# ---------------------------------------------------------------------------
human <- readRDS(HUMAN_RDS)
H <- run_species(human, "human", BHELPER_HUMAN, "cluster_code", "umap.tcell",
                 HUMAN_TARGET, HUMAN_ACCENT, 873,
                 "bhelper_featureumap_human.svg", 6.2, 5.0,
                 violin_file = "bhelper_score_violin_human.svg",
                 violin_w = 6.5, violin_h = 4.2,
                 heatmap_file = "bhelper_gene_heatmap_human.svg",
                 heatmap_file_plain = "bhelper_gene_heatmap_human_nolabels.svg",
                 heat_w = 6.8, heat_h = 5.2,
                 replicate_col = "donor", pt_size = 0.60, flip_y = FALSE)
rm(human); invisible(gc(verbose = FALSE))

mouse <- readRDS(MOUSE_RDS)
stopifnot(ncol(mouse) == 9224)
M <- run_species(mouse, "mouse", BHELPER_MOUSE, "celltype_merged", "umap.pca_TCRab",
                 MOUSE_TARGET, MOUSE_ACCENT, 1045,
                 "bhelper_featureumap_mouse.svg", 6.2, 5.0,
                 violin_file = "bhelper_score_violin_mouse.svg",
                 violin_w = 7.5, violin_h = 4.6,
                 heatmap_file = "bhelper_gene_heatmap_mouse.svg",
                 heatmap_file_plain = "bhelper_gene_heatmap_mouse_nolabels.svg",
                 heat_w = 7.8, heat_h = 5.2,
                 replicate_col = "orig.ident", pt_size = 0.85,
                 flip_y = TRUE)   # mouse UMAP flipped vertically for display
rm(mouse); invisible(gc(verbose = FALSE))

# ---------------------------------------------------------------------------
by_cluster <- bind_rows(H$by_cluster, M$by_cluster) %>%
  select(species, cluster, n_cells, mean_score, median_score, rank, is_target, kruskal_p)
write.csv(by_cluster, file.path(OBJ, "bhelper_score_by_cluster.csv"), row.names = FALSE)

per_cluster <- bind_rows(H$per_cluster, M$per_cluster)
write.csv(per_cluster, file.path(OBJ, "bhelper_target_vs_cluster_stats.csv"), row.names = FALSE)

rep_test <- bind_rows(H$rep_test, M$rep_test)
write.csv(rep_test, file.path(OBJ, "bhelper_replicate_level_test.csv"), row.names = FALSE)
cat("\n[tables] replicate-level target-vs-rest (the non-pseudo-replicated test):\n")
print(as.data.frame(rep_test %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)

heat <- bind_rows(H$heat, M$heat) %>%
  transmute(species, gene, cluster, mean_expr = round(mean_expr, 4),
            z_within_gene = round(z, 3), pct_expr = round(pct_expr, 1))
write.csv(heat, file.path(OBJ, "bhelper_gene_by_cluster_expression.csv"), row.names = FALSE)

detection <- bind_rows(H$detection, M$detection)
write.csv(detection, file.path(OBJ, "bhelper_gene_detection.csv"), row.names = FALSE)
cat("[tables] per-gene detection (% cells expressing):\n")
print(as.data.frame(detection), row.names = FALSE)
cat(sprintf("\n    genes enriched in the target cluster: human %d/15, mouse %d/15\n",
            sum(detection$enriched_in_target[detection$species == "human"]),
            sum(detection$enriched_in_target[detection$species == "mouse"])))
stopifnot(!any(toupper(detection$gene) == "CXCR3"))

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_04.txt"))
cat("\nDone — script 04.\n")
