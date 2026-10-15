#' ---------------------------------------------------------------------------
#' Module 06 — Cross-species Tph: reciprocal module scoring
#'
#' Score every human T cell for the mouse-derived Tph signature (ortholog-mapped)
#' and every mouse T cell for the human-derived one, plus a shared cross-species
#' anchor set and an orthogonal chronic-activation set. Test whether the Tph
#' cluster scores higher than the other clusters within its own dataset.
#'
#' Unlike script 01 (which compares whole marker LISTS), this works at the level
#' of CELLS and tests a hand-curated program. The two can legitimately disagree.
#'
#' Design — every set is scored in both species:
#'   human cells: TPH_MOUSE = cross-species probe | TPH_HUMAN = within-species
#'                positive control | ANCHOR = shared | CHRONIC = orthogonal
#'   mouse cells: TPH_HUMAN = cross-species probe | TPH_MOUSE = within-species
#'                positive control | ANCHOR = shared | CHRONIC = orthogonal
#'
#' Input:  ../../Human/human_mouse_Tph_chronic_activation_gene_table.csv
#'         ../../Human/tcell_final_seurat_object.rds
#'         ../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds
#' Output: ../objects/gene_sets_parsed.csv, gene_set_coverage.csv
#'         ../objects/module_scores_percell_{human,mouse}.csv.gz
#'         ../objects/kruskal_results.csv, dunn_results.csv
#'         ../objects/target_vs_rest_summary.csv
#'         ../objects/score_method_concordance.csv
#'         ../objects/donor_level_human.csv, sample_level_means_mouse.csv
#'         ../objects/chronic_by_condition_in_tph.csv
#'         ../figures/*.svg
#'         ../results/sessionInfo_02.txt
#'
#' Run from 06_CrossSpecies_Tph/scripts/:  Rscript 02_crossspecies_module_score.R
#' ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(FSA); library(UCell); library(svglite); library(patchwork)
})

set.seed(1234)

GENE_TABLE <- "../../Human/human_mouse_Tph_chronic_activation_gene_table.csv"
HUMAN_RDS  <- "../../Human/tcell_final_seurat_object.rds"
MOUSE_RDS  <- "../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds"
OBJ <- "../objects"; FIG <- "../figures"; RES <- "../results"

HUMAN_TARGET <- "C9"; MOUSE_TARGET <- "Tfh/Tph-like"
HUMAN_ACCENT <- "#6FA8B8"   # matches human_tcell_final_annotation_v4.R palette for C9
MOUSE_ACCENT <- "#E15759"   # matches Tcell_ab_analysis_reclustering.R for cluster 2
CTRL_N <- 100; NBIN <- 24
MIN_C9_CELLS_PER_DONOR <- 10

cat("=====================================================================\n")
cat(" Module 06 / script 02 — reciprocal cross-species module scoring\n")
cat("=====================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Parse the curated gene table
#
# Two rows (CD44, CFLAR) carry an unquoted comma inside `note`, which makes
# read.csv() spill them into phantom rows. Split on the FIRST FOUR commas only.
# ---------------------------------------------------------------------------
cat("[1] Parsing curated gene table\n")
raw <- readLines(GENE_TABLE)
raw <- raw[nzchar(trimws(raw))][-1]                 # drop header and blank lines
pat <- "^([^,]*),([^,]*),([^,]*),([^,]*),?(.*)$"
gt <- data.frame(
  human_gene       = sub(pat, "\\1", raw),
  mouse_gene       = sub(pat, "\\2", raw),
  category         = sub(pat, "\\3", raw),
  specificity_flag = sub(pat, "\\4", raw),
  note             = sub(pat, "\\5", raw),
  stringsAsFactors = FALSE
)

# --- parse gate --------------------------------------------------------------
stopifnot(nrow(gt) == 81)
stopifnot(sum(duplicated(gt$human_gene)) == 0, sum(duplicated(gt$mouse_gene)) == 0)
stopifnot(all(gt$specificity_flag %in% c("high", "medium", "low")))
expected_cats <- c("cross_species_anchor", "Tph_human_specific", "Tph_mouse_specific",
                   "chronic_activation", "Tph_human_specific_and_chronic_activation")
stopifnot(all(gt$category %in% expected_cats))
cat("    parse gate passed: 81 rows, 0 duplicate symbols\n")
print(table(gt$category))
write.csv(gt, file.path(OBJ, "gene_sets_parsed.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Build the four gene sets in both symbol spaces, in two filter variants
#
# PDCD4 is tagged Tph_human_specific_and_chronic_activation and so belongs to
# both TPH_HUMAN and CHRONIC by design.
# ---------------------------------------------------------------------------
BOTH <- "Tph_human_specific_and_chronic_activation"
set_rows <- list(
  ANCHOR    = gt$category == "cross_species_anchor",
  TPH_HUMAN = gt$category %in% c("Tph_human_specific", BOTH),
  TPH_MOUSE = gt$category == "Tph_mouse_specific",
  CHRONIC   = gt$category %in% c("chronic_activation", BOTH)
)
LOW_FLAG <- gt$human_gene[gt$specificity_flag == "low"]

build_sets <- function(variant) {
  keep <- if (variant == "strict") gt$specificity_flag != "low" else rep(TRUE, nrow(gt))
  lapply(set_rows, function(sel) list(human = gt$human_gene[sel & keep],
                                      mouse = gt$mouse_gene[sel & keep]))
}
SETS <- list(full = build_sets("full"), strict = build_sets("strict"))

cat("\n[2] Gene sets (n genes)\n")
for (v in names(SETS)) for (s in names(SETS[[v]]))
  cat(sprintf("    %-6s %-10s %d\n", v, s, length(SETS[[v]][[s]]$human)))
cat(sprintf("\n    'low' specificity_flag genes dropped by strict (%d): %s\n",
            length(LOW_FLAG), paste(LOW_FLAG, collapse = ", ")))
if (length(intersect(LOW_FLAG, SETS$full$ANCHOR$human)) == 0)
  cat("    NOTE: no anchor gene is flagged 'low' -> the ANCHOR score is IDENTICAL\n",
      "         between full and strict by construction. It is not independent\n",
      "         confirmation. strict only perturbs TPH_HUMAN and TPH_MOUSE.\n")

# ---------------------------------------------------------------------------
# 3. Score both objects
# ---------------------------------------------------------------------------
# Set names contain underscores (TPH_HUMAN), so never recover set/variant by
# splitting a column name. Give each combination an opaque key and join a lookup.
KEY_MAP <- expand.grid(set = names(set_rows), variant = names(SETS),
                       stringsAsFactors = FALSE)[, c("variant", "set")]
KEY_MAP <- KEY_MAP[order(match(KEY_MAP$variant, names(SETS)),
                         match(KEY_MAP$set, names(set_rows))), ]
KEY_MAP$key <- paste0("S", seq_len(nrow(KEY_MAP)))
ANCHOR_FULL_KEY <- KEY_MAP$key[KEY_MAP$set == "ANCHOR" & KEY_MAP$variant == "full"]

score_object <- function(obj, species, group_col, target) {
  cat(sprintf("\n[3] Scoring %s (%d cells)\n", species, ncol(obj)))
  expressed <- rownames(obj)[Matrix::rowSums(GetAssayData(obj, layer = "counts") > 0) > 0]
  cov <- list(); feats <- list()

  for (i in seq_len(nrow(KEY_MAP))) {
    v <- KEY_MAP$variant[i]; s <- KEY_MAP$set[i]
    req <- SETS[[v]][[s]][[species]]
    got <- intersect(req, expressed)
    cov[[length(cov) + 1]] <- data.frame(
      species = species, variant = v, set = s, key = KEY_MAP$key[i],
      n_requested = length(req), n_scored = length(got),
      dropped = paste(setdiff(req, got), collapse = ";"))
    feats[[KEY_MAP$key[i]]] <- got
  }
  coverage <- bind_rows(cov)

  # AddModuleScore: mean of the set minus mean of ctrl genes drawn from the same
  # expression bins, so highly expressed genes cannot dominate.
  obj <- AddModuleScore(obj, features = feats, name = "MS__",
                        ctrl = CTRL_N, nbin = NBIN, seed = 1234)
  ms <- obj@meta.data[, grep("^MS__", colnames(obj@meta.data)), drop = FALSE]
  colnames(ms) <- paste0("AddModuleScore.", names(feats))

  # UCell: rank-based, insensitive to sequencing depth -- the sensitivity check.
  obj <- AddModuleScore_UCell(obj, features = feats, name = "")
  uc <- obj@meta.data[, names(feats), drop = FALSE]
  colnames(uc) <- paste0("UCell.", names(feats))

  scores <- cbind(
    data.frame(cell = colnames(obj),
               cluster = as.character(obj@meta.data[[group_col]]),
               condition = as.character(obj$condition),
               donor = if ("donor" %in% colnames(obj@meta.data))
                         as.character(obj$donor) else as.character(obj$orig.ident),
               stringsAsFactors = FALSE),
    ms, uc)
  list(scores = scores, coverage = coverage, obj = obj)
}

human <- readRDS(HUMAN_RDS)
stopifnot(sum(human$cluster_code == HUMAN_TARGET) == 873)
H <- score_object(human, "human", "cluster_code", HUMAN_TARGET)
human <- H$obj

mouse <- readRDS(MOUSE_RDS)
stopifnot(sum(mouse$celltype_merged == MOUSE_TARGET) == 1045, ncol(mouse) == 9224)
M <- score_object(mouse, "mouse", "celltype_merged", MOUSE_TARGET)
mouse <- M$obj

coverage <- bind_rows(H$coverage, M$coverage)
write.csv(coverage, file.path(OBJ, "gene_set_coverage.csv"), row.names = FALSE)
cat("\n    gene set coverage:\n")
print(as.data.frame(coverage %>% select(-dropped, -key)), row.names = FALSE)
thin <- coverage %>% filter(n_scored < 5)
if (nrow(thin)) { cat("    WARNING — sets thinned below 5 genes:\n"); print(as.data.frame(thin)) }

write.csv(H$scores, gzfile(file.path(OBJ, "module_scores_percell_human.csv.gz")), row.names = FALSE)
write.csv(M$scores, gzfile(file.path(OBJ, "module_scores_percell_mouse.csv.gz")), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Statistics
# ---------------------------------------------------------------------------
# Cliff's delta from the Mann-Whitney U, via ranks: delta = 2*AUC - 1.
# Ranges -1..1; |delta| < 0.147 negligible, < 0.33 small, < 0.474 medium.
cliffs_delta <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(x, y))
  auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
  2 * auc - 1
}

long_scores <- function(sc, species) {
  sc %>%
    pivot_longer(cols = matches("^(AddModuleScore|UCell)\\."),
                 names_to = "col", values_to = "score") %>%
    # exactly one "." separates method from key; keys are S1..Sn with no dots
    separate(col, into = c("method", "key"), sep = "\\.", extra = "merge") %>%
    left_join(KEY_MAP, by = "key") %>%
    mutate(species = species)
}
L <- bind_rows(long_scores(H$scores, "human"), long_scores(M$scores, "mouse"))

cat("\n[4] Kruskal-Wallis + Dunn across clusters\n")
kw <- list(); dn <- list(); tvr <- list()
combos <- L %>% distinct(species, method, set, variant)

for (i in seq_len(nrow(combos))) {
  cb <- combos[i, ]
  d  <- L %>% filter(species == cb$species, method == cb$method,
                     set == cb$set, variant == cb$variant) %>%
        mutate(cluster = factor(cluster))
  target <- if (cb$species == "human") HUMAN_TARGET else MOUSE_TARGET

  k <- kruskal.test(score ~ cluster, data = d)
  kw[[i]] <- data.frame(cb, statistic = unname(k$statistic), df = unname(k$parameter),
                        p_value = k$p.value, n_cells = nrow(d),
                        n_clusters = nlevels(d$cluster))

  # Dunn twice: BH (this module) and Bonferroni (matches Module 08 convention)
  dbh <- dunnTest(score ~ cluster, data = d, method = "bh")$res
  dbf <- dunnTest(score ~ cluster, data = d, method = "bonferroni")$res
  dd  <- dbh %>% select(Comparison, Z, P.unadj, P.adj.BH = P.adj) %>%
         left_join(dbf %>% select(Comparison, P.adj.bonferroni = P.adj), by = "Comparison")
  dn[[i]] <- data.frame(cb, dd)

  # Target-vs-each-other-cluster, with direction and effect size. Dunn's Z sign
  # follows the alphabetical order of the two group labels in `Comparison`, so
  # recover direction from the medians rather than trusting the sign.
  tv <- lapply(setdiff(levels(d$cluster), target), function(other) {
    a <- d$score[d$cluster == target]; b <- d$score[d$cluster == other]
    row <- dd[dd$Comparison %in% c(paste(target, other, sep = " - "),
                                   paste(other, target, sep = " - ")), ][1, ]
    data.frame(cb, other_cluster = other,
               median_target = median(a), median_other = median(b),
               median_diff = median(a) - median(b),
               cliffs_delta = cliffs_delta(a, b),
               Z = row$Z, P.unadj = row$P.unadj,
               P.adj.BH = row$P.adj.BH, P.adj.bonferroni = row$P.adj.bonferroni)
  }) %>% bind_rows()
  tvr[[i]] <- tv
}

kruskal_results <- bind_rows(kw)
dunn_results    <- bind_rows(dn)
target_vs_rest  <- bind_rows(tvr) %>%
  mutate(higher_in_target = median_diff > 0,
         significant_BH   = P.adj.BH < 0.05,
         sig_and_higher   = higher_in_target & significant_BH)

write.csv(kruskal_results, file.path(OBJ, "kruskal_results.csv"), row.names = FALSE)
write.csv(dunn_results,    file.path(OBJ, "dunn_results.csv"),    row.names = FALSE)
write.csv(target_vs_rest,  file.path(OBJ, "target_vs_rest_summary.csv"), row.names = FALSE)

# --- headline table ----------------------------------------------------------
headline <- target_vs_rest %>%
  group_by(species, method, set, variant) %>%
  summarise(n_contrasts = n(), n_sig_and_higher = sum(sig_and_higher),
            median_cliffs_delta = round(median(cliffs_delta), 3),
            max_cliffs_delta = round(max(cliffs_delta), 3), .groups = "drop") %>%
  mutate(verdict = paste0(n_sig_and_higher, "/", n_contrasts))
cat("\n    target-vs-rest summary (target = C9 / Tfh-Tph-like):\n")
print(as.data.frame(headline %>% filter(method == "AddModuleScore")), row.names = FALSE)
cat("\n    (UCell sensitivity check:)\n")
print(as.data.frame(headline %>% filter(method == "UCell")), row.names = FALSE)
write.csv(headline, file.path(OBJ, "target_vs_rest_headline.csv"), row.names = FALSE)

# --- positive-control gate ---------------------------------------------------
cat("\n    POSITIVE-CONTROL GATE (within-species set must peak at its own target)\n")
peaks <- L %>% filter(method == "AddModuleScore", variant == "full") %>%
  group_by(species, set, cluster) %>% summarise(m = mean(score), .groups = "drop") %>%
  group_by(species, set) %>% slice_max(m, n = 1) %>% ungroup()
h_ctrl <- peaks %>% filter(species == "human", set == "TPH_HUMAN") %>% pull(cluster)
m_ctrl <- peaks %>% filter(species == "mouse", set == "TPH_MOUSE") %>% pull(cluster)
cat(sprintf("      TPH_HUMAN peaks in human at: %s  (expect %s) -> %s\n",
            h_ctrl, HUMAN_TARGET, ifelse(h_ctrl == HUMAN_TARGET, "PASS", "FAIL")))
cat(sprintf("      TPH_MOUSE peaks in mouse at: %s  (expect %s) -> %s\n",
            m_ctrl, MOUSE_TARGET, ifelse(m_ctrl == MOUSE_TARGET, "PASS", "FAIL")))
cat("\n    where every set peaks:\n"); print(as.data.frame(peaks), row.names = FALSE)
write.csv(peaks, file.path(OBJ, "score_peak_cluster.csv"), row.names = FALSE)

# --- method concordance ------------------------------------------------------
conc <- L %>%
  pivot_wider(id_cols = c(species, cell, set, variant),
              names_from = method, values_from = score) %>%
  group_by(species, set, variant) %>%
  summarise(spearman_rho = round(cor(AddModuleScore, UCell, method = "spearman"), 3),
            .groups = "drop")
write.csv(conc, file.path(OBJ, "score_method_concordance.csv"), row.names = FALSE)
cat("\n    AddModuleScore vs UCell concordance (Spearman rho):\n")
print(as.data.frame(conc), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Pseudo-replication checks
#
# Cells within a donor are not independent. Human has 18 donors (PD164a/b/c and
# PD161a/b are repeat samples from one patient), so collapse to donor means.
# Mouse has 8 samples, 2 per condition -- descriptive only, no test.
# ---------------------------------------------------------------------------
cat("\n[5] Pseudo-replication checks\n")
h_long <- L %>% filter(species == "human")
c9_per_donor <- H$scores %>% filter(cluster == HUMAN_TARGET) %>% count(donor, name = "n_C9")
usable <- c9_per_donor$donor[c9_per_donor$n_C9 >= MIN_C9_CELLS_PER_DONOR]
cat(sprintf("    human donors with >= %d C9 cells: %d of %d\n",
            MIN_C9_CELLS_PER_DONOR, length(usable), nrow(c9_per_donor)))

donor_res <- h_long %>%
  filter(donor %in% usable) %>%
  mutate(grp = ifelse(cluster == HUMAN_TARGET, "target", "rest")) %>%
  group_by(method, set, variant, donor, grp) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  pivot_wider(names_from = grp, values_from = mean_score) %>%
  group_by(method, set, variant) %>%
  summarise(n_donors = n(),
            median_target = median(target), median_rest = median(rest),
            median_paired_diff = median(target - rest),
            p_paired_wilcoxon = wilcox.test(target, rest, paired = TRUE)$p.value,
            .groups = "drop") %>%
  mutate(p_paired_wilcoxon = signif(p_paired_wilcoxon, 3))
write.csv(donor_res, file.path(OBJ, "donor_level_human.csv"), row.names = FALSE)
cat("    donor-level paired Wilcoxon, C9 vs rest (AddModuleScore):\n")
print(as.data.frame(donor_res %>% filter(method == "AddModuleScore") %>%
                    mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)

mouse_sample_means <- L %>% filter(species == "mouse") %>%
  mutate(grp = ifelse(cluster == MOUSE_TARGET, "target", "rest")) %>%
  group_by(method, set, variant, sample = donor, condition, grp) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
write.csv(mouse_sample_means, file.path(OBJ, "sample_level_means_mouse.csv"), row.names = FALSE)
cat("\n    mouse: 8 samples, 2 per condition -> sample-level means written\n")
cat("    descriptively; NO test run (n = 2 per group is not testable).\n")

# ---------------------------------------------------------------------------
# 6. Chronic-activation score by condition, within the Tph cluster
# ---------------------------------------------------------------------------
cat("\n[6] Chronic-activation score by condition, within the Tph cluster\n")
cond_rows <- list()

for (v in c("full", "strict")) {
  # --- human C9: Control vs PD, three ways ---------------------------------
  d <- h_long %>% filter(set == "CHRONIC", variant == v, method == "AddModuleScore",
                         cluster == HUMAN_TARGET)
  a <- d$score[d$condition == "Control"]; b <- d$score[d$condition == "PD"]
  cond_rows[[length(cond_rows) + 1]] <- data.frame(
    species = "human", variant = v, level = "cell-level (all C9)",
    group_a = sprintf("Control n=%d", length(a)), group_b = sprintf("PD n=%d", length(b)),
    statistic = "Wilcoxon", effect = round(cliffs_delta(b, a), 3),
    p_value = signif(wilcox.test(b, a)$p.value, 3),
    note = "cells within a donor are not independent")

  dm <- d %>% filter(donor %in% usable) %>%
        group_by(donor, condition) %>% summarise(m = mean(score), .groups = "drop")
  ca <- dm$m[dm$condition == "Control"]; cb <- dm$m[dm$condition == "PD"]
  cond_rows[[length(cond_rows) + 1]] <- data.frame(
    species = "human", variant = v, level = "donor-mean",
    group_a = sprintf("Control n=%d", length(ca)), group_b = sprintf("PD n=%d", length(cb)),
    statistic = "Wilcoxon", effect = round(cliffs_delta(cb, ca), 3),
    p_value = signif(wilcox.test(cb, ca)$p.value, 3),
    note = sprintf("donors with >= %d C9 cells", MIN_C9_CELLS_PER_DONOR))

  # PD164 alone contributes 343 of 616 PD C9 cells (56%) -- leave it out.
  dl <- dm %>% filter(donor != "PD164")
  la <- dl$m[dl$condition == "Control"]; lb <- dl$m[dl$condition == "PD"]
  cond_rows[[length(cond_rows) + 1]] <- data.frame(
    species = "human", variant = v, level = "donor-mean, PD164 dropped",
    group_a = sprintf("Control n=%d", length(la)), group_b = sprintf("PD n=%d", length(lb)),
    statistic = "Wilcoxon", effect = round(cliffs_delta(lb, la), 3),
    p_value = if (length(lb) > 1 && length(la) > 1)
                signif(wilcox.test(lb, la)$p.value, 3) else NA_real_,
    note = "PD164 supplies 56% of PD C9 cells")

  # --- mouse Tfh/Tph-like: BAS/PD/REC/RCU ----------------------------------
  dmz <- L %>% filter(species == "mouse", set == "CHRONIC", variant == v,
                      method == "AddModuleScore", cluster == MOUSE_TARGET) %>%
         mutate(condition = factor(condition))
  kk <- kruskal.test(score ~ condition, data = dmz)
  cond_rows[[length(cond_rows) + 1]] <- data.frame(
    species = "mouse", variant = v, level = "cell-level (Tfh/Tph-like)",
    group_a = paste(levels(dmz$condition), collapse = "/"),
    group_b = paste(as.integer(table(dmz$condition)), collapse = "/"),
    statistic = "Kruskal-Wallis", effect = round(unname(kk$statistic), 2),
    p_value = signif(kk$p.value, 3),
    note = "only 8 samples (2/condition); cells not independent")
  dnm <- dunnTest(score ~ condition, data = dmz, method = "bh")$res %>%
    mutate(species = "mouse", variant = v, set = "CHRONIC")
  write.csv(dnm, file.path(OBJ, sprintf("chronic_condition_dunn_mouse_%s.csv", v)),
            row.names = FALSE)
}
chronic_cond <- bind_rows(cond_rows)
write.csv(chronic_cond, file.path(OBJ, "chronic_by_condition_in_tph.csv"), row.names = FALSE)
print(as.data.frame(chronic_cond), row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Figures — theme_classic, no grid (repo convention)
# ---------------------------------------------------------------------------
cat("\n[7] Figures\n")
no_grid <- theme_classic(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

violin_panel <- function(species, target, accent, file, w, h) {
  d <- L %>% filter(species == !!species, method == "AddModuleScore") %>%
       mutate(is_target = cluster == target,
              set = factor(set, levels = c("ANCHOR","TPH_HUMAN","TPH_MOUSE","CHRONIC")))
  ord <- d %>% filter(set == "ANCHOR", variant == "full") %>%
         group_by(cluster) %>% summarise(m = median(score)) %>% arrange(desc(m)) %>% pull(cluster)
  d$cluster <- factor(d$cluster, levels = ord)
  p <- ggplot(d, aes(cluster, score, fill = is_target)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.2, fill = "white") +
    scale_fill_manual(values = c(`TRUE` = accent, `FALSE` = "grey80"), guide = "none") +
    facet_grid(set ~ variant, scales = "free_y") +
    labs(title = paste0(species, " — AddModuleScore by cluster"),
         subtitle = paste0("target = ", target,
                           " (accent). ANCHOR is identical between full and strict by construction."),
         x = NULL, y = "module score") +
    no_grid
  svglite(file.path(FIG, file), width = w, height = h); print(p); invisible(dev.off())
  cat("    saved:", file, "\n")
}
violin_panel("human", HUMAN_TARGET, HUMAN_ACCENT, "violin_scores_by_cluster_human.svg", 8, 10)
violin_panel("mouse", MOUSE_TARGET, MOUSE_ACCENT, "violin_scores_by_cluster_mouse.svg", 9, 10)

# --- anchor score on each species' UMAP -------------------------------------
# AddModuleScore wrote MS__1..MS__n in the order feats were supplied (= KEY_MAP
# row order), so resolve the ANCHOR/full column by position rather than assuming.
ANCHOR_MS_COL <- paste0("MS__", which(KEY_MAP$key == ANCHOR_FULL_KEY))
stopifnot(ANCHOR_MS_COL %in% colnames(human@meta.data),
          ANCHOR_MS_COL %in% colnames(mouse@meta.data))
fp <- (FeaturePlot(human, features = ANCHOR_MS_COL, reduction = "umap.tcell", order = TRUE) +
         scale_colour_viridis_c(option = "magma") +
         ggtitle("human — ANCHOR score (umap.tcell)") + theme_classic(base_size = 9) +
         theme(panel.grid = element_blank())) |
      (FeaturePlot(mouse, features = ANCHOR_MS_COL, reduction = "umap.pca_TCRab", order = TRUE) +
         scale_colour_viridis_c(option = "magma") +
         ggtitle("mouse — ANCHOR score (umap.pca_TCRab)") + theme_classic(base_size = 9) +
         theme(panel.grid = element_blank()))
svglite(file.path(FIG, "featureplot_anchor_score.svg"), width = 11, height = 5)
print(fp); invisible(dev.off())
cat("    saved: featureplot_anchor_score.svg\n")

# --- mean score heatmap, z-scored within set --------------------------------
hm_df <- L %>% filter(method == "AddModuleScore", variant == "full") %>%
  group_by(species, set, cluster) %>% summarise(m = mean(score), .groups = "drop") %>%
  group_by(species, set) %>% mutate(z = as.numeric(scale(m))) %>% ungroup() %>%
  mutate(set = factor(set, levels = c("ANCHOR","TPH_HUMAN","TPH_MOUSE","CHRONIC")),
         is_target = (species == "human" & cluster == HUMAN_TARGET) |
                     (species == "mouse" & cluster == MOUSE_TARGET))
ph <- ggplot(hm_df, aes(cluster, set, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = subset(hm_df, is_target), colour = "black", linewidth = 0.8, fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "z (within set)") +
  # space = "free_x" so each panel gets width proportional to its cluster count;
  # facet_wrap forces equal panel widths and clips the 10-cluster mouse panel.
  facet_grid(~ species, scales = "free_x", space = "free_x") +
  labs(title = "Mean module score by cluster, z-scored within each gene set",
       subtitle = "Black outline = the Tph cluster of that species", x = NULL, y = NULL) +
  no_grid
svglite(file.path(FIG, "heatmap_mean_score_by_cluster.svg"), width = 13, height = 4.6)
print(ph); invisible(dev.off())
cat("    saved: heatmap_mean_score_by_cluster.svg\n")

# --- anchor genes, gene by gene ---------------------------------------------
anchors_h <- SETS$full$ANCHOR$human; anchors_m <- SETS$full$ANCHOR$mouse
Idents(human) <- "cluster_code"; Idents(mouse) <- "celltype_merged"
pd <- (DotPlot(human, features = anchors_h) + coord_flip() +
         ggtitle("human — 13 anchor genes") + no_grid) |
      (DotPlot(mouse, features = anchors_m) + coord_flip() +
         ggtitle("mouse — 13 anchor genes") + no_grid)
svglite(file.path(FIG, "dotplot_anchor_genes.svg"), width = 13, height = 5.5)
print(pd); invisible(dev.off())
cat("    saved: dotplot_anchor_genes.svg\n")

# --- full vs strict, for the two sets that actually change ------------------
fs <- L %>% filter(method == "AddModuleScore", set %in% c("TPH_HUMAN", "TPH_MOUSE")) %>%
  pivot_wider(id_cols = c(species, cell, set, cluster), names_from = variant,
              values_from = score) %>%
  mutate(is_target = (species == "human" & cluster == HUMAN_TARGET) |
                     (species == "mouse" & cluster == MOUSE_TARGET))
pfs <- ggplot(fs, aes(full, strict, colour = is_target)) +
  geom_point(size = 0.25, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70"),
                      labels = c("other clusters", "Tph cluster"), name = NULL) +
  facet_grid(species ~ set) +
  labs(title = "full vs strict gene sets (only TPH_HUMAN and TPH_MOUSE change)",
       subtitle = "ANCHOR and CHRONIC are unaffected — no 'low'-flagged genes in them") +
  no_grid + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
svglite(file.path(FIG, "full_vs_strict_scatter.svg"), width = 8, height = 6)
print(pfs); invisible(dev.off())
cat("    saved: full_vs_strict_scatter.svg\n")

# --- chronic score by condition within the Tph cluster ----------------------
cc <- L %>% filter(set == "CHRONIC", method == "AddModuleScore", variant == "full") %>%
  filter((species == "human" & cluster == HUMAN_TARGET) |
         (species == "mouse" & cluster == MOUSE_TARGET))
pcc <- ggplot(cc, aes(condition, score, fill = species)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.2, fill = "white") +
  scale_fill_manual(values = c(human = HUMAN_ACCENT, mouse = MOUSE_ACCENT), guide = "none") +
  facet_wrap(~species, scales = "free_x") +
  labs(title = "Chronic-activation score within the Tph cluster",
       subtitle = "human C9 (Control vs PD) | mouse Tfh/Tph-like (BAS/PD/REC/RCU)",
       x = NULL, y = "chronic-activation module score") +
  no_grid
svglite(file.path(FIG, "chronic_by_condition_in_tph.svg"), width = 7, height = 4)
print(pcc); invisible(dev.off())
cat("    saved: chronic_by_condition_in_tph.svg\n")

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_02.txt"))
cat("\nDone — script 02.\n")
