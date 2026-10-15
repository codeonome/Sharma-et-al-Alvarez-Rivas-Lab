# =====================================================================
# 07 — T-cell level: cleanup, annotation, differential abundance, figures
# VERSION 4 — consolidated. Supersedes v3 and the loose 10e-bis / 11d files.
#
# Input : T-cell object clustered at res 0.4 (11 clusters, 21,899 cells)
# Output: cleaned object (9 clusters), UMAPs, composition CSVs,
#         differential abundance results, figures
#
# Run strictly top-to-bottom. Nothing after section 8 reads the res-0.4
# cluster column; everything reads `cluster_id` / `cluster_code`.
#
# CHANGES FROM v3
#   - 10e: topTable() now called with confint = TRUE, and ci_low / ci_high
#          are carried into the results tibble (needed by 11d)
#   - 10h: NEW. Variance components + rho sensitivity sweep.
#          Fixes the as.data.frame(table) bug: getTransformedProps returns
#          a `table`, which as.data.frame melts to long format instead of
#          keeping the matrix shape. Fixed with as.matrix(unclass(.)) and
#          check.names = FALSE.
#   - 11d: REWRITTEN. Bar plot now shows the log-odds effect from the
#          blocked model with 95% CIs and NUMERIC FDR labels, replacing
#          the old proportion-scale log2FC bars with asterisks.
#
# PRIMARY REPORTED MODEL
#   propeller on logit proportions, sample level, donor as blocking
#   factor, rho from duplicateCorrelation. Section 10e.
#
# SECTION MAP
#   1  Libraries
#   2  Config            <- edit here, nowhere else
#   3  Helper functions
#   4  Load + validate
#   5  Annotation reference (record only)
#   6  Fixes: Ig doublets, cluster removal, cluster merge, marker filter
#   7  Changelog
#   8  Final labels, display codes, palette
#   9  Composition tables (sample + donor)
#  10  DIFFERENTIAL ABUNDANCE
#        10a Wilcoxon + BH, sample level
#        10b Wilcoxon + BH, donor level
#        10c propeller, sample level
#        10d propeller, donor level
#        10e propeller, donor-blocked        <- PRIMARY
#        10f combined results table
#        10g leave-one-donor-out sensitivity
#        10h variance components + rho sensitivity sweep
#  11  Figures
#  12  Save + manifest + sessionInfo
# =====================================================================


# =====================================================================
# 1. LIBRARIES
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

set.seed(42)


# =====================================================================
# 2. CONFIG — the only block you should need to edit
# =====================================================================
# Upstream Harmony-clustered human T-cell object (not part of this repository;
# see Human/README.md). Place it next to this script before running.
IN_RDS    <- "tcell_seurat_object_harmony_after_clustering.rds"
RES_COL   <- "tcell_res.0.4"
REDUCTION <- "umap.tcell"

GM_LABEL  <- "Control"      # value in meta.data$condition for healthy
PD_LABEL  <- "PD"           # value in meta.data$condition for periodontitis

FIG_DIR <- file.path("qc_plots3", "final")
RES_DIR <- "results3"

MARKER_CSVS <- file.path(RES_DIR, "tcell_subset_markers.csv")

EXPECTED_CLUSTERS <- as.character(0:10)
DROP_CLUSTERS     <- "9"                     # low-complexity ribosomal artifact
MERGE_MAP         <- c("10" = "0")           # names = donor cluster, values = target

DO_IG_DOUBLET_REMOVAL <- TRUE
IG_DOUBLET_CUTOFF     <- 10                  # percent of counts from Ig genes
DO_MARKER_FILTER      <- TRUE

FDR_CUTOFF       <- 0.05
PROP_TRANSFORM   <- "logit"                  # "logit" or "asin"
SAMPLE_SUFFIX_RE <- "[a-z]+$"                # PD164b -> PD164 (repeat biopsy)

# rho values for the sensitivity sweep (section 10h). 0.99 is retained
# deliberately: it is a numerical boundary, not a plausible correlation,
# and deleting the one grid point that does not support the conclusion
# would be selective reporting. It is explained in the caption instead.
RHO_GRID_BASE <- c(0, 0.25, 0.5, 0.85, 0.95, 0.99)

# Figure annotation source. The blocked model is primary; the chain is a
# fallback so the script still finishes if speckle/limma are missing.
PRIMARY_TEST  <- "propeller_blocked_sample"
TEST_FALLBACK <- c("propeller_blocked_sample", "propeller_sample",
                   "wilcox_sample", "wilcox_donor")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)


# =====================================================================
# 3. HELPER FUNCTIONS
# =====================================================================
# Anything done more than once lives here. Writing the test ONCE and
# calling it twice (sample level, donor level) is the point of a function:
# change the effect-size definition in one place and both levels stay
# consistent. The single most useful habit to build.

# --- 3a. Wilcoxon rank-sum on a long proportion table ----------------
# `prop_df` needs columns: cluster_id, condition, proportion.
# One row per unit (sample or donor) per cluster.
wilcox_da <- function(prop_df, level_label) {
  # eps: half the smallest non-zero proportion anywhere, so a cluster
  # absent from one condition gives a finite fold change instead of Inf,
  # and the offset never dominates a real value.
  eps <- min(prop_df$proportion[prop_df$proportion > 0], na.rm = TRUE) / 2

  prop_df %>%
    group_by(cluster_id) %>%
    summarise(
      n_GM      = sum(condition == GM_LABEL),
      n_PD      = sum(condition == PD_LABEL),
      mean_GM   = mean(proportion[condition == GM_LABEL]),
      mean_PD   = mean(proportion[condition == PD_LABEL]),
      median_GM = median(proportion[condition == GM_LABEL]),
      median_PD = median(proportion[condition == PD_LABEL]),
      # Called twice (statistic, then p) for readability. Tie warnings are
      # expected and do not invalidate the p-value.
      statistic = tryCatch(
        unname(wilcox.test(proportion[condition == PD_LABEL],
                           proportion[condition == GM_LABEL])$statistic),
        error = function(e) NA_real_),
      p_value   = tryCatch(
        wilcox.test(proportion[condition == PD_LABEL],
                    proportion[condition == GM_LABEL])$p.value,
        error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      effect_type = "log2FC_proportion",
      effect      = log2((mean_PD + eps) / (mean_GM + eps)),
      ci_low      = NA_real_,
      ci_high     = NA_real_,
      FDR         = p.adjust(p_value, method = "BH"),
      sig         = ifelse(!is.na(FDR) & FDR < FDR_CUTOFF, "*", ""),
      method      = "Wilcoxon_BH",
      level       = level_label,
      cluster_id  = as.character(cluster_id)
    ) %>%
    arrange(p_value)
}

# --- 3b. propeller on cell-level vectors -----------------------------
# propeller logit-transforms the per-unit proportions to stabilise
# variance, then fits a limma model with EMPIRICAL BAYES moderation: each
# cluster's variance estimate is shrunk toward the average across all 9
# clusters. That shrinkage is where the power gain over Wilcoxon comes
# from. Note this is a SEPARATE mechanism from the rho blocking in 10e.
propeller_da <- function(clusters, unit, group, level_label) {
  res <- speckle::propeller(clusters  = clusters,
                            sample    = unit,
                            group     = group,
                            transform = PROP_TRANSFORM)

  # propeller returns clusters as ROWNAMES (usually also as a
  # BaselineProp.clusters column). Tibbles drop rownames silently, so pull
  # them out before any dplyr touches the object.
  cl <- if ("BaselineProp.clusters" %in% names(res)) {
    as.character(res$BaselineProp.clusters)
  } else {
    rownames(res)
  }

  # PropMean columns are named PropMean.<group level>; find by pattern
  # rather than hard-coding "PropMean.Control".
  pm     <- grep("^PropMean", names(res), value = TRUE)
  gm_col <- pm[grepl(GM_LABEL, pm, fixed = TRUE)]
  pd_col <- pm[grepl(PD_LABEL, pm, fixed = TRUE)]
  stopifnot(length(gm_col) == 1, length(pd_col) == 1)

  n_units <- tapply(as.character(group), as.character(unit), function(x) x[1])

  tibble(
    cluster_id  = cl,
    n_GM        = sum(n_units == GM_LABEL),
    n_PD        = sum(n_units == PD_LABEL),
    mean_GM     = res[[gm_col]],
    mean_PD     = res[[pd_col]],
    median_GM   = NA_real_,
    median_PD   = NA_real_,
    statistic   = if ("Tstatistic" %in% names(res)) res$Tstatistic else NA_real_,
    p_value     = res$P.Value,
    FDR         = res$FDR,
    effect_type = "log2_propRatio",
    effect      = log2(res[[pd_col]] / res[[gm_col]]),
    ci_low      = NA_real_,
    ci_high     = NA_real_,
    method      = paste0("propeller_", PROP_TRANSFORM),
    level       = level_label
  ) %>%
    mutate(sig = ifelse(!is.na(FDR) & FDR < FDR_CUTOFF, "*", "")) %>%
    arrange(p_value)
}

# --- 3c. Seurat-version-safe expression fetch ------------------------
fetch_expr <- function(obj, vars) {
  tryCatch(FetchData(obj, vars = vars, layer = "data"),
           error = function(e) FetchData(obj, vars = vars, slot = "data"))
}

# --- 3d. getTransformedProps -> samples x clusters data frame --------
# THE v3 BUG, FIXED HERE.
# getTransformedProps() returns `TransformedProps` as a `table`, not a
# plain matrix. as.data.frame() is GENERIC: as.data.frame.table melts a
# table into long format (clusters, sample, Freq) rather than preserving
# the 9 x 21 shape. Downstream, tp[[cluster_name]] then returned NULL and
# tp$donor <- smeta$donor silently RECYCLED 21 values into 189 rows,
# because 189 is divisible by 21. Silent recycling is the dangerous part.
#
# unclass() strips the table class so as.matrix() keeps the shape;
# check.names = FALSE stops R mangling "CD8 TEM (GZMK/GZMH)" into
# "CD8.TEM..GZMK.GZMH.". Both are needed.
props_to_df <- function(tprops, smeta) {
  m  <- as.matrix(unclass(tprops))                       # clusters x samples
  df <- as.data.frame(t(m), check.names = FALSE)         # samples x clusters
  df$donor     <- smeta$donor
  df$condition <- as.character(smeta$condition)
  # Guard rail: catches the recycling failure mode immediately.
  stopifnot(nrow(df) == ncol(m),
            all(rownames(m) %in% colnames(df)))
  df
}

# --- 3e. Column order for the combined results table -----------------
DA_COLS <- c("test_id", "method", "level", "cluster_id", "display_code",
             "n_GM", "n_PD", "mean_GM", "mean_PD", "median_GM", "median_PD",
             "effect_type", "effect", "ci_low", "ci_high",
             "statistic", "p_value", "FDR", "sig")


# =====================================================================
# 4. LOAD AND VALIDATE
# =====================================================================
tcell <- readRDS(IN_RDS)

stopifnot(
  RES_COL   %in% colnames(tcell@meta.data),
  REDUCTION %in% Reductions(tcell),
  all(c("sample", "condition") %in% colnames(tcell@meta.data))
)

cat("Loaded:", ncol(tcell), "cells,",
    nlevels(factor(tcell@meta.data[[RES_COL]])), "clusters at", RES_COL, "\n")
cat("Condition values present:", paste(unique(tcell$condition), collapse = ", "), "\n")

if (!all(c(GM_LABEL, PD_LABEL) %in% unique(tcell$condition))) {
  stop("GM_LABEL / PD_LABEL do not match the `condition` column. Fix section 2.")
}

# Work off an explicit column, not the active identity. Idents() gets
# silently overwritten when several resolutions live in the metadata.
tcell$res04 <- factor(as.character(tcell@meta.data[[RES_COL]]),
                      levels = EXPECTED_CLUSTERS)
stopifnot(!any(is.na(tcell$res04)))

# --- Donor column ----------------------------------------------------
# Trailing lowercase letter marks a repeat biopsy: PD164b -> PD164.
# ALWAYS eyeball this mapping. One naming exception silently merges two
# real donors and nothing downstream would reveal it.
tcell$donor <- sub(SAMPLE_SUFFIX_RE, "", tcell$sample)

donor_key <- tcell@meta.data %>%
  count(donor, sample, condition, name = "n_cells") %>%
  arrange(condition, donor, sample)
write.csv(donor_key, file.path(RES_DIR, "sample_to_donor_map.csv"), row.names = FALSE)

cat("\n--- sample -> donor mapping (CHECK THIS) ---\n")
print(as.data.frame(donor_key))
cat("Samples:", n_distinct(tcell$sample), "| Donors:", n_distinct(tcell$donor), "\n")

n_start      <- ncol(tcell)
counts_start <- table(tcell$res04)


# =====================================================================
# 5. ANNOTATION REFERENCE — res 0.4, 11 clusters, 21,899 cells
#    (record only; the code that acts on it is section 6)
# =====================================================================
# UMAP layout: CD8 lobe = 0, 2, 6, 10 (left) | CD4 lobe = 1, 3, 4, 5, 8 (right)
#              MAIT (7) its own island | 9 dead-centre between the lobes
#
# Cluster 0  : CD8+ effector memory (GZMK/GZMH)          [HIGH]
#              CD8A, CD8B, GZMH, GZMK, NKG7, CCL4/CCL5, CST7, CRTAM, HLA-DR
#
# Cluster 1  : CD4+ naive / central memory               [HIGH]
#              IL7R, CCR7, TCF7, LTB, ANXA1, KLF2, SELL
#
# Cluster 2  : CD8+ effector memory (ZEB2-hi, late-diff) [MODERATE]
#              GZMK, ZEB2, AOAH, CRTAM, KLRK1, DTHD1
#              !! ISSUE: ambient Ig. 4 Ig genes in top markers at
#                 pct.1 ~ pct.2 (IGKC 0.648 vs 0.632) -> inflated logFC,
#                 not expression. Highest percent.ig of any cluster (1.645).
#              !! ISSUE (DISPROVEN): suspected NK contamination -- CD3D/E/TRAC
#                 ~half other clusters. KLRF1 <= 0.063 everywhere. Low CD3
#                 tracks library complexity, not lost lineage identity. RETAINED.
#
# Cluster 3  : Th17 / CCR6+ CD4 effector                 [HIGH, resolved]
#              IL17A 19.7% (vs 2.4% next), CCR6 53.7%, IL1R1 37.8%,
#              IL23R 16.5%, CCL20, RORA, KLRB1, BATF, MAF
#              !! CORRECTION: first read as Th17-like effector Treg because
#                 TNFRSF18/GITR, TNFRSF4/OX40, IL2RA/CD25 were in its markers.
#                 Those are ACTIVATION markers, not Treg-lineage markers.
#                 FOXP3+RORC+ double-positives = 17 cells vs ~8 by chance:
#                 real direction, too few to call a population.
#
# Cluster 4  : Regulatory T cells (Treg)                 [HIGH]
#              IKZF2 51.4% vs 11.6%, FOXP3 22.6% vs 10.7%, RTKN2,
#              IL2RA, CTLA4, TIGIT, TNFRSF4/18
#
# Cluster 5  : CD4+ resting / central memory (TSHZ2+)    [MODERATE]
#              TSHZ2, LDLRAD4, RORA, SERINC5, PHTF2, BCL2
#              !! ISSUE (DISPROVEN): same NK suspicion as cluster 2. RETAINED.
#                 Second-highest percent.ig (1.304) -- again complexity-driven.
#
# Cluster 6  : CD8+ cytotoxic effector (GZMB/GNLY, TEMRA-like) [HIGH]
#              GNLY, GZMB, HOPX, CTSW, CD7, CD8A/B, NELL2, ID2
#
# Cluster 7  : MAIT cells                                [HIGH]
#              SLC4A10 log2FC 5.77 (38.1% vs 0.8%), KLRB1 in 90%,
#              ZBTB16, NCR3, IL23R, IL4I1
#              -> NOT contamination. MR1-restricted + semi-invariant but
#                 still alpha-beta (TRAV1-2+). Retained, labelled separately.
#                 Exclude only from analyses restricted to conventional T cells.
#
# Cluster 8  : CXCL13+ Tph-like CD4 (peripheral helper)  [HIGH identity]
#              CXCL13 32.3%, TOX2 46.1%, PTPN13 55.1%, TOX, PDE7B, ICOS, MAF
#              !! CAVEAT: PDCD1 only 9.1% -- below canonical Tph definition.
#                 Tph vs Tfh still needs CXCR5 (neg = Tph). Keep the caveat in text.
#
# Cluster 9  : LOW-COMPLEXITY / RIBOSOMAL ARTIFACT       [NOT A CELL TYPE]
#              all 118 DEGs are RPL/RPS, max log2FC 0.82, no lineage marker
#              nFeature 570 vs 989 | percent.ribo 24.2% vs 16.7%
#              percent.mt 4.06% vs 4.44% -> NOT elevated, so shallow sampling,
#              not cell death. In 21/21 samples, max 16% from one sample,
#              condition-balanced (462 GM / 402 PD) -> safe to drop.
#              -> ACTION: REMOVE (section 6c)
#
# Cluster 10 : Activated/transitional T                  [LOW CONFIDENCE]
#              GZMK, GZMA, DUSP4, CTLA4, ICOS, MAF, RGS1 -- only 46 sig genes
#              Sits inside cluster 0's UMAP territory; mixes CD8 cytotoxic
#              with CD4 activation programmes -> over-clustering at res 0.4
#              or a doublet pocket, not a population.
#              -> ACTION: MERGE INTO CLUSTER 0 (section 6d)
#
# Dataset-wide, not cluster-specific:
#   * Ambient Ig: IGKC detected in ~63% of T cells, but median percent.ig
#     0.17%, mean 0.88%. Real but tiny. -> filter Ig from MARKER TABLES only
#     (section 6f). No SoupX, no gene removal, no reclustering: at 0.17%
#     median it cannot move variable-feature selection or the embedding,
#     and rebuilding would reshuffle every cluster number.
#   * PD vs GM ambient Ig is NOT a condition effect: raw means 1.019 vs 0.391
#     but Wilcoxon on per-sample means p = 0.27, carried by one donor
#     (PD164b 3.808, PD164 1.936). Drop that donor -> PD 0.331 vs GM 0.391.
#     GM289 (1.872) shows controls reach high values too. DO NOT REPORT.
#   * T/plasma doublets: max percent.ig 68.6%, 1.2% of cells > 10%.
#     -> ACTION: REMOVE at 10% (section 6b)
#   * TCR-driven clustering: ruled out. Only constant regions appear as
#     markers (TRAC, TRBC1/2; logFC 0.39-0.73, ranks 82-335), ZERO
#     TRAV/TRBV/TRDV. Constant regions are shared by all alpha-beta T cells
#     and cannot create clonotype structure. Dropped as uninformative.
#   * gamma-delta / NK: already excluded upstream at the immune level
#     (immune cluster 7, TRDC+). No TRDC/TRGC signal in any subcluster here.


# =====================================================================
# 6. FIXES
# =====================================================================

# ---------------------------------------------------------------------
# 6a. Compute percent.ig (needed for the doublet filter)
# ---------------------------------------------------------------------
# grep() on rownames is the workhorse for gene-set selection.
# `value = TRUE` returns the gene NAMES rather than their positions --
# that is the difference between grep() and grepl()/which(). Always print
# the matched set: regex on gene symbols over-matches, hence the setdiff()
# removing IGHMBP2 and IGLON5.
ig_pattern <- "^IGH[ADEGM]|^IGHG[1-4]|^IGHGP|^IGK[CV]|^IGL[CV]"
ig_genes   <- grep(ig_pattern, rownames(tcell), value = TRUE)
ig_genes   <- setdiff(ig_genes, grep("^IGHMBP|^IGLON", ig_genes, value = TRUE))

cat("\nIg genes matched (", length(ig_genes), "):\n", sep = "")
print(ig_genes)

tcell <- PercentageFeatureSet(tcell, features = ig_genes, col.name = "percent.ig")

cat("\npercent.ig: median", round(median(tcell$percent.ig), 3),
    "| mean", round(mean(tcell$percent.ig), 3),
    "| max", round(max(tcell$percent.ig), 1), "\n")
cat("Cells >", IG_DOUBLET_CUTOFF, "% Ig:", sum(tcell$percent.ig > IG_DOUBLET_CUTOFF),
    sprintf("(%.2f%%)\n", 100 * mean(tcell$percent.ig > IG_DOUBLET_CUTOFF)))

# ---------------------------------------------------------------------
# 6b. Remove T-cell / plasma-cell doublets
# ---------------------------------------------------------------------
ig_high <- tcell$percent.ig > IG_DOUBLET_CUTOFF

# Record WHERE the doublets sat before deleting them. Once the cells are
# gone this table cannot be reconstructed, and it is what you cite if a
# reviewer asks whether the filter hit one cluster disproportionately.
ig_doublet_tbl <- tcell@meta.data[ig_high, ] %>%
  count(res04, condition, name = "n_doublets") %>%
  tidyr::complete(res04, condition, fill = list(n_doublets = 0))
write.csv(ig_doublet_tbl, file.path(RES_DIR, "ig_doublets_by_cluster.csv"),
          row.names = FALSE)

if (DO_IG_DOUBLET_REMOVAL) {
  doublet_cells <- colnames(tcell)[ig_high]
  cat("\nDoublet removal by condition:\n")
  print(table(tcell$condition[ig_high]))
  tcell <- subset(tcell, cells = setdiff(colnames(tcell), doublet_cells))
  cat("Removed", length(doublet_cells), "Ig-high doublets ->", ncol(tcell), "cells\n")
} else {
  doublet_cells <- character(0)
  cat("\nSkipping doublet removal (DO_IG_DOUBLET_REMOVAL = FALSE)\n")
}

# ---------------------------------------------------------------------
# 6c. Remove cluster 9 (low-complexity ribosomal artifact)
# ---------------------------------------------------------------------
n_dropped <- sum(tcell$res04 %in% DROP_CLUSTERS)
cat("\nCluster", DROP_CLUSTERS, "removal by condition:\n")
print(table(tcell$condition[tcell$res04 %in% DROP_CLUSTERS]))

tcell <- subset(tcell, cells = colnames(tcell)[!tcell$res04 %in% DROP_CLUSTERS])
cat("Removed", n_dropped, "cells ->", ncol(tcell), "cells\n")

# ---------------------------------------------------------------------
# 6d. Merge cluster 10 into cluster 0
# ---------------------------------------------------------------------
# Two things worth internalising.
# (1) Convert to character BEFORE editing. Factors store values as integer
#     codes pointing at a `levels` vector; assigning a value that is not an
#     existing level gives NA, silently. as.character() sidesteps it.
# (2) droplevels() removes the now-empty "9" and "10" levels. Without it,
#     empty levels persist as phantom rows in table(), phantom facets in
#     ggplot, and phantom legend entries.
res_chr  <- as.character(tcell$res04)
n_merged <- sum(res_chr %in% names(MERGE_MAP))

res_chr <- ifelse(res_chr %in% names(MERGE_MAP),
                  unname(MERGE_MAP[res_chr]),
                  res_chr)

tcell$cluster_res04_merged <- droplevels(factor(res_chr, levels = EXPECTED_CLUSTERS))
cat("\nMerged", paste(names(MERGE_MAP), collapse = ", "), "into",
    paste(unname(MERGE_MAP), collapse = ", "), ":", n_merged, "cells reassigned\n")
cat("Clusters remaining:", paste(levels(tcell$cluster_res04_merged), collapse = ", "), "\n")

# ---------------------------------------------------------------------
# 6e. NK / gamma-delta sanity check — the DISPROVEN hypothesis
# ---------------------------------------------------------------------
# Kept in the script (not deleted) so the record shows the check was run
# and what it returned. A disproven hypothesis you can point to is worth
# more than a clean script.
nk_genes <- intersect(c("KLRF1", "NCAM1", "TRDC", "TRGC1", "TRGC2",
                        "CD3D", "CD3E", "TRAC"), rownames(tcell))

nk_check <- fetch_expr(tcell, c(nk_genes, "cluster_res04_merged")) %>%
  group_by(cluster = cluster_res04_merged) %>%
  summarise(across(all_of(nk_genes), ~ round(mean(.x > 0), 3)), .groups = "drop")

write.csv(nk_check, file.path(RES_DIR, "nk_check_by_cluster.csv"), row.names = FALSE)
cat("\n--- NK / gamma-delta check: fraction of cells detecting each gene ---\n")
print(as.data.frame(nk_check))
cat("Expect KLRF1 and TRDC/TRGC near zero in every cluster.\n")

# ---------------------------------------------------------------------
# 6f. Filter Ig and TCR constant regions from marker tables
# ---------------------------------------------------------------------
# CSV-only operation: it does not touch the object, so nothing downstream
# needs redoing. This is the WHOLE ambient-Ig fix at the analysis level.
filter_marker_csv <- function(path) {
  if (!file.exists(path)) {
    warning("Marker file not found, skipping: ", path)
    return(invisible(NULL))
  }
  m  <- read.csv(path, stringsAsFactors = FALSE)
  gc <- if ("gene" %in% names(m)) "gene" else names(m)[1]

  drop <- grepl(ig_pattern, m[[gc]]) |
    m[[gc]] %in% c("JCHAIN", "MZB1", "TRAC", "TRBC1", "TRBC2")

  out <- sub("\\.csv$", "_filtered.csv", path)
  write.csv(m[!drop, ], out, row.names = FALSE)

  cat(sprintf("  %s: dropped %d of %d rows -> %s\n",
              basename(path), sum(drop), nrow(m), basename(out)))
  invisible(unique(m[[gc]][drop]))
}

if (DO_MARKER_FILTER) {
  cat("\n--- Filtering marker tables ---\n")
  dropped <- lapply(MARKER_CSVS, filter_marker_csv)
  cat("Genes removed:", paste(unique(unlist(dropped)), collapse = ", "), "\n")
}


# =====================================================================
# 7. WHAT CHANGED — before / after
# =====================================================================
counts_end <- table(tcell$cluster_res04_merged)

change_tbl <- tibble(
  cluster  = names(counts_start),
  n_before = as.integer(counts_start),
  n_after  = as.integer(counts_end[names(counts_start)])
) %>%
  mutate(
    n_after = tidyr::replace_na(n_after, 0L),
    delta   = n_after - n_before,
    note    = case_when(
      cluster %in% DROP_CLUSTERS     ~ "REMOVED - low-complexity ribosomal artifact",
      cluster %in% names(MERGE_MAP)  ~ paste("MERGED into cluster", MERGE_MAP[cluster]),
      cluster %in% unname(MERGE_MAP) ~ "gained merged cluster; lost Ig-high doublets",
      TRUE                           ~ "lost Ig-high doublets only"
    )
  )

write.csv(change_tbl, file.path(RES_DIR, "tcell_cleanup_changelog.csv"), row.names = FALSE)

cat("\n=====================================================================\n")
cat("CLEANUP SUMMARY\n")
cat("---------------------------------------------------------------------\n")
cat("Cells:            ", n_start, "->", ncol(tcell),
    sprintf("(-%d, %.1f%%)\n", n_start - ncol(tcell),
            100 * (n_start - ncol(tcell)) / n_start))
cat("  Ig-high doublets removed:", length(doublet_cells), "\n")
cat("  Artifact cluster removed:", n_dropped, "\n")
cat("Clusters:          ", length(EXPECTED_CLUSTERS), "->",
    nlevels(tcell$cluster_res04_merged), "\n")
cat("Cells reassigned by merge:", n_merged, "\n")
cat("Object modified: YES. UMAP/PCA/Harmony NOT recomputed -- coordinates\n")
cat("  are from the res-0.4 embedding. Dropping 1-4% of cells does not\n")
cat("  justify re-embedding, and re-running would renumber everything.\n")
cat("=====================================================================\n\n")
print(as.data.frame(change_tbl))


# =====================================================================
# 8. FINAL LABELS, DISPLAY CODES, PALETTE
# =====================================================================
# Named vector as a lookup table: names = old cluster IDs, values = new
# labels. map[as.character(old)] is the idiomatic R recode -- one
# vectorised pass, no loop.
label_map <- c(
  "0" = "CD8 TEM (GZMK/GZMH)",
  "1" = "CD4 naive/CM",
  "2" = "CD8 TEM (ZEB2-hi)",
  "3" = "Th17",
  "4" = "Treg",
  "5" = "CD4 resting/CM (TSHZ2+)",
  "6" = "CD8 cytotoxic (GZMB/GNLY)",
  "7" = "MAIT",
  "8" = "CD4 Tph-like (CXCL13+)"
)

new_labels <- unname(label_map)   # plot order = the order above

if (!setequal(names(label_map), levels(tcell$cluster_res04_merged))) {
  stop("label_map keys do not match surviving clusters: ",
       paste(setdiff(levels(tcell$cluster_res04_merged), names(label_map)),
             collapse = ", "))
}

# unname() is REQUIRED here. Subsetting a named vector by a character
# index returns names taken from the INDEX, not from the input -- so
# label_map[c("0","0","3")] comes back named "0","0","3". factor() copies
# names through, Seurat then tries to match them against cell barcodes,
# finds none, and errors with "No cell overlap between new meta data and
# Seurat object". unname() strips them so Seurat falls back on position.
tcell$cluster_id <- factor(
  unname(label_map[as.character(tcell$cluster_res04_merged)]),
  levels = new_labels
)
stopifnot(!any(is.na(tcell$cluster_id)))

# Short display codes C1-C9 for figures. Kept as a SEPARATE column:
# cluster_id stays the full annotation because the CSVs, stats tables and
# palette are all keyed to it.
code_map <- setNames(paste0("C", seq_along(new_labels)), new_labels)

tcell$cluster_code <- factor(
  unname(code_map[as.character(tcell$cluster_id)]),
  levels = unname(code_map)
)
stopifnot(!any(is.na(tcell$cluster_code)))

renaming_map <- tibble(
  res04_original = names(label_map),
  merged_from    = ifelse(names(label_map) %in% unname(MERGE_MAP),
                          paste(names(label_map), "+",
                                paste(names(MERGE_MAP), collapse = "+")),
                          names(label_map)),
  display_code   = unname(code_map[unname(label_map)]),
  final_label    = unname(label_map),
  n_cells        = as.integer(table(tcell$cluster_res04_merged)[names(label_map)])
)
write.csv(renaming_map, file.path(RES_DIR, "cluster_renaming_map.csv"), row.names = FALSE)
cat("\n--- Cluster renaming map ---\n")
print(as.data.frame(renaming_map))

# Palette keyed by LABEL, not by position, so colours stay bound to cell
# types even if new_labels is reordered later.
pal <- c(
  "CD8 TEM (GZMK/GZMH)"       = "#8B5E83",  # plum
  "CD4 naive/CM"              = "#A99B8E",  # taupe
  "CD8 TEM (ZEB2-hi)"         = "#E38A3B",  # orange
  "Th17"                      = "#E8C84A",  # gold
  "Treg"                      = "#C0392B",  # red
  "CD4 resting/CM (TSHZ2+)"   = "#EFA3AC",  # pink
  "CD8 cytotoxic (GZMB/GNLY)" = "#3E7A4E",  # green
  "MAIT"                      = "#3B5CA8",  # blue
  "CD4 Tph-like (CXCL13+)"    = "#6FA8B8"   # dusty teal
)
OTHER_COL <- "#D3D3D3"
stopifnot(all(new_labels %in% names(pal)))


# =====================================================================
# 9. COMPOSITION TABLES
# =====================================================================

# ---------------------------------------------------------------------
# 9a. Sample level — each biopsy is its own row
# ---------------------------------------------------------------------
# complete() comes BEFORE the totals: it fills sample x cluster
# combinations that had zero cells. Without it those rows are simply
# absent, and an absent row is not the same as a 0 -- it silently biases
# every mean and test that follows.
comp_sample <- tcell@meta.data %>%
  count(sample, donor, condition, cluster_id, name = "n_cells") %>%
  complete(nesting(sample, donor, condition), cluster_id,
           fill = list(n_cells = 0)) %>%
  group_by(sample) %>%
  mutate(unit_total = sum(n_cells),
         proportion = n_cells / unit_total) %>%
  ungroup()

write.csv(comp_sample, file.path(RES_DIR, "sample_x_cluster_long.csv"), row.names = FALSE)

comp_wide <- comp_sample %>%
  select(sample, condition, cluster_id, n_cells) %>%
  pivot_wider(names_from = cluster_id, values_from = n_cells, values_fill = 0) %>%
  arrange(condition, sample)
write.csv(comp_wide, file.path(RES_DIR, "sample_x_cluster_counts.csv"), row.names = FALSE)

comp_wide_prop <- comp_sample %>%
  select(sample, condition, cluster_id, proportion) %>%
  mutate(proportion = round(proportion, 4)) %>%
  pivot_wider(names_from = cluster_id, values_from = proportion, values_fill = 0) %>%
  arrange(condition, sample)
write.csv(comp_wide_prop, file.path(RES_DIR, "sample_x_cluster_proportions.csv"),
          row.names = FALSE)

cat("\n--- Sample x cluster counts ---\n")
print(as.data.frame(comp_wide))

# ---------------------------------------------------------------------
# 9b. Donor level — UNWEIGHTED mean across a donor's biopsies
# ---------------------------------------------------------------------
# Unweighted (mean of the biopsy proportions) rather than pooled (sum of
# counts / sum of totals): pooling would let a donor's largest biopsy
# dominate their own value, reintroducing exactly the size-driven
# weighting that moving to donor level is meant to remove.
comp_donor <- comp_sample %>%
  group_by(donor, condition, cluster_id) %>%
  summarise(proportion = mean(proportion),
            n_samples  = n_distinct(sample),
            n_cells    = sum(n_cells),
            unit_total = sum(unit_total),
            .groups    = "drop")

write.csv(comp_donor, file.path(RES_DIR, "donor_x_cluster_long.csv"), row.names = FALSE)

comp_donor_wide <- comp_donor %>%
  select(donor, condition, cluster_id, proportion) %>%
  mutate(proportion = round(proportion, 4)) %>%
  pivot_wider(names_from = cluster_id, values_from = proportion, values_fill = 0) %>%
  arrange(condition, donor)
write.csv(comp_donor_wide, file.path(RES_DIR, "donor_x_cluster_proportions.csv"),
          row.names = FALSE)

cat("\n--- Donor level ---\n")
print(as.data.frame(comp_donor %>% distinct(donor, condition) %>%
                      count(condition, name = "n_donors")))
print(as.data.frame(comp_donor_wide))


# =====================================================================
# 10. DIFFERENTIAL ABUNDANCE
# =====================================================================
# The unit of replication is the sample or the donor -- NEVER the cell.
# Testing pooled cell counts would treat ~21,000 cells as 21,000
# independent observations. (An earlier version of this analysis used SRR
# run IDs, which are sequencing runs from the same biopsy, i.e. technical
# replicates. Sample level is the floor.)
#
# 10e is the PRIMARY reported model: it keeps all 21 samples while
# modelling the fact that repeat biopsies from one donor are correlated.
# 10a-10d are reported as supplementary concordance checks.

da_list <- list()   # accumulator; one standardised tibble per test

# ---------------------------------------------------------------------
# 10a. Wilcoxon, sample level
# ---------------------------------------------------------------------
wilcox_sample <- wilcox_da(comp_sample, "sample")
write.csv(wilcox_sample, file.path(RES_DIR, "DA_wilcox_sample_level.csv"),
          row.names = FALSE)
cat("\n--- 10a. Wilcoxon, SAMPLE level (n =",
    n_distinct(comp_sample$sample), "samples), BH-adjusted ---\n")
print(as.data.frame(wilcox_sample), digits = 3)

da_list$wilcox_sample <- wilcox_sample %>% mutate(test_id = "wilcox_sample")

# ---------------------------------------------------------------------
# 10b. Wilcoxon, donor level
# ---------------------------------------------------------------------
wilcox_donor <- wilcox_da(comp_donor, "donor")
write.csv(wilcox_donor, file.path(RES_DIR, "DA_wilcox_donor_level.csv"),
          row.names = FALSE)
cat("\n--- 10b. Wilcoxon, DONOR level (n =",
    n_distinct(comp_donor$donor), "donors), BH-adjusted ---\n")
print(as.data.frame(wilcox_donor), digits = 3)

da_list$wilcox_donor <- wilcox_donor %>% mutate(test_id = "wilcox_donor")

# ---------------------------------------------------------------------
# 10c / 10d / 10e. propeller
# ---------------------------------------------------------------------
have_speckle <- requireNamespace("speckle", quietly = TRUE)
have_limma   <- requireNamespace("limma",   quietly = TRUE)

blocked <- NULL   # defined here so sections 10h / 11d can test for it

if (!have_speckle) {
  message("\nspeckle not installed -- skipping 10c/10d/10e/10h. Install with:\n",
          "  BiocManager::install('speckle')")
} else {

  # --- 10c. propeller, sample level ---------------------------------
  prop_sample <- propeller_da(clusters = tcell$cluster_id,
                              unit     = tcell$sample,
                              group    = tcell$condition,
                              level_label = "sample")
  write.csv(prop_sample, file.path(RES_DIR, "DA_propeller_sample_level.csv"),
            row.names = FALSE)
  cat("\n--- 10c. propeller (", PROP_TRANSFORM,
      ", empirical Bayes), SAMPLE level ---\n", sep = "")
  print(as.data.frame(prop_sample), digits = 3)

  da_list$propeller_sample <- prop_sample %>% mutate(test_id = "propeller_sample")

  # --- 10d. propeller, donor level ----------------------------------
  # NOTE: passing donor as the unit POOLS cells within a donor, so a
  # donor's biopsies are weighted by their cell numbers. That differs from
  # 10b, which averages the biopsy proportions unweighted. propeller needs
  # cell-level input to build its own count matrix, so this is the honest
  # way to run it at donor level -- just be consistent about which
  # definition you quote in the methods.
  prop_donor <- propeller_da(clusters = tcell$cluster_id,
                             unit     = tcell$donor,
                             group    = tcell$condition,
                             level_label = "donor")
  write.csv(prop_donor, file.path(RES_DIR, "DA_propeller_donor_level.csv"),
            row.names = FALSE)
  cat("\n--- 10d. propeller (", PROP_TRANSFORM,
      ", empirical Bayes), DONOR level (cells pooled within donor) ---\n", sep = "")
  print(as.data.frame(prop_donor), digits = 3)

  da_list$propeller_donor <- prop_donor %>% mutate(test_id = "propeller_donor")

  # --- 10e. PRIMARY: propeller, sample level, donor-blocked ----------
  # duplicateCorrelation() estimates one consensus within-donor
  # correlation (rho) by fitting a mixed model per cluster row, Fisher
  # z-averaging the nine estimates, and back-transforming. lmFit() then
  # uses it, so the 21 samples contribute as 21 observations with a known
  # dependency structure rather than 21 independent ones.
  #
  # rho enters at exactly ONE place: the `correlation` argument of lmFit,
  # which switches it from OLS to GLS with a block-diagonal correlation
  # matrix (1 on the diagonal, rho between samples sharing a donor, 0
  # elsewhere). It affects the STANDARD ERROR, not the effect estimate.
  if (!have_limma) {
    message("limma not installed -- skipping 10e.")
  } else {
    blocked <- tryCatch({
      props <- speckle::getTransformedProps(clusters  = tcell$cluster_id,
                                            sample    = tcell$sample,
                                            transform = PROP_TRANSFORM)

      # Align sample metadata to the COLUMN ORDER of the transformed
      # matrix. match() returns the positions of its first argument within
      # its second, so this reorders the rows to match. Getting it wrong
      # scrambles condition labels silently -- hence the identical() check.
      smeta <- tcell@meta.data %>% distinct(sample, condition, donor)
      smeta <- smeta[match(colnames(props$TransformedProps), smeta$sample), ]
      stopifnot(identical(as.character(smeta$sample),
                          colnames(props$TransformedProps)))

      grp    <- factor(as.character(smeta$condition), levels = c(GM_LABEL, PD_LABEL))
      design <- model.matrix(~ 0 + grp)
      colnames(design) <- c("GM", "PD")

      dupcor <- limma::duplicateCorrelation(props$TransformedProps,
                                            design = design,
                                            block  = smeta$donor)
      rho <- dupcor$consensus.correlation

      fit   <- limma::lmFit(props$TransformedProps, design,
                            block = smeta$donor, correlation = rho)
      contr <- limma::makeContrasts(contrasts = "PD-GM", levels = design)
      fit2  <- limma::eBayes(limma::contrasts.fit(fit, contr), robust = TRUE)

      # confint = TRUE adds CI.L / CI.R: moderated-t intervals built from
      # the POSTERIOR variance, so they are consistent with the empirical
      # Bayes p-value rather than an ordinary least-squares one.
      tt <- limma::topTable(fit2, coef = 1, number = Inf,
                            sort.by = "P", confint = TRUE)

      # Natural-scale means per condition, unweighted across samples.
      pmat  <- as.matrix(unclass(props$Proportions))
      mu_gm <- rowMeans(pmat[, smeta$condition == GM_LABEL, drop = FALSE])
      mu_pd <- rowMeans(pmat[, smeta$condition == PD_LABEL, drop = FALSE])

      # sort.by = "P" orders tt by p-value and rownames(tt) picks up that
      # same order, so labels stay paired with their statistics.
      cl <- rownames(tt)
      out <- tibble(
        cluster_id  = cl,
        n_GM        = sum(smeta$condition == GM_LABEL),
        n_PD        = sum(smeta$condition == PD_LABEL),
        mean_GM     = as.numeric(mu_gm[cl]),
        mean_PD     = as.numeric(mu_pd[cl]),
        median_GM   = NA_real_,
        median_PD   = NA_real_,
        statistic   = tt$t,
        p_value     = tt$P.Value,
        FDR         = tt$adj.P.Val,
        effect_type = paste0("logOdds_", PROP_TRANSFORM),
        effect      = tt$logFC,
        ci_low      = tt$CI.L,
        ci_high     = tt$CI.R,
        method      = "propeller_donor_blocked",
        level       = "sample (donor-blocked)"
      ) %>%
        mutate(sig = ifelse(!is.na(FDR) & FDR < FDR_CUTOFF, "*", ""))

      list(res = out, rho = rho, props = props, smeta = smeta,
           design = design, contr = contr)
    }, error = function(e) {
      message("10e failed: ", conditionMessage(e))
      NULL
    })

    if (!is.null(blocked)) {
      cat("\n--- 10e. PRIMARY MODEL: propeller, SAMPLE level, DONOR-BLOCKED ---\n")
      cat(sprintf("Consensus intra-donor correlation (duplicateCorrelation): %.3f\n",
                  blocked$rho))
      cat("  -> REPORT THIS. A correlation this size means repeat biopsies from\n")
      cat("     one donor are not independent, so the blocked model is required.\n")
      print(as.data.frame(blocked$res), digits = 3)

      write.csv(blocked$res, file.path(RES_DIR, "DA_propeller_donor_blocked.csv"),
                row.names = FALSE)
      writeLines(sprintf("consensus_intra_donor_correlation\t%.6f", blocked$rho),
                 file.path(RES_DIR, "intra_donor_correlation.txt"))

      da_list$propeller_blocked_sample <- blocked$res %>%
        mutate(test_id = "propeller_blocked_sample")
    }
  }
}

# ---------------------------------------------------------------------
# 10f. Combined results — every p-value in one CSV
# ---------------------------------------------------------------------
# bind_rows() stacks tibbles and fills missing columns with NA, which is
# why every helper returns the SAME column names. Standardising the output
# shape inside each helper is what makes this one line work.
da_all <- bind_rows(da_list) %>%
  mutate(display_code = unname(code_map[cluster_id])) %>%
  select(all_of(DA_COLS)) %>%
  arrange(test_id, p_value)

write.csv(da_all, file.path(RES_DIR, "DA_all_tests_combined.csv"), row.names = FALSE)

# Wide side-by-side view: one row per cluster, one p-value column per
# test. A cluster that moves in all five is a finding; one that only
# reaches significance in the most liberal test is not.
da_wide <- da_all %>%
  select(display_code, cluster_id, test_id, p_value, FDR) %>%
  pivot_wider(names_from = test_id, values_from = c(p_value, FDR),
              names_glue = "{test_id}_{.value}") %>%
  arrange(display_code)

write.csv(da_wide, file.path(RES_DIR, "DA_all_tests_wide.csv"), row.names = FALSE)

cat("\n=====================================================================\n")
cat("10f. ALL TESTS, SIDE BY SIDE\n")
cat("=====================================================================\n")
print(as.data.frame(da_wide), digits = 3)

# Effect sizes matter more than p-values at n = 18-21 with this much
# inter-donor variability. Print them together so the two are read jointly.
cat("\n--- Effect sizes (mean per-sample proportions, sample level) ---\n")
print(as.data.frame(
  wilcox_sample %>%
    transmute(display_code = unname(code_map[cluster_id]), cluster_id,
              mean_GM = round(mean_GM, 4), mean_PD = round(mean_PD, 4),
              fold_change = round(mean_PD / mean_GM, 2),
              log2FC = round(effect, 3)) %>%
    arrange(desc(log2FC))
))

# --- Which test supplies figure annotation ---------------------------
available <- unique(da_all$test_id)
chosen    <- TEST_FALLBACK[TEST_FALLBACK %in% available][1]
if (is.na(chosen)) stop("No differential abundance test produced results.")
if (chosen != PRIMARY_TEST) {
  message("PRIMARY_TEST '", PRIMARY_TEST, "' unavailable; using '", chosen, "'.")
}

primary_res  <- da_all %>% filter(test_id == chosen)
sig_clusters <- primary_res$cluster_id[primary_res$sig == "*"]

cat("\nFigure annotation from:", chosen, "\n")
cat("Significant at FDR <", FDR_CUTOFF, ":",
    if (length(sig_clusters)) paste(unname(code_map[sig_clusters]), collapse = ", ") else "none",
    "\n")

# ---------------------------------------------------------------------
# 10g. Leave-one-donor-out sensitivity
# ---------------------------------------------------------------------
# One donor contributes 3 of 8 PD samples. Any PD-enriched hit should
# survive dropping that donor, or be reported as donor-driven.
#
# The loop pattern (build a list inside the loop, bind_rows once at the
# end) is a fundamental worth writing yourself -- much faster and cleaner
# than growing a data frame with rbind() each iteration.
multi_donors <- comp_sample %>%
  distinct(sample, donor) %>%
  count(donor, name = "n_samples") %>%
  filter(n_samples > 1) %>%
  pull(donor)

cat("\n--- 10g. Sensitivity: donors with >1 biopsy:",
    paste(multi_donors, collapse = ", "), "---\n")

sens_list <- list()
for (d in multi_donors) {
  sens_list[[d]] <- wilcox_da(filter(comp_sample, donor != d), "sample") %>%
    transmute(cluster_id, dropped_donor = d,
              p_value_after_drop = p_value, FDR_after_drop = FDR)
}

loo_list <- list()
for (d in unique(comp_donor$donor)) {
  loo_list[[d]] <- wilcox_da(filter(comp_donor, donor != d), "donor") %>%
    transmute(cluster_id, dropped_donor = d, p_value_after_drop = p_value)
}

sens_all <- bind_rows(
  bind_rows(sens_list) %>% mutate(level = "sample"),
  bind_rows(loo_list)  %>% mutate(level = "donor", FDR_after_drop = NA_real_)
) %>%
  left_join(da_all %>% filter(test_id == "wilcox_sample") %>%
              select(cluster_id, p_value_full = p_value),
            by = "cluster_id") %>%
  mutate(display_code = unname(code_map[cluster_id])) %>%
  arrange(level, display_code, dropped_donor)

write.csv(sens_all, file.path(RES_DIR, "DA_leave_one_donor_out.csv"), row.names = FALSE)

loo_summary <- sens_all %>%
  group_by(level, display_code, cluster_id) %>%
  summarise(p_full = first(p_value_full),
            p_min  = min(p_value_after_drop, na.rm = TRUE),
            p_max  = max(p_value_after_drop, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(level, p_max)

write.csv(loo_summary, file.path(RES_DIR, "DA_leave_one_donor_out_summary.csv"),
          row.names = FALSE)
cat("\n--- Leave-one-donor-out p-value range per cluster ---\n")
print(as.data.frame(loo_summary), digits = 3)
cat("A hit whose p_max stays below 0.05 is robust to any single donor.\n")

# ---------------------------------------------------------------------
# 10h. Variance components + rho sensitivity sweep
# ---------------------------------------------------------------------
# rho rests on only 3 within-donor degrees of freedom (PD164 gives 2,
# PD161 gives 1, the other 16 donors give none), pooled across 9 clusters.
# There is no more precise estimator for this design. So rather than
# defending a point estimate, refit the model across the full range of
# possible rho and report that the conclusion does not depend on it.
if (!is.null(blocked)) {

  props  <- blocked$props
  smeta  <- blocked$smeta
  design <- blocked$design
  contr  <- blocked$contr
  cluster_names <- rownames(props$TransformedProps)

  RHO_GRID <- sort(unique(c(RHO_GRID_BASE, round(blocked$rho, 3))))

  # --- A. How much information does rho rest on? --------------------
  donor_sizes <- table(smeta$donor)
  within_df   <- sum(donor_sizes - 1)

  cat("\n=====================================================================\n")
  cat("10h-A. INFORMATION AVAILABLE FOR ESTIMATING rho\n")
  cat("---------------------------------------------------------------------\n")
  cat("Donors:", length(donor_sizes), "| Samples:", nrow(smeta), "\n")
  cat("Donors with repeat biopsies:",
      paste(names(donor_sizes)[donor_sizes > 1], collapse = ", "), "\n")
  print(as.data.frame(donor_sizes) %>%
          setNames(c("donor", "n_biopsies")) %>%
          mutate(within_donor_df = n_biopsies - 1) %>%
          filter(n_biopsies > 1))
  cat("Within-donor df per cluster:", within_df,
      "| pooled across", length(cluster_names), "clusters: ~",
      within_df * length(cluster_names), "df\n")
  cat("-> This is the ceiling. More precision needs more donors with repeat\n")
  cat("   biopsies, not a different estimator.\n")

  # --- B. Variance components, computed by hand ---------------------
  # Shows the raw scatter rho is derived from.
  #   within-donor variance  = pooled across biopsies of the same donor
  #   between-donor variance = variance of donor means
  #   ICC                    = between / (between + within)
  #
  # CAVEAT: var_between here pools BOTH conditions, so a real PD-vs-GM
  # difference inflates it and the ICC is biased upward for
  # condition-associated clusters. duplicateCorrelation does not have this
  # problem because it fits the design matrix first and works on
  # residuals. Quote 0.73 from duplicateCorrelation; treat these
  # per-cluster ICCs as supporting evidence only.
  tp <- props_to_df(props$TransformedProps, smeta)
  rep_donors <- names(donor_sizes)[donor_sizes > 1]

  vc <- lapply(cluster_names, function(cl) {
    v <- tp[[cl]]

    ss_w <- 0; df_w <- 0
    for (d in rep_donors) {
      x <- v[tp$donor == d]
      ss_w <- ss_w + sum((x - mean(x))^2)
      df_w <- df_w + length(x) - 1
    }
    var_within <- if (df_w > 0) ss_w / df_w else NA_real_

    donor_means <- tapply(v, tp$donor, mean)
    var_between <- var(donor_means)

    tibble(cluster_id   = cl,
           display_code = unname(code_map[cl]),
           sd_within    = sqrt(var_within),
           sd_between   = sqrt(var_between),
           var_within   = var_within,
           var_between  = var_between,
           ICC          = var_between / (var_between + var_within))
  }) %>% bind_rows() %>% arrange(display_code)

  write.csv(vc, file.path(RES_DIR, "intra_donor_variance_components.csv"),
            row.names = FALSE)

  cat("\n=====================================================================\n")
  cat("10h-B. VARIANCE COMPONENTS (", PROP_TRANSFORM, " scale)\n", sep = "")
  cat("---------------------------------------------------------------------\n")
  print(as.data.frame(vc), digits = 3)
  cat("\nPer-cluster ICC: median", round(median(vc$ICC, na.rm = TRUE), 3),
      "| range", paste(round(range(vc$ICC, na.rm = TRUE), 3), collapse = " - "), "\n")
  cat("Compare with duplicateCorrelation consensus:", round(blocked$rho, 3), "\n")

  # --- C. Effective sample size as a function of rho ----------------
  # k / (1 + (k-1) * rho) is the standard design effect: how many
  # independent observations a group of k correlated ones is worth.
  eff_n   <- function(rho) sum(donor_sizes / (1 + (donor_sizes - 1) * rho))
  eff_tbl <- tibble(rho = RHO_GRID,
                    effective_N = round(sapply(RHO_GRID, eff_n), 2),
                    nominal_N   = nrow(smeta))

  cat("\n=====================================================================\n")
  cat("10h-C. EFFECTIVE SAMPLE SIZE vs rho\n")
  cat("---------------------------------------------------------------------\n")
  print(as.data.frame(eff_tbl))
  cat("\nWith only", length(rep_donors), "repeat donors, effective N moves very\n")
  cat("little across the plausible range of rho -- which is why the exact\n")
  cat("value does not need to be pinned down.\n")

  # --- D. Refit across the rho grid ---------------------------------
  # lmFit() accepts `correlation` as a FIXED input, so rho can be set by
  # hand rather than estimated. Sweeping it replaces a fragile point
  # estimate with a robustness statement.
  fit_at_rho <- function(rho) {
    fit <- if (rho == 0) {
      limma::lmFit(props$TransformedProps, design)   # rho = 0 IS independence
    } else {
      limma::lmFit(props$TransformedProps, design,
                   block = smeta$donor, correlation = rho)
    }
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr), robust = TRUE)
    tt   <- limma::topTable(fit2, coef = 1, number = Inf,
                            sort.by = "none", confint = TRUE)

    tibble(rho          = rho,
           cluster_id   = rownames(tt),
           display_code = unname(code_map[rownames(tt)]),
           logFC        = tt$logFC,
           ci_low       = tt$CI.L,
           ci_high      = tt$CI.R,
           t            = tt$t,
           p_value      = tt$P.Value,
           FDR          = tt$adj.P.Val)
  }

  rho_sweep <- bind_rows(lapply(RHO_GRID, fit_at_rho)) %>%
    left_join(eff_tbl %>% select(rho, effective_N), by = "rho")

  write.csv(rho_sweep, file.path(RES_DIR, "DA_propeller_rho_sensitivity.csv"),
            row.names = FALSE)

  rho_wide_p <- rho_sweep %>%
    select(display_code, cluster_id, rho, p_value) %>%
    pivot_wider(names_from = rho, values_from = p_value, names_prefix = "p_rho_") %>%
    arrange(display_code)

  rho_wide_fdr <- rho_sweep %>%
    select(display_code, cluster_id, rho, FDR) %>%
    pivot_wider(names_from = rho, values_from = FDR, names_prefix = "FDR_rho_") %>%
    arrange(display_code)

  write.csv(left_join(rho_wide_p, rho_wide_fdr, by = c("display_code", "cluster_id")),
            file.path(RES_DIR, "DA_propeller_rho_sensitivity_wide.csv"),
            row.names = FALSE)

  cat("\n=====================================================================\n")
  cat("10h-D. RAW p-VALUE ACROSS THE rho GRID\n")
  cat("---------------------------------------------------------------------\n")
  print(as.data.frame(rho_wide_p), digits = 3)
  cat("\n--- BH-adjusted (FDR) across the rho grid ---\n")
  print(as.data.frame(rho_wide_fdr), digits = 3)

  rho_summary <- rho_sweep %>%
    group_by(display_code, cluster_id) %>%
    summarise(p_min = min(p_value), p_max = max(p_value),
              fdr_min = min(FDR),   fdr_max = max(FDR),
              logFC_min = min(logFC), logFC_max = max(logFC),
              .groups = "drop") %>%
    arrange(p_max)

  write.csv(rho_summary, file.path(RES_DIR, "DA_propeller_rho_sensitivity_summary.csv"),
            row.names = FALSE)
  cat("\n--- Range across the whole rho grid ---\n")
  print(as.data.frame(rho_summary), digits = 3)
  cat("\nlogFC barely moves: rho affects the STANDARD ERROR, not the effect\n")
  cat("estimate. rho = 0.99 is a numerical boundary (the block-correlation\n")
  cat("matrix approaches singularity), not a plausible correlation. It is\n")
  cat("retained deliberately and explained in the figure caption.\n")

  # --- E. Plot the sweep --------------------------------------------
  p_rho <- ggplot(rho_sweep, aes(x = rho, y = p_value, group = display_code)) +
    geom_hline(yintercept = FDR_CUTOFF, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_line(aes(colour = display_code), linewidth = 0.7) +
    geom_point(aes(colour = display_code), size = 1.6) +
    geom_vline(xintercept = blocked$rho, linetype = "dotted", colour = "grey20") +
    annotate("text", x = blocked$rho, y = 1,
             label = sprintf("estimated \u03c1 = %.2f", blocked$rho),
             hjust = 1.05, vjust = 1, size = 3, colour = "grey20") +
    scale_y_log10() +
    scale_colour_discrete(name = NULL) +
    theme_classic(base_size = 11) +
    labs(x = "assumed intra-donor correlation (\u03c1)",
         y = "raw p-value (log scale)",
         title = "Sensitivity of differential abundance to intra-donor correlation",
         subtitle = "\u03c1 = 0: samples independent | \u03c1 \u2192 1: biopsies interchangeable",
         caption = paste0("\u03c1 = 0.99 is a numerical boundary where the block-",
                          "correlation matrix approaches singularity;\nit lies ",
                          "outside the range supported by the data (per-cluster ",
                          "ICC 0.51-0.88)."))

  ggsave(file.path(FIG_DIR, "rho_sensitivity.png"), p_rho,
         width = 8, height = 5.5, dpi = 300)
  ggsave(file.path(FIG_DIR, "rho_sensitivity.svg"), p_rho,
         width = 8, height = 5.5)
}


# =====================================================================
# 11. FIGURES
# =====================================================================

# ---------------------------------------------------------------------
# 11a. UMAP — codes C1-C9 on the plot, colours keyed to annotation
# ---------------------------------------------------------------------
umap_df <- as.data.frame(Embeddings(tcell, REDUCTION))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster   <- tcell$cluster_id
umap_df$code      <- tcell$cluster_code
umap_df$condition <- tcell$condition

# cluster and code are one-to-one, so grouping by both keeps `code`
# available in the summarised frame without needing a join.
centroids <- umap_df %>%
  group_by(cluster, code) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# guide = "none" AT THE SCALE is the fix for the ggplot2 3.5+ guide crash.
# NoLegend() and theme(legend.position = "none") only hide the guide after
# it has been built, so they still trigger it.
p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_colour_manual(values = pal, guide = "none") +
  geom_text(data = centroids, aes(label = code), color = "black",
            size = 5, fontface = "bold") +
  theme_classic() +
  theme(axis.text = element_text(size = 8)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "All subjects")

ggsave(file.path(FIG_DIR, "umap_tcell_final.png"), p_umap,
       width = 7, height = 6.5, dpi = 300)
ggsave(file.path(FIG_DIR, "umap_tcell_final.svg"), p_umap,
       width = 7, height = 6.5)

# Legend version as a SEPARATE object. If the guide bug fires you still
# have the clean p_umap saved above. override.aes is needed because the
# legend would otherwise inherit size = 0.35.
legend_labels <- setNames(paste0(unname(code_map), " \u2014 ", new_labels), new_labels)

invisible(tryCatch({
  pl <- p_umap +
    scale_colour_manual(values = pal, labels = legend_labels, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.text = element_text(size = 8))
  ggsave(file.path(FIG_DIR, "umap_tcell_final_legend.png"), pl,
         width = 9.5, height = 6.5, dpi = 300)
  ggsave(file.path(FIG_DIR, "umap_tcell_final_legend.svg"), pl,
         width = 9.5, height = 6.5)
  pl
}, error = function(e) {
  message("Legend version failed (ggplot2 guide bug?): ", conditionMessage(e))
  NULL
}))

cat("\n--- Figure caption key ---\n")
cat(paste0(unname(code_map), " = ", new_labels, collapse = "; "), "\n")

# Facet on a COLUMN of umap_df, not on tcell$condition. Faceting on an
# external vector only works while its order matches the plotted rows --
# exactly the assumption subset() breaks. Everything a ggplot needs should
# live inside the data frame you pass it.
p_umap_split <- p_umap + facet_wrap(~ condition, ncol = 2) +
  ggtitle("T-cell clusters by condition")
ggsave(file.path(FIG_DIR, "umap_tcell_final_by_condition.png"), p_umap_split,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "umap_tcell_final_by_condition.svg"), p_umap_split,
       width = 12, height = 6)

# ---------------------------------------------------------------------
# 11b. UMAP density by condition
# ---------------------------------------------------------------------
# Background points WITHOUT a condition column -> drawn identically in
# every facet, giving all cells as spatial context behind each panel.
background_df <- umap_df[, c("UMAP_1", "UMAP_2")]

p_density <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = background_df, inherit.aes = FALSE,
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey62", size = 0.35, alpha = 0.55) +
  stat_density_2d_filled(aes(fill  = after_stat(level_mid),
                             alpha = after_stat(level_mid)),
                         color = NA, bins = 10, adjust = 0.4, n = 200) +
  stat_density_2d(color = "black", linewidth = 0.1, alpha = 0.25,
                  bins = 10, adjust = 0.4, n = 200) +
  scale_fill_gradientn(colours = viridisLite::turbo(256),
                       name = "Density", na.value = NA) +
  # range starting at 0 (not 0.05) so the lowest band is fully
  # transparent -- otherwise it paints a faint sheet whose edges follow
  # the DATA extent, not the panel rectangle. That is the tell for a drawn
  # layer versus a theme element.
  scale_alpha_continuous(range = c(0, 0.85), guide = "none") +
  facet_wrap(~ condition, ncol = 2) +
  labs(title = "T-cell spatial density by condition", x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 12) +
  theme(axis.title       = element_text(face = "bold", size = 11),
        axis.text        = element_text(size = 9),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 12),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position  = "right")

ggsave(file.path(FIG_DIR, "umap_density_by_condition.png"), p_density,
       width = 11, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "umap_density_by_condition.svg"), p_density,
       width = 11, height = 6, device = "svg")

# ---------------------------------------------------------------------
# 11c. Donut plots (GM and PD) — ALL 9 clusters
# ---------------------------------------------------------------------
# Donut percentages POOL all cells within a condition: "what fraction of
# PD T cells sit in this cluster". The statistics use per-sample /
# per-donor proportions. The two can differ when samples vary in size --
# expected, not a contradiction. Say which is which in the legend.
#
# Wedge order is FIXED to code_map order (C1 -> C9), not sorted by size,
# so the same cluster occupies the same angular position in both donuts.
# Sorting each independently would make them uncomparable by eye.
# complete() guards against a cluster being absent from one condition:
# an absent row would silently shift every wedge after it.
make_donut <- function(cond_label, display_name) {
  d <- tcell@meta.data %>%
    filter(condition == cond_label) %>%
    count(cluster_id, name = "n") %>%
    complete(cluster_id = factor(new_labels, levels = new_labels),
             fill = list(n = 0)) %>%
    mutate(pct   = 100 * n / sum(n),
           label = factor(as.character(cluster_id), levels = new_labels),
           code  = unname(code_map[as.character(cluster_id)])) %>%
    arrange(label)
  
  d$display <- paste0(d$code, "\n", sprintf("%.1f%%", d$pct))
  
  # geom_col stacks in REVERSE factor order, so a plain cumsum() puts
  # labels on the wrong wedges. Mirroring the cumulative midpoint fixes it.
  d <- d %>% mutate(ypos = sum(pct) - (cumsum(pct) - pct / 2))
  
  ggplot(d, aes(x = 2, y = pct, fill = label)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 3.6) +
    scale_fill_manual(values = pal, guide = "none") +
    geom_text(aes(x = 3.02, y = ypos, label = display),
              size = 2.9, lineheight = 0.9, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0, label = display_name,
             size = 6, fontface = "bold") +
    theme_void()
}

p_gm <- make_donut(GM_LABEL, "GM")
p_pd <- make_donut(PD_LABEL, "PD")

ggsave(file.path(FIG_DIR, "donut_GM.png"), p_gm, width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG_DIR, "donut_GM.svg"), p_gm, width = 6.5, height = 5.5)
ggsave(file.path(FIG_DIR, "donut_PD.png"), p_pd, width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG_DIR, "donut_PD.svg"), p_pd, width = 6.5, height = 5.5)
# Saved separately rather than combined: patchwork's plot assembly is the
# component crashing in this R install. Stack them in Illustrator/Inkscape.

# Numbers behind the donuts, so the figure legend can quote them
donut_tbl <- tcell@meta.data %>%
  count(condition, cluster_id, name = "n_cells") %>%
  group_by(condition) %>%
  mutate(pct_pooled = round(100 * n_cells / sum(n_cells), 2)) %>%
  ungroup() %>%
  mutate(display_code = unname(code_map[as.character(cluster_id)])) %>%
  select(condition, display_code, cluster_id, n_cells, pct_pooled) %>%
  arrange(condition, display_code)

write.csv(donut_tbl, file.path(RES_DIR, "donut_pooled_percentages.csv"), row.names = FALSE)
cat("\n--- Donut percentages (pooled cells within condition) ---\n")
print(as.data.frame(donut_tbl))

# =====================================================================
# 11d. DIFFERENTIAL ABUNDANCE — publication figure
# =====================================================================
# Layout notes, because this is the part that makes a figure look
# typeset rather than plotted:
#
# (1) Bars are drawn on the LOG-ODDS scale but the axis is LABELLED in
#     odds ratio. That keeps the geometry honest -- a 2-fold depletion
#     (OR 0.5) and a 2-fold enrichment (OR 2.0) sit equidistant from
#     OR 1. On a linear OR axis, depletion is squeezed into the 0-1
#     interval and reads as half as important as the equivalent
#     enrichment. This is the single most common way abundance figures
#     mislead.
#
# (2) The numeric columns live OUTSIDE the panel. coord_cartesian(clip =
#     "off") lets geoms draw beyond the panel limits; the right plot
#     margin makes room. Because the panel limits stop at the data, the
#     axis line does not run underneath the text.
#
# (3) y is a NUMERIC position, not a factor. Factors cannot carry a
#     header row -- with numbers, the headers are just another text layer
#     at y = n + 0.8. This is the trick worth remembering for any
#     table-beside-plot figure.

# --- Config for this figure -------------------------------------------
BAR_FULL_LABELS <- FALSE    # y axis: "C9 - CD4 Tph-like (CXCL13+)" vs just "C9"
BAR_SHOW_HEADERS <- TRUE

# --- Source the primary model ----------------------------------------
if (!is.null(blocked)) {
  da_src <- blocked$res
  rho_fig <- blocked$rho
} else if (file.exists(file.path(RES_DIR, "DA_propeller_donor_blocked.csv"))) {
  da_src <- read.csv(file.path(RES_DIR, "DA_propeller_donor_blocked.csv"),
                     stringsAsFactors = FALSE)
  rho_fig <- tryCatch(
    as.numeric(strsplit(readLines(
      file.path(RES_DIR, "intra_donor_correlation.txt"))[1], "\t")[[1]][2]),
    error = function(e) NA_real_)
} else {
  da_src <- NULL
  message("Blocked model unavailable -- skipping 11d.")
}

if (!is.null(da_src)) {
  
  stopifnot(all(c("ci_low", "ci_high") %in% names(da_src)))
  
  # p-value formatter. Journals want "<0.001" rather than "0.000", and a
  # consistent number of decimals down the column.
  fmt_p <- function(p, digits = 3) {
    ifelse(is.na(p), "-",
           ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = digits)))
  }
  
  da_plot <- da_src %>%
    mutate(
      display_code = unname(code_map[cluster_id]),
      odds_ratio   = exp(effect),
      or_low       = exp(ci_low),
      or_high      = exp(ci_high),
      direction    = ifelse(effect > 0, "PD", "GM"),
      # CI excluding zero on the log-odds scale == unadjusted p < 0.05.
      # Drives opacity only; the actual numbers are printed.
      ci_excl_zero = (ci_low > 0) | (ci_high < 0),
      or_txt  = sprintf("%.2f (%.2f\u2013%.2f)", odds_ratio, or_low, or_high),
      p_txt   = fmt_p(p_value),
      fdr_txt = fmt_p(FDR),
      y_lab   = if (BAR_FULL_LABELS) paste0(display_code, " \u2014 ", cluster_id)
      else display_code
    ) %>%
    # Ascending sort, then numeric y in that order, puts the LARGEST
    # effect at the top of a horizontal chart. The rev() normally needed
    # on a y-axis variable is unnecessary here because the sort is
    # ascending rather than descending -- same rule, one step earlier.
    arrange(effect) %>%
    mutate(ypos = row_number())
  
  n_cl <- nrow(da_plot)
  
  # --- Panel extent: data + CI only ----------------------------------
  e_rng <- range(c(da_plot$ci_low, da_plot$ci_high))
  pad   <- 0.06 * diff(e_rng)
  x_lo  <- min(0, e_rng[1] - pad)
  x_hi  <- e_rng[2] + pad
  span  <- x_hi - x_lo
  
  # Column x positions, in data units, to the right of the panel
  x_or  <- x_hi + 0.14 * span
  x_p   <- x_hi + 0.60 * span
  x_fdr <- x_hi + 0.80 * span
  
  # Odds-ratio ticks. Over-complete on purpose: ggplot silently drops any
  # break outside the panel, so the same vector works whatever the range.
  or_breaks <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6)
  keep_br   <- log(or_breaks) >= x_lo & log(or_breaks) <= x_hi
  
  p_da <- ggplot(da_plot, aes(x = effect, y = ypos)) +
    # Reference line first so bars sit on top of it
    geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
    geom_col(aes(fill = direction, alpha = ci_excl_zero),
             width = 0.62, orientation = "y") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.20, linewidth = 0.45, colour = "grey20") +
    # Numeric columns
    geom_text(aes(x = x_or,  label = or_txt),  hjust = 0.5, size = 3.0,
              colour = "grey15", family = "sans") +
    geom_text(aes(x = x_p,   label = p_txt),   hjust = 0.5, size = 3.0,
              colour = "grey15", family = "sans") +
    geom_text(aes(x = x_fdr, label = fdr_txt),
              fontface = ifelse(da_plot$FDR < 0.10, "bold", "plain"),
              hjust = 0.5, size = 3.0, colour = "grey15", family = "sans") +
    scale_fill_manual(values = c(PD = "#A5231C", GM = "#1E7A32"), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.42), guide = "none") +
    scale_x_continuous(name   = "odds ratio (PD vs GM), log scale",
                       breaks = log(or_breaks[keep_br]),
                       labels = format(or_breaks[keep_br], trim = TRUE,
                                       drop0trailing = TRUE)) +
    scale_y_continuous(breaks = da_plot$ypos, labels = da_plot$y_lab,
                       expand = expansion(add = c(0.6, 1.1))) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y   = element_text(size = 10, colour = "black", hjust = 1),
      axis.text.x   = element_text(size = 9,  colour = "black"),
      axis.title.x  = element_text(size = 10, margin = margin(t = 6)),
      axis.line.y   = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y  = element_line(colour = "black", linewidth = 0.4),
      plot.title    = element_text(face = "bold", size = 12,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, colour = "grey35",
                                   margin = margin(b = 10)),
      plot.caption  = element_text(size = 7.2, colour = "grey35", hjust = 0,
                                   margin = margin(t = 10)),
      plot.margin   = margin(t = 8, r = 128, b = 8, l = 8)
    ) +
    labs(
      y = NULL,
      title = "T-cell cluster abundance in periodontitis",
      subtitle = sprintf(paste0("propeller on logit proportions, donor-blocked ",
                                "(\u03c1 = %.2f); n = 21 samples / 18 donors"),
                         rho_fig),
      caption = paste0(
        "Bars: odds ratio for each cluster in PD relative to healthy gingiva. ",
        "Whiskers: unadjusted 95% moderated-t CI.\n",
        "p = nominal; FDR = Benjamini\u2013Hochberg across the 9 clusters. ",
        "Faded bars have a CI spanning OR = 1.\n",
        "Red = enriched in PD, green = enriched in GM."
      )
    )
  
  # Column headers as a separate text layer above the top row.
  if (BAR_SHOW_HEADERS) {
    hdr <- tibble(
      x     = c(x_or, x_p, x_fdr),
      y     = n_cl + 0.85,
      lab   = c("OR (95% CI)", "p", "FDR"),
      face  = c("bold", "bold.italic", "bold")
    )
    p_da <- p_da +
      geom_text(data = hdr, aes(x = x, y = y, label = lab),
                fontface = hdr$face, hjust = 0.5, size = 3.1,
                colour = "black", inherit.aes = FALSE)
  }
  
  # 300 dpi TIFF/PNG for submission, SVG for editing in Illustrator.
  ggsave(file.path(FIG_DIR, "DA_bar_blocked.png"), p_da,
         width = 9.2, height = 4.6, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, "DA_bar_blocked.svg"), p_da,
         width = 9.2, height = 4.6, bg = "white")
  ggsave(file.path(FIG_DIR, "DA_bar_blocked.pdf"), p_da,
         width = 9.2, height = 4.6, bg = "white")
  
  # Values behind the figure
  da_table <- da_plot %>%
    transmute(display_code, cluster_label = cluster_id,
              log_odds   = round(effect, 3),
              odds_ratio = round(odds_ratio, 3),
              OR_ci_low  = round(or_low, 3),
              OR_ci_high = round(or_high, 3),
              mean_prop_GM = round(mean_GM, 4),
              mean_prop_PD = round(mean_PD, 4),
              p_nominal  = signif(p_value, 4),
              FDR_BH     = signif(FDR, 4)) %>%
    arrange(p_nominal)
  
  write.csv(da_table, file.path(RES_DIR, "DA_bar_blocked_values.csv"),
            row.names = FALSE)
  cat("\n--- 11d. Values behind the DA figure ---\n")
  print(as.data.frame(da_table))
  cat("NOTE: log_odds is on the LOGIT scale; exp(log_odds) = odds_ratio.\n")
  cat("      Never label either as log2 fold change.\n")
}

# ---------------------------------------------------------------------
# 11e. Per-sample and per-donor proportion plots
# ---------------------------------------------------------------------
# The figure showing the DISTRIBUTION behind every p-value above. At
# n = 21 samples / 18 donors this is more informative than any bar, and it
# makes donor-driven signals visible by eye. For C9 specifically this is
# the more persuasive main-figure panel -- reviewers trust points they can
# count more than they trust a bar.
prop_plot_df <- comp_sample %>%
  mutate(display_code = factor(unname(code_map[as.character(cluster_id)]),
                               levels = unname(code_map)),
         condition     = factor(condition, levels = c(GM_LABEL, PD_LABEL)))

cond_fill <- setNames(c("#1E7A32", "#A5231C"), c(GM_LABEL, PD_LABEL))

p_props <- ggplot(prop_plot_df, aes(x = condition, y = 100 * proportion,
                                    fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85, shape = 21, colour = "grey20") +
  facet_wrap(~ display_code, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = cond_fill, guide = "none") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold")) +
  labs(x = NULL, y = "% of sample's T cells",
       title = "Per-sample cluster proportions",
       subtitle = "One point per biopsy (n = 21). Repeat biopsies not collapsed.")

ggsave(file.path(FIG_DIR, "per_sample_proportions.png"), p_props,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "per_sample_proportions.svg"), p_props,
       width = 9, height = 7)

p_props_donor <- prop_plot_df %>%
  group_by(donor, condition, display_code) %>%
  summarise(proportion = mean(proportion), .groups = "drop") %>%
  ggplot(aes(x = condition, y = 100 * proportion, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.85, shape = 21, colour = "grey20") +
  facet_wrap(~ display_code, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = cond_fill, guide = "none") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold")) +
  labs(x = NULL, y = "% of donor's T cells",
       title = "Per-donor cluster proportions",
       subtitle = "One point per donor (n = 18), repeat biopsies averaged unweighted.")

ggsave(file.path(FIG_DIR, "per_donor_proportions.png"), p_props_donor,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "per_donor_proportions.svg"), p_props_donor,
       width = 9, height = 7)

# ---------------------------------------------------------------------
# 11f. Grouped marker dot plot
# ---------------------------------------------------------------------
cluster_col   <- "cluster_code"
cluster_order <- paste0("C", c(1, 3, 7, 8, 2, 4, 5, 6, 9))  # CD8 lobe, then CD4
stopifnot(setequal(cluster_order, levels(tcell$cluster_code)))

gene_groups <- list(
  "Lineage"                          = c("CD4", "CD8A", "CD3D"),
  "Naive, Memory\n& Differentiation" = c("IL7R", "CCR7", "KLF2", "LEF1", "TCF7", "SELL"),
  "Activation &\nCostimulation"      = c("CD69", "ICOS", "CD40LG", "TNFRSF4", "IL2RB"),
  "Migration &\nTissue Residency"    = c("CCR6", "CXCR6", "ITGAE", "CX3CR1"),
  "Immune Checkpoint\n& Exhaustion"  = c("PDCD1", "TIGIT", "TOX2", "CXCL13"),
  "Regulatory"                       = c("FOXP3", "CTLA4", "IL2RA", "TNFRSF18"),
  "Innate-like\n& Cytotoxic"         = c("KLRB1", "GZMB", "GZMK", "NKG7", "ZBED2")
)

# stack() turns a named list into a two-column data frame (values, ind) --
# the quickest way to get a gene -> group lookup out of a list.
gene_map <- utils::stack(gene_groups) %>% rename(gene = values, gene_group = ind)

missing <- setdiff(gene_map$gene, rownames(tcell))
if (length(missing)) message("Not in object, dropped: ", paste(missing, collapse = ", "))
gene_map <- gene_map %>% filter(gene %in% rownames(tcell))

expr <- fetch_expr(tcell, c(gene_map$gene, cluster_col))
names(expr)[names(expr) == cluster_col] <- "cluster"

dot_df <- expr %>%
  pivot_longer(all_of(gene_map$gene), names_to = "gene", values_to = "value") %>%
  group_by(cluster, gene) %>%
  summarise(pct_exp = 100 * mean(value > 0),
            avg_exp = mean(expm1(value)), .groups = "drop") %>%
  group_by(gene) %>%
  # Zero-variance guard: scale() on a constant vector returns NaN and the
  # dots silently disappear rather than erroring.
  mutate(avg_scaled = if (sd(log1p(avg_exp)) > 0) as.numeric(scale(log1p(avg_exp))) else 0) %>%
  ungroup() %>%
  left_join(gene_map, by = "gene") %>%
  mutate(cluster    = factor(cluster, levels = cluster_order),
         gene       = factor(gene, levels = rev(gene_map$gene)),
         gene_group = factor(gene_group, levels = names(gene_groups)))

p_dot <- ggplot(dot_df, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct_exp, colour = avg_scaled)) +
  scale_size_continuous(name = "% Expressing", limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100), range = c(0, 5)) +
  scale_colour_distiller(name = "Avg Expression", palette = "YlOrRd",
                         direction = 1, limits = c(-2, 2),
                         breaks = c(-2, -1, 0, 1, 2), oob = scales::squish) +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(position = "right") +
  theme_bw(base_size = 10) +
  theme(panel.grid        = element_blank(),
        panel.background  = element_rect(fill = "white", colour = NA),
        plot.background   = element_rect(fill = "white", colour = NA),
        panel.spacing.y   = unit(2, "pt"),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left = element_text(angle = 90, size = 6.5),
        axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y.right = element_text(face = "italic", size = 8),
        axis.title        = element_blank(),
        legend.position   = "bottom")

ggsave(file.path(FIG_DIR, "dotplot_grouped.png"), p_dot,
       width = 5.5, height = 9.5, dpi = 300)
ggsave(file.path(FIG_DIR, "dotplot_grouped.svg"), p_dot,
       width = 5.5, height = 9.5)


# =====================================================================
# 12. SAVE
# =====================================================================
saveRDS(tcell, file.path(RES_DIR, "tcell_final_seurat_object.rds"))

cat("\n=====================================================================\n")
cat("FILES WRITTEN\n")
cat("---------------------------------------------------------------------\n")
cat("Figures  (", FIG_DIR, "):\n", sep = "")
cat(paste0("  ", sort(list.files(FIG_DIR)), collapse = "\n"), "\n")
cat("Results  (", RES_DIR, "):\n", sep = "")
cat(paste0("  ", sort(list.files(RES_DIR)), collapse = "\n"), "\n")
cat("=====================================================================\n")
cat("PRIMARY RESULT -> DA_propeller_donor_blocked.csv (+ DA_bar_blocked_values.csv)\n")
cat("SUPPLEMENTARY  -> DA_all_tests_wide.csv, DA_propeller_rho_sensitivity*.csv,\n")
cat("                  DA_leave_one_donor_out_summary.csv,\n")
cat("                  intra_donor_variance_components.csv\n")
cat("=====================================================================\n\n")

sessionInfo()


# =====================================================================
# 13. DIFFERENTIAL EXPRESSION: PD vs GM in T cells
#     Pseudobulk + limma-voom, donor-blocked. Volcano plot.
#
# Append to human_tcell_final_annotation_v4.R after section 12 (or before the final
# saveRDS if you prefer the object saved last). Needs from earlier:
#   tcell, GM_LABEL, PD_LABEL, FIG_DIR, RES_DIR, code_map, ig_pattern
#
# WHY PSEUDOBULK AND NOT FindMarkers()
# Seurat's default Wilcoxon DE treats every CELL as an independent
# observation. With ~21,600 cells from 18 donors that inflates the
# effective N by three orders of magnitude and returns p-values in the
# 1e-50 range for genes with trivial effects. It is the same
# pseudoreplication error as testing cell-level proportions, and it is
# now widely flagged in review (Squair et al. 2021, Nat Commun; Zimmerman
# et al. 2021; Murphy & Skene 2022).
#
# Aggregating counts to one profile per SAMPLE, then running limma-voom
# with donor as a blocking factor, keeps the DE analysis on the same
# statistical footing as the differential abundance model in 10e. Same
# unit of replication, same rho, same defensibility.
# =====================================================================

DE_ASSAY        <- "RNA"     # raw counts, NOT integrated / SCT
DE_MIN_SAMPLES  <- 6         # filterByExpr min number of samples
DE_LFC_CUTOFF   <- 1         # log2FC line on the volcano (annotation only)
DE_FDR_CUTOFF   <- 0.05
DE_N_LABEL      <- 20        # genes to label on the volcano

# Restrict to one cluster instead of all T cells by setting this to a
# display code, e.g. "C9". NULL = whole T-cell compartment.
DE_SUBSET_CODE  <- NULL

# Ig genes are ambient in this dataset (median percent.ig 0.17%, but
# donor-driven: PD164b reaches 3.8%). They will surface as apparent
# condition effects that are really one donor's soup. Flagged and written
# to a separate file rather than silently dropped, so the decision is
# visible and reversible.
DE_FLAG_IG      <- TRUE

stopifnot(requireNamespace("limma", quietly = TRUE),
          requireNamespace("edgeR", quietly = TRUE))

# ---------------------------------------------------------------------
# 13a. Subset (optional)
# ---------------------------------------------------------------------
if (!is.null(DE_SUBSET_CODE)) {
  keep_cells <- colnames(tcell)[as.character(tcell$cluster_code) == DE_SUBSET_CODE]
  stopifnot(length(keep_cells) > 0)
  de_obj <- subset(tcell, cells = keep_cells)
  de_tag <- paste0("_", DE_SUBSET_CODE)
  cat("\nDE restricted to", DE_SUBSET_CODE, ":", ncol(de_obj), "cells\n")
} else {
  de_obj <- tcell
  de_tag <- "_allTcells"
  cat("\nDE across all T cells:", ncol(de_obj), "cells\n")
}

# ---------------------------------------------------------------------
# 13b. Pseudobulk: sum raw counts within each sample
# ---------------------------------------------------------------------
# Done by matrix multiplication rather than AggregateExpression() because
# the Seurat API for this has changed repeatedly across v4/v5. A sparse
# indicator matrix (cells x samples, 1 where a cell belongs to a sample)
# turns the aggregation into one matrix product -- fast, and immune to
# version drift.
#
# This is a fundamental worth understanding: counts %*% indicator sums
# columns by group. The same trick does group means if you divide each
# indicator column by its column sum.
cnt <- tryCatch(
  SeuratObject::GetAssayData(de_obj, assay = DE_ASSAY, layer = "counts"),
  error = function(e) SeuratObject::GetAssayData(de_obj, assay = DE_ASSAY, slot = "counts")
)

sample_vec <- factor(de_obj$sample)
ind <- Matrix::sparse.model.matrix(~ 0 + sample_vec)
colnames(ind) <- levels(sample_vec)

pb <- as.matrix(cnt %*% ind)          # genes x samples
stopifnot(ncol(pb) == nlevels(sample_vec))

cat("Pseudobulk matrix:", nrow(pb), "genes x", ncol(pb), "samples\n")
cat("Library sizes (total counts per sample):\n")
print(summary(colSums(pb)))

# Sample metadata, aligned to the COLUMNS of pb
smeta_de <- de_obj@meta.data %>%
  distinct(sample, condition, donor) %>%
  mutate(sample = as.character(sample))
smeta_de <- smeta_de[match(colnames(pb), smeta_de$sample), ]
stopifnot(identical(smeta_de$sample, colnames(pb)))

# Cells contributing per sample -- a sample with very few cells makes a
# noisy pseudobulk profile. Report it; consider dropping samples below
# ~30 cells if you subset to a small cluster.
cells_per_sample <- as.integer(table(sample_vec)[colnames(pb)])
cat("\nCells per sample: min", min(cells_per_sample),
    "| median", median(cells_per_sample), "| max", max(cells_per_sample), "\n")

# ---------------------------------------------------------------------
# 13c. Filter, normalise, voom
# ---------------------------------------------------------------------
grp_de <- factor(as.character(smeta_de$condition), levels = c(GM_LABEL, PD_LABEL))
design_de <- model.matrix(~ 0 + grp_de)
colnames(design_de) <- c("GM", "PD")

dge <- edgeR::DGEList(counts = pb, samples = smeta_de, group = grp_de)

# filterByExpr keeps genes with enough counts in enough samples, using the
# design so it does not require expression in BOTH groups. Removing
# low-count genes matters: they dominate the multiple-testing burden and
# have unstable variance estimates.
keep <- edgeR::filterByExpr(dge, design = design_de,
                            min.count = 10, min.total.count = 15)
cat("\nGenes retained by filterByExpr:", sum(keep), "of", length(keep), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalisation: corrects for composition bias, i.e. a few very highly
# expressed genes soaking up library depth in some samples.
dge <- edgeR::calcNormFactors(dge, method = "TMM")

# voom estimates the mean-variance relationship of log-CPM and turns it
# into precision weights, which is what lets a linear model be used on
# count data. Save the plot -- a badly behaved voom curve is the first
# sign the pseudobulk profiles are too shallow.
png(file.path(FIG_DIR, paste0("DE_voom_meanvar", de_tag, ".png")),
    width = 1400, height = 1100, res = 180)
v <- limma::voom(dge, design_de, plot = TRUE)
dev.off()

# ---------------------------------------------------------------------
# 13d. Donor blocking — same structure as 10e
# ---------------------------------------------------------------------
# Two-pass voom is the standard limma recipe when blocking: estimate the
# correlation from a first-pass voom, then re-run voom aware of it, then
# fit. Skipping the second pass leaves the weights slightly wrong.
dupcor_de <- limma::duplicateCorrelation(v$E, design_de, block = smeta_de$donor)
rho_de    <- dupcor_de$consensus.correlation
cat("\nIntra-donor correlation (expression, pseudobulk):",
    round(rho_de, 3), "\n")
cat("Compare with the abundance model's rho. They need not match -- one is\n")
cat("about composition, the other about transcriptional state.\n")

v2 <- limma::voom(dge, design_de, block = smeta_de$donor,
                  correlation = rho_de, plot = FALSE)

fit_de  <- limma::lmFit(v2, design_de, block = smeta_de$donor,
                        correlation = rho_de)
contr_de <- limma::makeContrasts(contrasts = "PD-GM", levels = design_de)
fit2_de  <- limma::eBayes(limma::contrasts.fit(fit_de, contr_de), robust = TRUE)

res_de <- limma::topTable(fit2_de, coef = 1, number = Inf,
                          sort.by = "P", confint = TRUE) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    is_ig   = grepl(ig_pattern, gene),
    is_ribo = grepl("^RP[LS]", gene),
    is_mito = grepl("^MT-", gene),
    direction = case_when(
      adj.P.Val < DE_FDR_CUTOFF & logFC >  DE_LFC_CUTOFF ~ paste("up in", PD_LABEL),
      adj.P.Val < DE_FDR_CUTOFF & logFC < -DE_LFC_CUTOFF ~ paste("up in", GM_LABEL),
      adj.P.Val < DE_FDR_CUTOFF                          ~ "significant, small effect",
      TRUE                                               ~ "not significant"
    )
  )

write.csv(res_de, file.path(RES_DIR, paste0("DE_PD_vs_GM", de_tag, ".csv")),
          row.names = FALSE)

cat("\n--- DE summary (FDR <", DE_FDR_CUTOFF, ") ---\n")
print(table(res_de$direction))
cat("\nTop 20 by p-value:\n")
print(as.data.frame(res_de %>%
                      select(gene, logFC, CI.L, CI.R, P.Value, adj.P.Val,
                             is_ig, is_ribo) %>%
                      head(20)), digits = 3)

if (DE_FLAG_IG) {
  ig_hits <- res_de %>% filter(is_ig, adj.P.Val < DE_FDR_CUTOFF)
  if (nrow(ig_hits) > 0) {
    write.csv(ig_hits, file.path(RES_DIR, paste0("DE_ig_flagged", de_tag, ".csv")),
              row.names = FALSE)
    cat("\n!! ", nrow(ig_hits), " significant Ig genes. These track ambient\n",
        "   contamination, which in this dataset is DONOR-driven (PD164b\n",
        "   3.8% vs cohort median 0.17%), not condition-driven. Excluded from\n",
        "   the volcano below and from interpretation. See DE_ig_flagged*.csv\n", sep = "")
    print(as.data.frame(ig_hits %>% select(gene, logFC, adj.P.Val)), digits = 3)
  }
}

# ---------------------------------------------------------------------
# 13e. Volcano plot
# ---------------------------------------------------------------------
# y = -log10(raw p), colour = FDR-based call. Plotting raw p on y with an
# FDR-based colour is the convention: it preserves the spread of the
# y-axis while the significance call remains corrected.
volc_df <- res_de %>%
  { if (DE_FLAG_IG) filter(., !is_ig) else . } %>%
  mutate(neg_log10_p = -log10(P.Value))

# Label the strongest hits, excluding ribosomal and mitochondrial genes,
# which are technical here (see the cluster 9 diagnosis in section 5).
lab_df <- volc_df %>%
  filter(adj.P.Val < DE_FDR_CUTOFF, abs(logFC) > DE_LFC_CUTOFF,
         !is_ribo, !is_mito) %>%
  slice_head(n = DE_N_LABEL)

de_cols <- setNames(
  c("#A5231C", "#1E7A32", "#8C8C8C", "#D3D3D3"),
  c(paste("up in", PD_LABEL), paste("up in", GM_LABEL),
    "significant, small effect", "not significant")
)

p_volc <- ggplot(volc_df, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(colour = direction), size = 1.1, alpha = 0.8) +
  geom_vline(xintercept = c(-DE_LFC_CUTOFF, DE_LFC_CUTOFF),
             linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_hline(yintercept = -log10(max(volc_df$P.Value[volc_df$adj.P.Val < DE_FDR_CUTOFF],
                                     na.rm = TRUE)),
             linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  scale_colour_manual(values = de_cols, name = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        plot.caption  = element_text(size = 7.5, colour = "grey35", hjust = 0)) +
  labs(
    x = expression(log[2]~"fold change (PD / GM)"),
    y = expression(-log[10]~"(p-value)"),
    title = paste0("Differential expression, PD vs GM",
                   if (is.null(DE_SUBSET_CODE)) " (all T cells)"
                   else paste0(" (", DE_SUBSET_CODE, ")")),
    subtitle = sprintf("Pseudobulk limma-voom, donor-blocked (\u03c1 = %.2f)", rho_de),
    caption = paste0(
      "One profile per sample (n = ", ncol(pb), " samples, ",
      n_distinct(smeta_de$donor), " donors). TMM-normalised, voom-weighted.\n",
      "Dashed lines: |log2FC| = ", DE_LFC_CUTOFF, " and the p-value at which ",
      "FDR = ", DE_FDR_CUTOFF, ".\n",
      if (DE_FLAG_IG) "Immunoglobulin genes excluded (ambient, donor-driven).\n" else "",
      "Colour reflects BH-adjusted significance; the y-axis shows unadjusted p."
    )
  )

# ggrepel keeps labels off the points. Optional dependency -- the plot
# still saves without it, just unlabelled.
if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_df) > 0) {
  p_volc <- p_volc +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene),
                             size = 2.9, fontface = "italic",
                             max.overlaps = 30, box.padding = 0.35,
                             min.segment.length = 0.1, segment.colour = "grey60")
} else if (nrow(lab_df) == 0) {
  message("No genes pass both FDR and logFC cutoffs -- volcano left unlabelled.")
} else {
  message("ggrepel not installed; volcano saved without gene labels.")
}

ggsave(file.path(FIG_DIR, paste0("DE_volcano_PD_vs_GM", de_tag, ".png")),
       p_volc, width = 8, height = 6.5, dpi = 300)
ggsave(file.path(FIG_DIR, paste0("DE_volcano_PD_vs_GM", de_tag, ".svg")),
       p_volc, width = 8, height = 6.5)

cat("\nWritten:\n")
cat("  results/DE_PD_vs_GM", de_tag, ".csv\n", sep = "")
cat("  qc_plots/final/DE_volcano_PD_vs_GM", de_tag, ".[png|svg]\n", sep = "")
cat("  qc_plots/final/DE_voom_meanvar", de_tag, ".png\n", sep = "")








