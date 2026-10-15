#' ---------------------------------------------------------------------------
#' Module 06 — Cross-species Tph: marker overlap and specificity matrix
#'
#' Question: do the human C9 "CD4 Tph-like (CXCL13+)" marker list and the mouse
#' "Tfh/Tph-like" marker list share more genes than chance -- and is that pair
#' the BEST match in the full cluster-by-cluster matrix, or just one significant
#' cell among many?
#'
#' The hypergeometric test is the urn problem. A bag holds every gene that could
#' have been a marker (N). K of them are human C9 markers. Draw n genes (the
#' mouse list size). Then expected = K*n/N, fold = observed/expected, and p is
#' the chance of drawing >= observed by chance.
#'
#' The answer depends almost entirely on N, so THREE universes are reported side
#' by side rather than one being chosen. Two marker-table versions are also
#' reported (as-published, and recomputed at matched parameters) because the two
#' saved tables used different min.pct and the mouse one was top-500 truncated.
#'
#' Input:  ../../Human/tcell_final_seurat_object.rds
#'         ../../Human/tcell_subset_markers.csv
#'         ../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds
#'         ../../03_Tcell_Analysis/objects/markers_TCRab_res0.8_top500.csv
#' Output: ../objects/markers_matched_{human,mouse}.csv
#'         ../objects/overlap_matrix_<version>_<universe>.csv
#'         ../objects/reciprocal_best_hits.csv
#'         ../objects/anchor_gene_specificity.csv
#'         ../objects/universe_sizes.csv
#'         ../figures/overlap_fold_heatmap.svg
#'         ../results/sessionInfo_01.txt
#'
#' Run from 06_CrossSpecies_Tph/scripts/:  Rscript 01_marker_overlap_specificity.R
#' ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(babelgene)
  library(svglite)
})

set.seed(1234)

HUMAN_RDS    <- "../../Human/tcell_final_seurat_object.rds"
HUMAN_MARK   <- "../../Human/tcell_subset_markers.csv"
MOUSE_RDS    <- "../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds"
MOUSE_MARK   <- "../../03_Tcell_Analysis/objects/markers_TCRab_res0.8_top500.csv"
OBJ <- "../objects"; FIG <- "../figures"; RES <- "../results"

# Matched recomputation parameters. logfc.threshold 0.25 = 2^0.25 = 1.19x, a 19%
# difference -- appropriate for T-cell-subset-vs-other-T-cells. Raising it to 0.5
# drops RORA, JUNB and DOCK8 (3 of the 13 anchors) from the human C9 list on
# threshold alone, so it is deliberately NOT raised.
MATCH_MINPCT <- 0.25
MATCH_LFC    <- 0.25
PADJ         <- 0.05

# Human raw cluster (0-8, post contamination removal) -> display code.
# Human-data-analysis step 06 drops cluster 9 (ribosomal artifact) and merges 10 into 0.
HUMAN_CODE <- c("0"="C1","1"="C2","2"="C3","3"="C4","4"="C5",
                "5"="C6","6"="C7","7"="C8","8"="C9")
HUMAN_LAB  <- c(C1="CD8 TEM (GZMK/GZMH)", C2="CD4 naive/CM", C3="CD8 TEM (ZEB2-hi)",
                C4="Th17", C5="Treg", C6="CD4 resting/CM (TSHZ2+)",
                C7="CD8 cytotoxic (GZMB/GNLY)", C8="MAIT", C9="CD4 Tph-like (CXCL13+)")

HUMAN_TARGET <- "C9"
MOUSE_TARGET <- "Tfh/Tph-like"

# Cell-cycle cluster: 1,326 markers untruncated. A proliferation program attracts
# spurious best-hits from every other cluster, so it is excluded from best-hit
# assignment (but kept in the matrices, flagged).
MOUSE_EXCLUDE_BESTHIT <- "Prolif T"

ANCHORS <- data.frame(
  human = c("BATF","STAT3","RORA","JUNB","ARAP2","DOCK8","ETV6",
            "INPP4B","LDLRAD4","NR3C1","PTPN13","RNF19A","TRPS1"),
  mouse = c("Batf","Stat3","Rora","Junb","Arap2","Dock8","Etv6",
            "Inpp4b","Ldlrad4","Nr3c1","Ptpn13","Rnf19a","Trps1"),
  stringsAsFactors = FALSE
)

cat("=====================================================================\n")
cat(" Module 06 / script 01 — marker overlap and specificity matrix\n")
cat("=====================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Load objects
# ---------------------------------------------------------------------------
cat("[1] Loading Seurat objects\n")
human <- readRDS(HUMAN_RDS)
mouse <- readRDS(MOUSE_RDS)

# --- identity gates: catch a wrong object silently substituted ---------------
stopifnot("cluster_code" %in% colnames(human@meta.data))
stopifnot("celltype_merged" %in% colnames(mouse@meta.data))
n_C9 <- sum(human$cluster_code == HUMAN_TARGET)
n_mt <- sum(mouse$celltype_merged == MOUSE_TARGET)
cat(sprintf("    human: %d cells, %d genes | C9 = %d cells\n",
            ncol(human), nrow(human), n_C9))
cat(sprintf("    mouse: %d cells, %d genes | %s = %d cells\n",
            ncol(mouse), nrow(mouse), MOUSE_TARGET, n_mt))
if (n_C9 != 873)    stop("IDENTITY GATE: human C9 should be 873 cells, got ", n_C9)
if (n_mt != 1045)   stop("IDENTITY GATE: mouse Tfh/Tph-like should be 1045 cells, got ", n_mt)
if (ncol(mouse) != 9224) stop("IDENTITY GATE: mouse object should be 9224 cells, got ", ncol(mouse))
cat("    identity gates passed\n\n")

# ---------------------------------------------------------------------------
# 2. Ortholog space
#
# Gene symbols differ between species. Uppercasing mouse symbols misses true
# orthologs with different names and invents pairs where names coincide by
# accident, so use a curated table instead:
#   babelgene::orthologs() -> support_n >= 3 (>=3 supporting databases)
#                          -> strict 1:1 both directions
#                          -> detected in >= MIN_CELLS cells in BOTH datasets
# ---------------------------------------------------------------------------
cat("[2] Building 1:1 ortholog space\n")
MIN_CELLS <- 3
h_det <- rownames(human)[Matrix::rowSums(GetAssayData(human, layer = "counts") > 0) >= MIN_CELLS]
m_det <- rownames(mouse)[Matrix::rowSums(GetAssayData(mouse, layer = "counts") > 0) >= MIN_CELLS]
h_det <- h_det[!grepl("^ENSG[0-9]+$", h_det)]   # babelgene rejects mixed symbol/Ensembl input

orth <- orthologs(genes = h_det, species = "mouse", human = TRUE) %>%
  filter(support_n >= 3) %>%
  group_by(human_symbol) %>% filter(n() == 1) %>%     # one mouse gene per human gene
  group_by(symbol)       %>% filter(n() == 1) %>%     # and vice versa
  ungroup() %>%
  filter(human_symbol %in% h_det, symbol %in% m_det) %>%
  select(human_symbol, mouse_symbol = symbol)

m2h <- setNames(orth$human_symbol, orth$mouse_symbol)  # mouse symbol -> human symbol
write.csv(orth, file.path(OBJ, "ortholog_map_1to1.csv"), row.names = FALSE)  # reused by script 03
cat(sprintf("    human detected: %d | mouse detected: %d\n", length(h_det), length(m_det)))
cat(sprintf("    strict 1:1 orthologs detected in both: %d\n", nrow(orth)))
cat(sprintf("    (naive uppercase symbol matching would give %d, for reference)\n\n",
            length(intersect(h_det, toupper(m_det)))))

# ---------------------------------------------------------------------------
# 3. Marker tables: as-published and recomputed at matched parameters
# ---------------------------------------------------------------------------
cat("[3] Marker tables\n")

# --- as published ------------------------------------------------------------
hm_pub <- read.csv(HUMAN_MARK, stringsAsFactors = FALSE) %>%
  filter(p_val_adj < PADJ, as.character(cluster) %in% names(HUMAN_CODE)) %>%
  mutate(grp = unname(HUMAN_CODE[as.character(cluster)]))

# Mouse published table is keyed by RNA_snn_res.0.8 (0-10). Map each numeric
# cluster to its dominant celltype_merged label so both versions are readable.
mouse_num_lab <- mouse@meta.data %>%
  count(RNA_snn_res.0.8, celltype_merged) %>%
  group_by(RNA_snn_res.0.8) %>% slice_max(n, n = 1) %>% ungroup() %>%
  transmute(cluster = as.character(RNA_snn_res.0.8),
            lab = paste0(celltype_merged, " (", cluster, ")"))
MOUSE_NUM_MAP <- setNames(mouse_num_lab$lab, mouse_num_lab$cluster)

mm_pub <- read.csv(MOUSE_MARK, stringsAsFactors = FALSE) %>%
  filter(p_val_adj < PADJ) %>%
  mutate(grp = ifelse(as.character(cluster) %in% names(MOUSE_NUM_MAP),
                      unname(MOUSE_NUM_MAP[as.character(cluster)]),
                      paste0("cluster ", cluster)))

cat("    published human list sizes (sig):\n"); print(table(hm_pub$grp))
cat("    published mouse list sizes (sig, top-500 TRUNCATED):\n"); print(table(mm_pub$grp))

# --- matched recomputation ---------------------------------------------------
matched_h_path <- file.path(OBJ, "markers_matched_human.csv")
matched_m_path <- file.path(OBJ, "markers_matched_mouse.csv")

if (file.exists(matched_h_path) && file.exists(matched_m_path)) {
  cat("    matched marker tables already on disk -- reusing\n")
  hm_mat <- read.csv(matched_h_path, stringsAsFactors = FALSE)
  mm_mat <- read.csv(matched_m_path, stringsAsFactors = FALSE)
} else {
  cat(sprintf("    recomputing both at min.pct=%.2f, logfc.threshold=%.2f, wilcox, untruncated\n",
              MATCH_MINPCT, MATCH_LFC))
  Idents(human) <- "cluster_code"
  hm_mat <- FindAllMarkers(human, only.pos = TRUE, test.use = "wilcox",
                           min.pct = MATCH_MINPCT, logfc.threshold = MATCH_LFC,
                           verbose = FALSE) %>%
    filter(p_val_adj < PADJ) %>% mutate(grp = as.character(cluster))
  write.csv(hm_mat, matched_h_path, row.names = FALSE)

  Idents(mouse) <- "celltype_merged"
  mm_mat <- FindAllMarkers(mouse, only.pos = TRUE, test.use = "wilcox",
                           min.pct = MATCH_MINPCT, logfc.threshold = MATCH_LFC,
                           verbose = FALSE) %>%
    filter(p_val_adj < PADJ) %>% mutate(grp = as.character(cluster))
  write.csv(mm_mat, matched_m_path, row.names = FALSE)
}
cat("    matched human list sizes:\n"); print(table(hm_mat$grp))
cat("    matched mouse list sizes:\n");  print(table(mm_mat$grp))
cat("\n")

# ---------------------------------------------------------------------------
# 4. Universes
# ---------------------------------------------------------------------------
cat("[4] Universes\n")

# "testable": genes FindAllMarkers could actually have returned -- expressed in
# >= min.pct of cells in at least one cluster in EACH dataset.
pct_by_group <- function(obj, group_col, genes) {
  cnt <- GetAssayData(obj, layer = "counts")[genes, , drop = FALSE]
  g   <- as.character(obj@meta.data[[group_col]])
  sapply(unique(g), function(k) Matrix::rowMeans(cnt[, g == k, drop = FALSE] > 0))
}
h_testable <- {
  p <- pct_by_group(human, "cluster_code", orth$human_symbol)
  orth$human_symbol[apply(p, 1, max) >= MATCH_MINPCT]
}
m_testable_mouse <- {
  p <- pct_by_group(mouse, "celltype_merged", orth$mouse_symbol)
  orth$mouse_symbol[apply(p, 1, max) >= MATCH_MINPCT]
}
m_testable <- unname(m2h[m_testable_mouse])

universes <- list(
  permissive   = orth$human_symbol,
  testable     = intersect(h_testable, m_testable),
  marker_space = NULL   # version-specific; filled per marker-table version
)
cat(sprintf("    permissive (all 1:1 orthologs detected in both) : %d\n", length(universes$permissive)))
cat(sprintf("    testable   (passes min.pct in >=1 cluster, both) : %d\n", length(universes$testable)))

# ---------------------------------------------------------------------------
# 5. Overlap matrices
# ---------------------------------------------------------------------------
overlap_matrix <- function(hmark, mmark, U) {
  hg_all <- split(hmark$gene, hmark$grp)
  mg_all <- split(mmark$gene, mmark$grp)
  N <- length(U)
  out <- expand.grid(human = names(hg_all), mouse = names(mg_all),
                     stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(out)), function(i) {
    hg <- intersect(unique(hg_all[[out$human[i]]]), U)
    mg <- intersect(unique(unname(na.omit(m2h[unique(mg_all[[out$mouse[i]]])]))), U)
    k  <- length(intersect(hg, mg))
    e  <- length(hg) * length(mg) / N
    data.frame(K = length(hg), n = length(mg), N = N, k = k,
               expected = e, fold = if (e > 0) k / e else NA_real_,
               p = phyper(k - 1, length(hg), N - length(hg), length(mg),
                          lower.tail = FALSE))
  })
  bind_cols(out, bind_rows(res)) %>%
    mutate(p_BH = p.adjust(p, "BH"), p_bonferroni = p.adjust(p, "bonferroni"))
}

cat("\n[5] Overlap matrices\n")
versions <- list(published = list(h = hm_pub, m = mm_pub),
                 matched   = list(h = hm_mat, m = mm_mat))
all_mats <- list(); uni_sizes <- list()

for (v in names(versions)) {
  hmark <- versions[[v]]$h; mmark <- versions[[v]]$m
  ms <- intersect(intersect(unique(hmark$gene), orth$human_symbol),
                  unique(unname(na.omit(m2h[unique(mmark$gene)]))))
  unis <- universes; unis$marker_space <- ms
  for (u in names(unis)) {
    U <- unis[[u]]
    m <- overlap_matrix(hmark, mmark, U) %>% mutate(version = v, universe = u)
    all_mats[[paste(v, u, sep = "_")]] <- m
    write.csv(m, file.path(OBJ, sprintf("overlap_matrix_%s_%s.csv", v, u)),
              row.names = FALSE)
    uni_sizes[[paste(v, u, sep = "_")]] <-
      data.frame(version = v, universe = u, N = length(U))
    tgt <- m %>% filter(human == HUMAN_TARGET, grepl("^Tfh/Tph-like", mouse))
    cat(sprintf("    %-9s / %-12s  N=%5d  |  C9 x mouse-Tph: k=%d exp=%.1f fold=%.2f p=%.3g\n",
                v, u, length(U), tgt$k, tgt$expected, tgt$fold, tgt$p))
  }
}
write.csv(bind_rows(uni_sizes), file.path(OBJ, "universe_sizes.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Reciprocal best hits + positive controls
#
# The two pairings that MUST work if marker overlap is a usable instrument:
#   C2 (CD4 naive/CM) -> mouse Naive T-like
#   C4 (Th17)         -> mouse Th17-like
# Whether they are recovered is recorded per version x universe. If they are
# not, the instrument is unreliable on this data and the Tph result -- positive
# or negative -- must be reported with that caveat.
# ---------------------------------------------------------------------------
cat("\n[6] Reciprocal best hits and positive controls\n")
POS_CTRL <- c(C2 = "Naive T-like", C4 = "Th17-like")

bh <- lapply(names(all_mats), function(nm) {
  m <- all_mats[[nm]] %>% filter(!grepl(MOUSE_EXCLUDE_BESTHIT, mouse, fixed = TRUE))
  best_h <- m %>% group_by(human) %>% slice_max(fold, n = 1, with_ties = FALSE) %>%
    ungroup() %>% transmute(version, universe, from = human, best_match = mouse,
                            fold = round(fold, 2), p = signif(p, 3),
                            direction = "human -> mouse")
  best_m <- m %>% group_by(mouse) %>% slice_max(fold, n = 1, with_ties = FALSE) %>%
    ungroup() %>% transmute(version, universe, from = mouse, best_match = human,
                            fold = round(fold, 2), p = signif(p, 3),
                            direction = "mouse -> human")
  bind_rows(best_h, best_m)
}) %>% bind_rows() %>%
  mutate(positive_control = ifelse(
    from %in% names(POS_CTRL) & direction == "human -> mouse",
    # grepl() does not vectorise over `pattern`, so pair each control with its
    # own expected match explicitly.
    mapply(function(f, b) {
      if (!f %in% names(POS_CTRL)) return(NA_character_)
      if (grepl(POS_CTRL[[f]], b, fixed = TRUE)) "RECOVERED" else "FAILED"
    }, from, best_match, USE.NAMES = FALSE),
    NA_character_))
write.csv(bh, file.path(OBJ, "reciprocal_best_hits.csv"), row.names = FALSE)

pc <- bh %>% filter(!is.na(positive_control)) %>%
  select(version, universe, control = from, best_match, positive_control)
cat("    positive-control recovery:\n"); print(as.data.frame(pc), row.names = FALSE)

cat("\n    where the Tph clusters actually point:\n")
tph_bh <- bh %>%
  filter((from == HUMAN_TARGET & direction == "human -> mouse") |
         (grepl("^Tfh/Tph-like", from) & direction == "mouse -> human"))
print(as.data.frame(tph_bh %>% select(version, universe, from, best_match, fold)),
      row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Anchor-gene audit
#
# Are the 13 "cross-species anchor" genes specific to the Tph clusters, or do
# they mark many subsets in both species? This is what decides whether the
# anchor module score reads as Tph identity or as generic T-cell activation.
# ---------------------------------------------------------------------------
cat("\n[7] Anchor-gene audit\n")
h_tgt <- hm_pub %>% filter(grp == HUMAN_TARGET) %>%
  select(gene, h_lfc = avg_log2FC, h_pct1 = pct.1, h_pct2 = pct.2, h_padj = p_val_adj)
m_tgt <- mm_pub %>% filter(grepl("^Tfh/Tph-like", grp)) %>%
  select(gene, m_lfc = avg_log2FC, m_pct1 = pct.1, m_pct2 = pct.2, m_padj = p_val_adj)

anchor_audit <- ANCHORS %>%
  left_join(h_tgt, by = c("human" = "gene")) %>%
  left_join(m_tgt, by = c("mouse" = "gene")) %>%
  mutate(
    n_human_clusters_marked = sapply(human, function(g) sum(hm_pub$gene == g)),
    n_mouse_clusters_marked = sapply(mouse, function(g) sum(mm_pub$gene == g)),
    human_target_exclusive  = n_human_clusters_marked == 1,
    across(where(is.numeric), ~ signif(.x, 4))
  ) %>%
  arrange(desc(h_lfc))
write.csv(anchor_audit, file.path(OBJ, "anchor_gene_specificity.csv"), row.names = FALSE)
print(as.data.frame(anchor_audit %>%
        select(human, h_lfc, h_pct1, h_pct2, m_lfc, m_pct1, m_pct2,
               n_human_clusters_marked, n_mouse_clusters_marked)), row.names = FALSE)
cat(sprintf("\n    anchors in human C9 marker list: %d / 13\n", sum(!is.na(anchor_audit$h_lfc))))
cat(sprintf("    anchors in mouse Tph marker list: %d / 13\n", sum(!is.na(anchor_audit$m_lfc))))
cat(sprintf("    anchors EXCLUSIVE to C9 among 9 human clusters: %d / 13  (%s)\n",
            sum(anchor_audit$human_target_exclusive),
            paste(anchor_audit$human[anchor_audit$human_target_exclusive], collapse = ", ")))

# ---------------------------------------------------------------------------
# 8. Figure — fold-enrichment heatmap
# ---------------------------------------------------------------------------
cat("\n[8] Figure\n")
plot_df <- bind_rows(all_mats) %>%
  mutate(universe = factor(universe, levels = c("permissive","testable","marker_space")),
         version  = factor(version,  levels = c("published","matched")),
         is_target = human == HUMAN_TARGET & grepl("^Tfh/Tph-like", mouse),
         # log2 fold, clipped so one extreme cell does not flatten the palette
         l2 = pmax(pmin(log2(pmax(fold, 1e-6)), 3), -3))

p <- ggplot(plot_df, aes(x = mouse, y = human, fill = l2)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = subset(plot_df, is_target),
            colour = "black", linewidth = 0.9, fill = NA) +
  geom_text(aes(label = k), size = 2.1, colour = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-3, 3),
                       name = expression(log[2]~fold)) +
  facet_grid(version ~ universe, scales = "free_x", space = "free_x") +
  labs(title = "Cross-species marker overlap — fold enrichment",
       subtitle = paste0("Cell labels = observed shared genes. Black outline = ",
                         "human C9 x mouse Tfh/Tph-like.\n",
                         "Fold depends almost entirely on the background universe N."),
       x = "mouse cluster", y = "human cluster") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank())

svglite(file.path(FIG, "overlap_fold_heatmap.svg"), width = 15, height = 7)
print(p); invisible(dev.off())
cat("    saved: ../figures/overlap_fold_heatmap.svg\n")

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_01.txt"))
cat("\nDone — script 01.\n")
