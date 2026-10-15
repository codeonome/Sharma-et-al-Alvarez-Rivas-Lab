#' ---------------------------------------------------------------------------
#' Module 06 — Cross-species Tph: MetaNeighbor
#'
#' The arbiter. Scripts 01 and 02 both failed to support a conserved Tph cell
#' type, but each has a known weakness: marker-list overlap is sensitive to
#' thresholds and background choice, and module scoring only tests a curated
#' gene list. MetaNeighbor (Crow et al. 2018, Nat Commun) is purpose-built for
#' this question -- it asks, using the whole transcriptome over 1:1 orthologs,
#' how well cells of one cluster in one dataset can be identified by their
#' nearest neighbours in the other dataset.
#'
#' Output is an AUROC per cluster pair:
#'   0.5  = chance (no correspondence)
#'   >0.8 = good correspondence
#'   >0.9 = strong / replicable cell type
#'
#' Two modes are run:
#'   one_vs_rest  -- each pair against all other cells. Generous; a cluster can
#'                   score high simply by being broadly distinct.
#'   one_vs_best  -- each cluster against its single closest competitor. The
#'                   stringent mode, and the one Crow et al. recommend for
#'                   deciding whether a cell type genuinely replicates.
#'
#' Input:  ../objects/ortholog_map_1to1.csv           (from script 01)
#'         ../../Human/tcell_final_seurat_object.rds
#'         ../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds
#' Output: ../objects/metaneighbor_auroc_{one_vs_rest,one_vs_best}.csv
#'         ../objects/metaneighbor_top_hits.csv
#'         ../objects/metaneighbor_tph_summary.csv
#'         ../objects/metaneighbor_variable_genes.txt
#'         ../figures/metaneighbor_auroc_heatmap.svg
#'         ../figures/metaneighbor_one_vs_best.svg
#'         ../results/sessionInfo_03.txt
#'
#' Run from 06_CrossSpecies_Tph/scripts/:  Rscript 03_metaneighbor.R
#' ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(MetaNeighbor)
  library(dplyr); library(tidyr); library(ggplot2); library(svglite)
})

set.seed(1234)

ORTH_MAP  <- "../objects/ortholog_map_1to1.csv"
HUMAN_RDS <- "../../Human/tcell_final_seurat_object.rds"
MOUSE_RDS <- "../../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds"
OBJ <- "../objects"; FIG <- "../figures"; RES <- "../results"

HUMAN_TARGET <- "C9"; MOUSE_TARGET <- "Tfh/Tph-like"
# Biologically unambiguous pairings -- if MetaNeighbor cannot find these, it is
# not working on this data and the Tph result cannot be interpreted either way.
POS_CTRL <- list(c("C2", "Naive T-like"), c("C4", "Th17-like"))

cat("=====================================================================\n")
cat(" Module 06 / script 03 — MetaNeighbor cross-species cell-type replicability\n")
cat("=====================================================================\n\n")

# ---------------------------------------------------------------------------
# 1. Ortholog map (written by script 01)
# ---------------------------------------------------------------------------
if (!file.exists(ORTH_MAP)) stop("Run 01_marker_overlap_specificity.R first — ", ORTH_MAP, " missing")
orth <- read.csv(ORTH_MAP, stringsAsFactors = FALSE)
cat(sprintf("[1] 1:1 orthologs detected in both datasets: %d\n", nrow(orth)))

# ---------------------------------------------------------------------------
# 2. Build a single merged matrix on the shared ortholog feature space
#
# Mouse rows are renamed to their human ortholog symbol so both datasets share
# an identical, identically ordered feature space. Log-normalised `data` is used;
# MetaNeighbor ranks within each dataset, so separate normalisation is fine.
# ---------------------------------------------------------------------------
cat("[2] Building merged expression matrix\n")
human <- readRDS(HUMAN_RDS)
stopifnot(sum(human$cluster_code == HUMAN_TARGET) == 873)
h_mat <- GetAssayData(human, layer = "data")[orth$human_symbol, ]
h_meta <- data.frame(study_id = "human",
                     cell_type = as.character(human$cluster_code),
                     stringsAsFactors = FALSE)
rm(human); invisible(gc(verbose = FALSE))

mouse <- readRDS(MOUSE_RDS)
stopifnot(sum(mouse$celltype_merged == MOUSE_TARGET) == 1045, ncol(mouse) == 9224)
m_mat <- GetAssayData(mouse, layer = "data")[orth$mouse_symbol, ]
rownames(m_mat) <- orth$human_symbol          # common feature space
m_meta <- data.frame(study_id = "mouse",
                     cell_type = as.character(mouse$celltype_merged),
                     stringsAsFactors = FALSE)
rm(mouse); invisible(gc(verbose = FALSE))

stopifnot(identical(rownames(h_mat), rownames(m_mat)))
colnames(h_mat) <- paste0("h_", seq_len(ncol(h_mat)))   # guarantee unique barcodes
colnames(m_mat) <- paste0("m_", seq_len(ncol(m_mat)))

merged <- cbind(h_mat, m_mat)
meta   <- rbind(h_meta, m_meta)
rm(h_mat, m_mat); invisible(gc(verbose = FALSE))

sce <- SingleCellExperiment(assays = list(logcounts = merged),
                            colData = DataFrame(meta))
cat(sprintf("    merged: %d genes x %d cells (human %d, mouse %d)\n",
            nrow(sce), ncol(sce), sum(meta$study_id == "human"),
            sum(meta$study_id == "mouse")))
cat("    cell types per study:\n")
print(table(meta$study_id, meta$cell_type))

# ---------------------------------------------------------------------------
# 3. Variable genes — genes variable in BOTH studies
# ---------------------------------------------------------------------------
cat("\n[3] Selecting variable genes (variable in both studies)\n")
var_genes <- variableGenes(dat = sce, exp_labels = sce$study_id)
cat(sprintf("    %d variable genes\n", length(var_genes)))
writeLines(var_genes, file.path(OBJ, "metaneighbor_variable_genes.txt"))
if (length(var_genes) < 200)
  warning("Fewer than 200 variable genes — AUROCs will be unstable")

# ---------------------------------------------------------------------------
# 4. MetaNeighborUS
# ---------------------------------------------------------------------------
cat("\n[4] Running MetaNeighborUS\n")
cat("    one_vs_rest ...\n")
auroc_rest <- MetaNeighborUS(var_genes = var_genes, dat = sce,
                             study_id = sce$study_id, cell_type = sce$cell_type,
                             fast_version = TRUE)
cat("    one_vs_best ...\n")
auroc_best <- MetaNeighborUS(var_genes = var_genes, dat = sce,
                             study_id = sce$study_id, cell_type = sce$cell_type,
                             fast_version = TRUE, one_vs_best = TRUE,
                             symmetric_output = FALSE)

tidy_auroc <- function(m, mode) {
  as.data.frame(as.table(as.matrix(m))) %>%
    setNames(c("row", "col", "auroc")) %>%
    filter(!is.na(auroc)) %>%
    separate(row, c("row_study", "row_type"), sep = "\\|", extra = "merge") %>%
    separate(col, c("col_study", "col_type"), sep = "\\|", extra = "merge") %>%
    filter(row_study != col_study) %>%           # cross-species pairs only
    mutate(mode = mode, auroc = round(auroc, 3))
}
cross_rest <- tidy_auroc(auroc_rest, "one_vs_rest")
cross_best <- tidy_auroc(auroc_best, "one_vs_best")
write.csv(cross_rest, file.path(OBJ, "metaneighbor_auroc_one_vs_rest.csv"), row.names = FALSE)
write.csv(cross_best, file.path(OBJ, "metaneighbor_auroc_one_vs_best.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Top hits and the Tph verdict
# ---------------------------------------------------------------------------
cat("\n[5] Results\n")
th <- topHits(cell_NV = auroc_rest, dat = sce, study_id = sce$study_id,
              cell_type = sce$cell_type, threshold = 0.9)
write.csv(th, file.path(OBJ, "metaneighbor_top_hits.csv"), row.names = FALSE)
cat("    topHits (AUROC >= 0.9):\n")
if (nrow(th)) print(as.data.frame(th), row.names = FALSE) else cat("      none\n")

hk <- paste0("human|", HUMAN_TARGET); mk <- paste0("mouse|", MOUSE_TARGET)

# Human C9's ranking over every mouse cluster, and vice versa
h_row <- cross_rest %>% filter(row_study == "human", row_type == HUMAN_TARGET) %>%
  arrange(desc(auroc)) %>% select(mouse_cluster = col_type, auroc)
m_row <- cross_rest %>% filter(row_study == "mouse", row_type == MOUSE_TARGET) %>%
  arrange(desc(auroc)) %>% select(human_cluster = col_type, auroc)
cat("\n    human C9 vs every mouse cluster (one_vs_rest AUROC):\n")
print(as.data.frame(h_row), row.names = FALSE)
cat("\n    mouse Tfh/Tph-like vs every human cluster (one_vs_rest AUROC):\n")
print(as.data.frame(m_row), row.names = FALSE)

get_auroc <- function(d, rs, rt, ct) {
  v <- d$auroc[d$row_study == rs & d$row_type == rt & d$col_type == ct]
  if (length(v)) v[1] else NA_real_
}

cat("\n    POSITIVE CONTROLS (one_vs_rest AUROC, human -> mouse):\n")
ctrl_rows <- lapply(POS_CTRL, function(p) {
  a <- get_auroc(cross_rest, "human", p[1], p[2])
  rank_of <- cross_rest %>% filter(row_study == "human", row_type == p[1]) %>%
    arrange(desc(auroc)) %>% mutate(r = row_number()) %>% filter(col_type == p[2]) %>% pull(r)
  cat(sprintf("      %-3s -> %-14s AUROC = %.3f   (rank %d of %d)   %s\n",
              p[1], p[2], a, rank_of,
              sum(cross_rest$row_study == "human" & cross_rest$row_type == p[1]),
              ifelse(a >= 0.8 & rank_of == 1, "PASS", "check")))
  data.frame(pair = paste(p[1], p[2], sep = " <-> "), auroc = a, rank = rank_of)
}) %>% bind_rows()

tph_rest <- get_auroc(cross_rest, "human", HUMAN_TARGET, MOUSE_TARGET)
tph_best <- get_auroc(cross_best, "human", HUMAN_TARGET, MOUSE_TARGET)
tph_rank_h <- h_row %>% mutate(r = row_number()) %>% filter(mouse_cluster == MOUSE_TARGET) %>% pull(r)
tph_rank_m <- m_row %>% mutate(r = row_number()) %>% filter(human_cluster == HUMAN_TARGET) %>% pull(r)

cat(sprintf("\n    >>> C9 <-> mouse Tfh/Tph-like\n"))
cat(sprintf("          one_vs_rest AUROC = %.3f  (rank %d of %d from human side, %d of %d from mouse side)\n",
            tph_rest, tph_rank_h, nrow(h_row), tph_rank_m, nrow(m_row)))
cat(sprintf("          one_vs_best AUROC = %.3f\n", tph_best))
cat(sprintf("          reciprocal top hit: %s\n",
            ifelse(tph_rank_h == 1 && tph_rank_m == 1, "YES", "NO")))

summary_tbl <- data.frame(
  comparison = c("C9 <-> mouse Tfh/Tph-like", ctrl_rows$pair),
  auroc_one_vs_rest = c(tph_rest, ctrl_rows$auroc),
  rank_from_human = c(tph_rank_h, ctrl_rows$rank),
  auroc_one_vs_best = c(tph_best,
                        sapply(POS_CTRL, function(p) get_auroc(cross_best, "human", p[1], p[2]))),
  role = c("test", rep("positive control", length(POS_CTRL))))
write.csv(summary_tbl, file.path(OBJ, "metaneighbor_tph_summary.csv"), row.names = FALSE)
cat("\n    summary:\n"); print(summary_tbl, row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------
cat("\n[6] Figures\n")
no_grid <- theme_classic(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

hm_df <- cross_rest %>% filter(row_study == "human") %>%
  mutate(is_tph = row_type == HUMAN_TARGET & col_type == MOUSE_TARGET)
p1 <- ggplot(hm_df, aes(col_type, row_type, fill = auroc)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = subset(hm_df, is_tph), colour = "black", linewidth = 0.9, fill = NA) +
  geom_text(aes(label = sprintf("%.2f", auroc)), size = 2.3, colour = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0.5, limits = c(0, 1), name = "AUROC") +
  labs(title = "MetaNeighbor — cross-species cell-type replicability (one_vs_rest)",
       subtitle = "0.5 = chance, >0.8 good, >0.9 strong. Black outline = C9 x mouse Tfh/Tph-like.",
       x = "mouse cluster", y = "human cluster") + no_grid
svglite(file.path(FIG, "metaneighbor_auroc_heatmap.svg"), width = 9.5, height = 5)
print(p1); invisible(dev.off())
cat("    saved: metaneighbor_auroc_heatmap.svg\n")

hb_df <- cross_best %>% filter(row_study == "human") %>%
  mutate(is_tph = row_type == HUMAN_TARGET & col_type == MOUSE_TARGET)
p2 <- ggplot(hb_df, aes(col_type, row_type, fill = auroc)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = subset(hb_df, is_tph), colour = "black", linewidth = 0.9, fill = NA) +
  geom_text(aes(label = sprintf("%.2f", auroc)), size = 2.3, colour = "grey15") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0.5, limits = c(0, 1), name = "AUROC") +
  labs(title = "MetaNeighbor — one_vs_best (stringent)",
       subtitle = "Each cluster against its single closest competitor, not against all others.",
       x = "mouse cluster", y = "human cluster") + no_grid
svglite(file.path(FIG, "metaneighbor_one_vs_best.svg"), width = 9.5, height = 5)
print(p2); invisible(dev.off())
cat("    saved: metaneighbor_one_vs_best.svg\n")

writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_03.txt"))
cat("\nDone — script 03.\n")
