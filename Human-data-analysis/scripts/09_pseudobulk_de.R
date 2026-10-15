# =====================================================================
# Step 09 -- Pseudobulk differential expression, PD vs GM
#   All T cells pooled, and each cluster C1-C9 separately.
#   One profile per sample (summed raw counts); samples with fewer than
#   DE_MIN_CELLS_PER_SAMPLE cells in the cell set are dropped.
#   filterByExpr -> TMM -> edgeR::voomLmFit (donor as blocking factor,
#   sample quality weights) -> eBayes(robust = TRUE), contrast PD - GM.
# Input  : objects/06_tcell_final.rds
# Output : results/09_pseudobulk_de/DE_PD_vs_GM_<allTcells|C1..C9>.csv
#          DE_per_cluster_all_genes.csv, DE_per_cluster_summary.csv,
#          DE_highlight_genes_per_cluster.csv
#          figures/09_pseudobulk_de/DE_volcano_PD_vs_GM_*.svg, DE_voom_meanvar_allTcells.png
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})
source("config.R"); source("R/helpers.R")
STEP <- "09_pseudobulk_de"
d <- step_dirs(STEP)
stopifnot(requireNamespace("limma", quietly = TRUE), requireNamespace("edgeR", quietly = TRUE))

tcell <- readRDS(OBJ$tcell_final)
counts_all <- SeuratObject::GetAssayData(tcell, assay = DE_ASSAY, layer = "counts")
meta <- tcell@meta.data
new_labels <- unname(LABEL_MAP)
code_map   <- setNames(paste0("C", seq_along(new_labels)), new_labels)

# ---------------------------------------------------------------------
# Fit one cell set
# ---------------------------------------------------------------------
run_de <- function(cells, tag) {
  m <- meta[cells, ]
  keep_s <- names(which(table(m$sample) >= DE_MIN_CELLS_PER_SAMPLE))
  cells <- cells[m$sample %in% keep_s]
  m <- meta[cells, ]

  # Pseudobulk: genes x samples, by sparse indicator matrix product
  sv  <- factor(m$sample)
  ind <- Matrix::sparse.model.matrix(~ 0 + sv)
  colnames(ind) <- levels(sv)
  pb <- as.matrix(counts_all[, cells] %*% ind)

  sm <- m %>% distinct(sample, condition, donor) %>% mutate(sample = as.character(sample))
  sm <- sm[match(colnames(pb), sm$sample), ]
  stopifnot(identical(sm$sample, colnames(pb)))
  grp    <- factor(as.character(sm$condition), levels = c(GM_LABEL, PD_LABEL))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c("GM", "PD")

  dge  <- edgeR::DGEList(counts = pb, samples = sm, group = grp)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = DE_MIN_COUNT, min.total.count = DE_MIN_TOTAL_COUNT)
  dge  <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")

  # Diagnostic only: consensus correlation of plain-voom logCPM within donors
  rho_diag <- limma::duplicateCorrelation(limma::voom(dge, design)$E, design,
                                          block = sm$donor)$consensus.correlation

  if (tag == "allTcells") png(file.path(d$figures, "DE_voom_meanvar_allTcells.png"),
                              width = 1400, height = 1100, res = 180)
  fit <- edgeR::voomLmFit(dge, design, block = sm$donor, sample.weights = TRUE, plot = tag == "allTcells")
  if (tag == "allTcells") dev.off()

  contr <- limma::makeContrasts(contrasts = "PD-GM", levels = design)
  fit2  <- limma::eBayes(limma::contrasts.fit(fit, contr), robust = TRUE)
  res <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P", confint = TRUE) %>%
    rownames_to_column("gene") %>%
    mutate(is_ig        = grepl(IG_PATTERN, gene),
           is_ribo      = grepl("^RP[LS]", gene),
           is_mito      = grepl("^MT-", gene),
           is_highlight = gene %in% DE_HIGHLIGHT_GENES,
           direction    = case_when(
             adj.P.Val < DE_FDR_CUTOFF & logFC >  DE_LFC_CUTOFF ~ paste("up in", PD_LABEL),
             adj.P.Val < DE_FDR_CUTOFF & logFC < -DE_LFC_CUTOFF ~ paste("up in", GM_LABEL),
             adj.P.Val < DE_FDR_CUTOFF                          ~ "significant, small effect",
             TRUE                                               ~ "not significant"))
  if (tag == "allTcells") res$is_C9_marker <- res$gene %in% DE_C9_SIGNATURE
  write.csv(res, file.path(d$results, paste0("DE_PD_vs_GM_", tag, ".csv")), row.names = FALSE)

  plot_volcano(res, tag, rho = fit$correlation, n_samples = ncol(pb), n_donors = n_distinct(sm$donor))
  list(res = res, n_GM = sum(grp == GM_LABEL), n_PD = sum(grp == PD_LABEL),
       genes_tested = nrow(dge), rho_diag = rho_diag, rho_fit = fit$correlation)
}

# ---------------------------------------------------------------------
# Volcano: y = -log10(raw p), colour = FDR-based call; Ig genes excluded
# (ambient, donor-driven); ribosomal/mitochondrial genes not labelled
# ---------------------------------------------------------------------
plot_volcano <- function(res, tag, rho, n_samples, n_donors) {
  volc_df <- res %>% filter(!is_ig) %>% mutate(neg_log10_p = -log10(P.Value))
  lab_df  <- volc_df %>% filter(adj.P.Val < DE_FDR_CUTOFF, abs(logFC) > DE_LFC_CUTOFF, !is_ribo, !is_mito) %>%
    slice_head(n = DE_N_LABEL)
  de_cols <- setNames(c("#A5231C", "#1E7A32", "#8C8C8C", "#D3D3D3"),
                      c(paste("up in", PD_LABEL), paste("up in", GM_LABEL),
                        "significant, small effect", "not significant"))
  sig_p <- volc_df$P.Value[volc_df$adj.P.Val < DE_FDR_CUTOFF]
  title_set <- if (tag == "allTcells") "all T cells" else
    paste0(tag, " ", names(code_map)[code_map == tag])
  p <- ggplot(volc_df, aes(x = logFC, y = neg_log10_p)) +
    geom_point(aes(colour = direction), size = 1.1, alpha = 0.8) +
    geom_vline(xintercept = c(-DE_LFC_CUTOFF, DE_LFC_CUTOFF), linetype = "dashed",
               colour = "grey55", linewidth = 0.35) +
    { if (length(sig_p)) geom_hline(yintercept = -log10(max(sig_p)), linetype = "dashed",
                                    colour = "grey55", linewidth = 0.35) } +
    scale_colour_manual(values = de_cols, name = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)) +
    labs(x = expression(log[2]~"fold change (PD / GM)"), y = expression(-log[10]~"(p-value)"),
         title = paste0("Differential expression, PD vs GM (", title_set, ")"),
         subtitle = sprintf("Pseudobulk voomLmFit with sample weights, donor-blocked (ρ = %.2f)", rho),
         caption = paste0("One profile per sample (n = ", n_samples, " samples, ", n_donors,
                          " donors). TMM-normalised.\nDashed lines: absolute log2FC = ", DE_LFC_CUTOFF,
                          " and the p-value at which FDR = ", DE_FDR_CUTOFF,
                          ".\nImmunoglobulin genes excluded (ambient, donor-driven).\n",
                          "Colour reflects BH-adjusted significance; the y-axis shows unadjusted p."))
  if (length(sig_p)) {
    p <- p + annotate("text", x = min(volc_df$logFC), y = -log10(max(sig_p)), label = "FDR = 0.05",
                      hjust = 0, vjust = -0.4, size = 2.6, colour = "grey45")
  }
  # Significant highlight genes (Tph/Tfh programme): yellow points, boxed labels
  hl_df  <- volc_df %>% filter(is_highlight, adj.P.Val < DE_FDR_CUTOFF)
  lab_df <- lab_df %>% filter(!is_highlight)
  p <- p + geom_point(data = hl_df, shape = 21, fill = "#F2C94C", colour = "black", size = 2.2, stroke = 0.4)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (nrow(lab_df) > 0) {
      p <- p + ggrepel::geom_text_repel(data = lab_df, aes(label = gene), size = 2.9, fontface = "italic",
                                        max.overlaps = 30, box.padding = 0.35, min.segment.length = 0.1,
                                        segment.colour = "grey60", seed = 42)
    }
    if (nrow(hl_df) > 0) {
      p <- p + ggrepel::geom_label_repel(data = hl_df, aes(label = gene), size = 2.9, fontface = "bold.italic",
                                         label.size = 0.3, label.padding = 0.15, box.padding = 0.5,
                                         min.segment.length = 0, segment.colour = "grey30", seed = 42)
    }
  }
  f <- file.path(d$figures, paste0("DE_volcano_PD_vs_GM_", tag))
  ggsave(paste0(f, ".svg"), p, width = 8, height = 6.5)
}

# ---------------------------------------------------------------------
# Run: pooled, then per cluster
# ---------------------------------------------------------------------
pooled <- run_de(rownames(meta), "allTcells")
log_step("All T cells: ", sum(pooled$res$adj.P.Val < DE_FDR_CUTOFF), " genes at FDR < ", DE_FDR_CUTOFF)

per_cluster <- list()
for (cc in unname(code_map)) {
  per_cluster[[cc]] <- run_de(rownames(meta)[as.character(meta$cluster_code) == cc], cc)
  log_step(cc, ": ", sum(per_cluster[[cc]]$res$adj.P.Val < DE_FDR_CUTOFF), " genes at FDR < ", DE_FDR_CUTOFF)
}

all_genes <- bind_rows(lapply(names(per_cluster), function(cc)
  cbind(display_code = cc, per_cluster[[cc]]$res)))
write.csv(all_genes, file.path(d$results, "DE_per_cluster_all_genes.csv"), row.names = FALSE)

summary_tbl <- bind_rows(lapply(names(per_cluster), function(cc) {
  x <- per_cluster[[cc]]; sig <- x$res$adj.P.Val < DE_FDR_CUTOFF
  data.frame(display_code = cc, cluster = sprintf("%s (%s)", cc, names(code_map)[code_map == cc]),
             n_GM = x$n_GM, n_PD = x$n_PD, genes_tested = x$genes_tested,
             rho = round(x$rho_diag, 3), fit_path = "voomLmFit",
             n_sig = sum(sig), n_sig_up_PD = sum(sig & x$res$logFC > 0),
             n_sig_up_GM = sum(sig & x$res$logFC < 0))
})) %>% arrange(desc(n_sig))
write.csv(summary_tbl, file.path(d$results, "DE_per_cluster_summary.csv"), row.names = FALSE)

highlight <- all_genes %>%
  filter(gene %in% DE_HIGHLIGHT_GENES) %>%
  transmute(display_code, gene, logFC = round(logFC, 3), CI_low = round(CI.L, 3), CI_high = round(CI.R, 3),
            AveExpr = round(AveExpr, 2), p_value = signif(P.Value, 3), FDR = signif(adj.P.Val, 3)) %>%
  arrange(gene, display_code)
write.csv(highlight, file.path(d$results, "DE_highlight_genes_per_cluster.csv"), row.names = FALSE)

save_session(STEP, d)
log_step("Done")
