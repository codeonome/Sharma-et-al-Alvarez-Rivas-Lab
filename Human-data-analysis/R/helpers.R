# =====================================================================
# Shared helper functions. Sourced after config.R by every step script.
# =====================================================================

# --- Logging ------------------------------------------------------------
log_step <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# Write sessionInfo() for the step next to its log.
save_session <- function(step, dirs) {
  writeLines(capture.output(sessionInfo()),
             file.path(dirs$logs, paste0(step, "_sessionInfo.txt")))
}

# --- Plotting -----------------------------------------------------------
# UMAP coloured by a metadata column. `guide = "none"` is set on the scale
# itself: hiding the legend afterwards still builds the guide, which
# crashes with some ggplot2 versions.
plot_umap <- function(obj, reduction, group_col, title, file_stub,
                      facet = FALSE, label_centroids = TRUE,
                      point_size = 0.3, width = 9, height = 7, facet_ncol = NULL) {
  df <- as.data.frame(Embeddings(obj, reduction))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$group <- as.character(obj[[group_col]][, 1])

  if (facet) {
    p <- ggplot() +
      geom_point(data = df[, c("UMAP_1", "UMAP_2")], aes(UMAP_1, UMAP_2),
                 color = "grey85", size = point_size) +
      geom_point(data = df, aes(UMAP_1, UMAP_2, color = group), size = point_size) +
      facet_wrap(~group, ncol = facet_ncol) +
      scale_colour_discrete(guide = "none") +
      theme_void() + ggtitle(title)
  } else {
    p <- ggplot(df, aes(UMAP_1, UMAP_2, color = group)) +
      geom_point(size = point_size, alpha = 0.6) +
      scale_colour_discrete(guide = "none") +
      theme_void() + ggtitle(title)
    if (label_centroids) {
      centroids <- df %>% group_by(group) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
      p <- p + geom_text(data = centroids, aes(label = group), color = "black",
                         size = 4, fontface = "bold")
    }
  }
  ggsave(paste0(file_stub, ".svg"), p, width = width, height = height)
  invisible(p)
}

# UMAP feature panels, one facet per gene; saved as SVG.
plot_marker_panel <- function(obj, reduction, genes, title, file_stub) {
  df <- as.data.frame(Embeddings(obj, reduction))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df <- cbind(df, FetchData(obj, vars = genes))
  long <- df %>% pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr")
  p <- ggplot(long, aes(UMAP_1, UMAP_2, color = expr)) +
    geom_point(size = 0.3) +
    scale_colour_gradient(low = "grey85", high = "blue4", guide = "none") +
    facet_wrap(~gene, ncol = 2) +
    theme_void() + ggtitle(title)
  h <- 4.5 * ceiling(length(genes) / 2)
  ggsave(paste0(file_stub, ".svg"), p, width = 9, height = h)
  invisible(p)
}

# Violins of a gene panel across clusters, `highlight` clusters filled blue.
plot_violin_panel <- function(obj, genes, cluster_col, title, file_stub,
                              highlight = NULL, ncol = 2) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    message("None of the requested genes found: skipping ", file_stub)
    return(invisible(NULL))
  }
  df <- FetchData(obj, vars = c(genes, cluster_col))
  colnames(df)[ncol(df)] <- "cluster"
  df_long <- df %>% pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr")
  df_long$focus <- if (is.null(highlight)) "no" else ifelse(df_long$cluster %in% highlight, "yes", "no")
  p <- ggplot(df_long, aes(x = cluster, y = expr, fill = focus)) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_fill_manual(values = c(no = "grey70", yes = "steelblue"), guide = "none") +
    facet_wrap(~gene, ncol = ncol, scales = "free_y") +
    theme_classic() + ggtitle(title)
  h <- 3 * ceiling(length(genes) / ncol)
  ggsave(paste0(file_stub, ".svg"), p, width = 10, height = h)
  invisible(p)
}

# Elbow plot of a PCA reduction.
plot_elbow <- function(obj, title, file) {
  sd <- Stdev(obj, reduction = "pca")
  p <- ggplot(data.frame(PC = seq_along(sd), stdev = sd), aes(PC, stdev)) +
    geom_point(size = 1.5) + theme_classic() +
    labs(x = "PC", y = "Standard deviation", title = title)
  ggsave(file, p, width = 8, height = 6)
}

# --- Clustering ----------------------------------------------------------
# FindClusters names columns RNA_snn_res.<res> at every level of
# subclustering; rename them immediately to a level-specific prefix.
rename_res_cols <- function(obj, prefix) {
  idx <- grep("^RNA_snn_res\\.", colnames(obj@meta.data))
  colnames(obj@meta.data)[idx] <- sub("^RNA_snn_res\\.", prefix, colnames(obj@meta.data)[idx])
  obj
}
res_cols <- function(obj, prefix) {
  cols <- grep(paste0("^", gsub(".", "\\.", prefix, fixed = TRUE)), colnames(obj@meta.data), value = TRUE)
  cols[order(as.numeric(sub(prefix, "", cols, fixed = TRUE)))]
}

# FindAllMarkers + top-10 per cluster (FDR < 0.05, ranked by log2FC).
markers_and_top10 <- function(obj, markers_csv, top10_csv) {
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(mk, markers_csv, row.names = FALSE)
  top10 <- mk %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
  write.csv(top10, top10_csv, row.names = FALSE)
  invisible(mk)
}

# Fraction of cells expressing each gene, per cluster.
frac_expressing <- function(obj, genes, cluster_col) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)
  df <- FetchData(obj, vars = c(genes, cluster_col))
  colnames(df)[ncol(df)] <- "cluster"
  df %>% group_by(cluster) %>%
    summarise(across(all_of(genes), ~round(mean(. > 0), 3)), n_cells = n(), .groups = "drop")
}

fetch_expr <- function(obj, vars) {
  tryCatch(FetchData(obj, vars = vars, layer = "data"),
           error = function(e) FetchData(obj, vars = vars, slot = "data"))
}

# --- Immunoglobulin genes -------------------------------------------------
ig_genes_in <- function(obj) {
  g <- grep(IG_PATTERN, rownames(obj), value = TRUE)
  setdiff(g, grep(IG_EXCLUDE, g, value = TRUE))
}

# --- Differential abundance ------------------------------------------------
DA_COLS <- c("test_id", "method", "level", "cluster_id", "display_code",
             "n_GM", "n_PD", "mean_GM", "mean_PD", "median_GM", "median_PD",
             "effect_type", "effect", "ci_low", "ci_high",
             "statistic", "p_value", "FDR", "sig")

# Wilcoxon rank-sum per cluster on a long proportion table
# (columns cluster_id, condition, proportion; one row per unit x cluster).
wilcox_da <- function(prop_df, level_label) {
  # eps: half the smallest non-zero proportion, so log2FC stays finite
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

# speckle::propeller on cell-level vectors (unit = sample or donor).
propeller_da <- function(clusters, unit, group, level_label) {
  res <- speckle::propeller(clusters = clusters, sample = unit, group = group,
                            transform = PROP_TRANSFORM)
  cl <- if ("BaselineProp.clusters" %in% names(res)) as.character(res$BaselineProp.clusters) else rownames(res)
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

# Primary model: propeller logit proportions at sample level, donor as a
# blocking factor with the consensus within-donor correlation from
# limma::duplicateCorrelation, eBayes(robust = TRUE).
propeller_blocked <- function(clusters, sample, meta) {
  props <- speckle::getTransformedProps(clusters = clusters, sample = sample,
                                        transform = PROP_TRANSFORM)
  smeta <- meta %>% distinct(sample, condition, donor)
  smeta <- smeta[match(colnames(props$TransformedProps), smeta$sample), ]
  stopifnot(identical(as.character(smeta$sample), colnames(props$TransformedProps)))

  grp    <- factor(as.character(smeta$condition), levels = c(GM_LABEL, PD_LABEL))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c("GM", "PD")

  rho <- limma::duplicateCorrelation(props$TransformedProps, design = design,
                                     block = smeta$donor)$consensus.correlation
  fit   <- limma::lmFit(props$TransformedProps, design, block = smeta$donor, correlation = rho)
  contr <- limma::makeContrasts(contrasts = "PD-GM", levels = design)
  fit2  <- limma::eBayes(limma::contrasts.fit(fit, contr), robust = TRUE)
  tt    <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P", confint = TRUE)

  pmat  <- as.matrix(unclass(props$Proportions))
  mu_gm <- rowMeans(pmat[, smeta$condition == GM_LABEL, drop = FALSE])
  mu_pd <- rowMeans(pmat[, smeta$condition == PD_LABEL, drop = FALSE])
  cl <- rownames(tt)
  res <- tibble(
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

  list(res = res, rho = rho, props = props, smeta = smeta, design = design, contr = contr)
}

# getTransformedProps() returns a `table`; as.matrix(unclass()) keeps the
# clusters x samples shape (as.data.frame() would melt it to long format).
props_to_df <- function(tprops, smeta) {
  m  <- as.matrix(unclass(tprops))
  df <- as.data.frame(t(m), check.names = FALSE)
  df$donor     <- smeta$donor
  df$condition <- as.character(smeta$condition)
  stopifnot(nrow(df) == ncol(m), all(rownames(m) %in% colnames(df)))
  df
}

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = digits)))
}
