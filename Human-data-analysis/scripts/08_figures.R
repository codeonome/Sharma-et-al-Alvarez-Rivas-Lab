# =====================================================================
# Step 08 -- T-cell figures
#   UMAPs (C1-C9), density by condition, donuts, primary DA bar chart,
#   per-sample / per-donor proportions, grouped marker dot plot,
#   rho sensitivity.
# Input  : objects/06_tcell_final.rds, objects/07_da_results.rds
# Output : figures/08_figures/*.svg (DA_bar_blocked also .pdf)
#          results/08_figures/donut_pooled_percentages.csv, DA_bar_blocked_values.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(scales); library(viridisLite)
})
source("config.R"); source("R/helpers.R")
STEP <- "08_figures"
d <- step_dirs(STEP)
fig <- function(f) file.path(d$figures, f)

tcell <- readRDS(OBJ$tcell_final)
da    <- readRDS(OBJ$da)
new_labels <- unname(LABEL_MAP)
code_map   <- setNames(paste0("C", seq_along(new_labels)), new_labels)
pal        <- PALETTE
stopifnot(all(new_labels %in% names(pal)))

# ---------------------------------------------------------------------
# 1. UMAP: codes on the plot, colours keyed to the annotation
# ---------------------------------------------------------------------
umap_df <- as.data.frame(Embeddings(tcell, TCELL_UMAP))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cluster   <- tcell$cluster_id
umap_df$code      <- tcell$cluster_code
umap_df$condition <- tcell$condition
centroids <- umap_df %>% group_by(cluster, code) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_colour_manual(values = pal, guide = "none") +
  geom_text(data = centroids, aes(label = code), color = "black", size = 5, fontface = "bold") +
  theme_classic() + theme(axis.text = element_text(size = 8)) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "All subjects")
ggsave(fig("umap_tcell_final.svg"), p_umap, width = 7, height = 6.5)

legend_labels <- setNames(paste0(unname(code_map), " — ", new_labels), new_labels)
p_leg <- p_umap +
  scale_colour_manual(values = pal, labels = legend_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.text = element_text(size = 8))
ggsave(fig("umap_tcell_final_legend.svg"), p_leg, width = 9.5, height = 6.5)

p_umap_split <- p_umap + facet_wrap(~ condition, ncol = 2) + ggtitle("T-cell clusters by condition")
ggsave(fig("umap_tcell_final_by_condition.svg"), p_umap_split, width = 12, height = 6)

# ---------------------------------------------------------------------
# 2. UMAP density by condition (all cells as grey background)
# ---------------------------------------------------------------------
background_df <- umap_df[, c("UMAP_1", "UMAP_2")]
p_density <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(data = background_df, inherit.aes = FALSE, aes(x = UMAP_1, y = UMAP_2),
             color = "grey62", size = 0.35, alpha = 0.55) +
  stat_density_2d_filled(aes(fill = after_stat(level_mid), alpha = after_stat(level_mid)),
                         color = NA, bins = 10, adjust = 0.4, n = 200) +
  stat_density_2d(color = "black", linewidth = 0.1, alpha = 0.25, bins = 10, adjust = 0.4, n = 200) +
  scale_fill_gradientn(colours = viridisLite::turbo(256), name = "Density", na.value = NA) +
  scale_alpha_continuous(range = c(0, 0.85), guide = "none") +
  facet_wrap(~ condition, ncol = 2) +
  labs(title = "T-cell spatial density by condition", x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(face = "bold", size = 11), axis.text = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), legend.position = "right")
ggsave(fig("umap_density_by_condition.svg"), p_density, width = 11, height = 6, device = "svg")

# ---------------------------------------------------------------------
# 3. Donuts: % of all T cells in a condition (pooled cells), wedges in
#    fixed C1-C9 order so the two donuts are comparable
# ---------------------------------------------------------------------
make_donut <- function(cond_label, display_name) {
  dd <- tcell@meta.data %>%
    filter(condition == cond_label) %>%
    count(cluster_id, name = "n") %>%
    complete(cluster_id = factor(new_labels, levels = new_labels), fill = list(n = 0)) %>%
    mutate(pct = 100 * n / sum(n), label = factor(as.character(cluster_id), levels = new_labels),
           code = unname(code_map[as.character(cluster_id)])) %>%
    arrange(label)
  dd$display <- paste0(dd$code, "\n", sprintf("%.1f%%", dd$pct))
  dd <- dd %>% mutate(ypos = sum(pct) - (cumsum(pct) - pct / 2))   # geom_col stacks in reverse
  ggplot(dd, aes(x = 2, y = pct, fill = label)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y", start = 0) + xlim(0.5, 3.6) +
    scale_fill_manual(values = pal, guide = "none") +
    geom_text(aes(x = 3.02, y = ypos, label = display), size = 2.9, lineheight = 0.9, hjust = 0.5) +
    annotate("text", x = 0.5, y = 0, label = display_name, size = 6, fontface = "bold") +
    theme_void()
}
for (cond in list(c(GM_LABEL, "GM"), c(PD_LABEL, "PD"))) {
  p <- make_donut(cond[1], cond[2])
  ggsave(fig(paste0("donut_", cond[2], ".svg")), p, width = 6.5, height = 5.5)
}
donut_tbl <- tcell@meta.data %>%
  count(condition, cluster_id, name = "n_cells") %>%
  group_by(condition) %>% mutate(pct_pooled = round(100 * n_cells / sum(n_cells), 2)) %>% ungroup() %>%
  mutate(display_code = unname(code_map[as.character(cluster_id)])) %>%
  select(condition, display_code, cluster_id, n_cells, pct_pooled) %>%
  arrange(condition, display_code)
write.csv(donut_tbl, file.path(d$results, "donut_pooled_percentages.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 4. Primary DA figure: odds ratio (log axis) with 95% CI, p and FDR columns
# ---------------------------------------------------------------------
da_plot <- da$blocked_res %>%
  mutate(display_code = unname(code_map[cluster_id]),
         odds_ratio = exp(effect), or_low = exp(ci_low), or_high = exp(ci_high),
         direction = ifelse(effect > 0, "PD", "GM"),
         ci_excl_zero = (ci_low > 0) | (ci_high < 0),
         or_txt = sprintf("%.2f (%.2f–%.2f)", odds_ratio, or_low, or_high),
         p_txt = fmt_p(p_value), fdr_txt = fmt_p(FDR), y_lab = display_code) %>%
  arrange(effect) %>%
  mutate(ypos = row_number())
n_cl  <- nrow(da_plot)
e_rng <- range(c(da_plot$ci_low, da_plot$ci_high))
pad   <- 0.06 * diff(e_rng)
x_lo  <- min(0, e_rng[1] - pad); x_hi <- e_rng[2] + pad; span <- x_hi - x_lo
x_or  <- x_hi + 0.14 * span; x_p <- x_hi + 0.60 * span; x_fdr <- x_hi + 0.80 * span
or_breaks <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6)
keep_br   <- log(or_breaks) >= x_lo & log(or_breaks) <= x_hi

p_da <- ggplot(da_plot, aes(x = effect, y = ypos)) +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.5) +
  geom_col(aes(fill = direction, alpha = ci_excl_zero), width = 0.62, orientation = "y") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.20, linewidth = 0.45, colour = "grey20") +
  geom_text(aes(x = x_or, label = or_txt), hjust = 0.5, size = 3.0, colour = "grey15", family = "sans") +
  geom_text(aes(x = x_p, label = p_txt), hjust = 0.5, size = 3.0, colour = "grey15", family = "sans") +
  geom_text(aes(x = x_fdr, label = fdr_txt), fontface = ifelse(da_plot$FDR < 0.10, "bold", "plain"),
            hjust = 0.5, size = 3.0, colour = "grey15", family = "sans") +
  scale_fill_manual(values = c(PD = "#A5231C", GM = "#1E7A32"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.42), guide = "none") +
  scale_x_continuous(name = "odds ratio (PD vs GM), log scale", breaks = log(or_breaks[keep_br]),
                     labels = format(or_breaks[keep_br], trim = TRUE, drop0trailing = TRUE)) +
  scale_y_continuous(breaks = da_plot$ypos, labels = da_plot$y_lab, expand = expansion(add = c(0.6, 1.1))) +
  coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 10, colour = "black", hjust = 1),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 10, margin = margin(t = 6)),
        axis.line.y = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 12, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 9, colour = "grey35", margin = margin(b = 10)),
        plot.caption = element_text(size = 7.2, colour = "grey35", hjust = 0, margin = margin(t = 10)),
        plot.margin = margin(t = 8, r = 128, b = 8, l = 8)) +
  labs(y = NULL, title = "T-cell cluster abundance in periodontitis",
       subtitle = sprintf("propeller on logit proportions, donor-blocked (ρ = %.2f); n = 21 samples / 18 donors", da$rho),
       caption = paste0("Bars: odds ratio for each cluster in PD relative to healthy gingiva. ",
                        "Whiskers: unadjusted 95% moderated-t CI.\n",
                        "p = nominal; FDR = Benjamini–Hochberg across the 9 clusters. ",
                        "Faded bars have a CI spanning OR = 1.\n",
                        "Red = enriched in PD, green = enriched in GM.")) +
  geom_text(data = tibble(x = c(x_or, x_p, x_fdr), y = n_cl + 0.85, lab = c("OR (95% CI)", "p", "FDR")),
            aes(x = x, y = y, label = lab), fontface = c("bold", "bold.italic", "bold"),
            hjust = 0.5, size = 3.1, colour = "black", inherit.aes = FALSE)
for (ext in c("svg", "pdf")) {
  ggsave(fig(paste0("DA_bar_blocked.", ext)), p_da, width = 9.2, height = 4.6, bg = "white", dpi = 300)
}
da_table <- da_plot %>%
  transmute(display_code, cluster_label = cluster_id, log_odds = round(effect, 3),
            odds_ratio = round(odds_ratio, 3), OR_ci_low = round(or_low, 3), OR_ci_high = round(or_high, 3),
            mean_prop_GM = round(mean_GM, 4), mean_prop_PD = round(mean_PD, 4),
            p_nominal = signif(p_value, 4), FDR_BH = signif(FDR, 4)) %>%
  arrange(p_nominal)
write.csv(da_table, file.path(d$results, "DA_bar_blocked_values.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 5. Per-sample and per-donor proportions (the data behind every test)
# ---------------------------------------------------------------------
prop_plot_df <- da$comp_sample %>%
  mutate(display_code = factor(unname(code_map[as.character(cluster_id)]), levels = unname(code_map)),
         condition = factor(condition, levels = c(GM_LABEL, PD_LABEL)))
prop_theme <- list(facet_wrap(~ display_code, scales = "free_y", nrow = 3),
                   scale_fill_manual(values = COND_COLOURS, guide = "none"),
                   theme_bw(base_size = 10),
                   theme(panel.grid.minor = element_blank(),
                         strip.background = element_rect(fill = "grey92", colour = NA),
                         strip.text = element_text(face = "bold")))
p_props <- ggplot(prop_plot_df, aes(x = condition, y = 100 * proportion, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85, shape = 21, colour = "grey20") + prop_theme +
  labs(x = NULL, y = "% of sample's T cells", title = "Per-sample cluster proportions",
       subtitle = "One point per biopsy (n = 21). Repeat biopsies not collapsed.")
ggsave(fig("per_sample_proportions.svg"), p_props, width = 9, height = 7)

p_props_donor <- prop_plot_df %>%
  group_by(donor, condition, display_code) %>% summarise(proportion = mean(proportion), .groups = "drop") %>%
  ggplot(aes(x = condition, y = 100 * proportion, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.85, shape = 21, colour = "grey20") + prop_theme +
  labs(x = NULL, y = "% of donor's T cells", title = "Per-donor cluster proportions",
       subtitle = "One point per donor (n = 18), repeat biopsies averaged unweighted.")
ggsave(fig("per_donor_proportions.svg"), p_props_donor, width = 9, height = 7)

# ---------------------------------------------------------------------
# 6. Grouped marker dot plot (CD8 lobe, then CD4)
# ---------------------------------------------------------------------
cluster_order <- paste0("C", c(1, 3, 7, 8, 2, 4, 5, 6, 9))
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
gene_map <- utils::stack(gene_groups) %>% rename(gene = values, gene_group = ind) %>%
  filter(gene %in% rownames(tcell))
expr <- fetch_expr(tcell, c(gene_map$gene, "cluster_code"))
names(expr)[names(expr) == "cluster_code"] <- "cluster"
dot_df <- expr %>%
  pivot_longer(all_of(gene_map$gene), names_to = "gene", values_to = "value") %>%
  group_by(cluster, gene) %>%
  summarise(pct_exp = 100 * mean(value > 0), avg_exp = mean(expm1(value)), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(avg_scaled = if (sd(log1p(avg_exp)) > 0) as.numeric(scale(log1p(avg_exp))) else 0) %>%
  ungroup() %>%
  left_join(gene_map, by = "gene") %>%
  mutate(cluster = factor(cluster, levels = cluster_order),
         gene = factor(gene, levels = rev(gene_map$gene)),
         gene_group = factor(gene_group, levels = names(gene_groups)))
p_dot <- ggplot(dot_df, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct_exp, colour = avg_scaled)) +
  scale_size_continuous(name = "% Expressing", limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), range = c(0, 5)) +
  scale_colour_distiller(name = "Avg Expression", palette = "YlOrRd", direction = 1, limits = c(-2, 2),
                         breaks = c(-2, -1, 0, 1, 2), oob = scales::squish) +
  facet_grid(gene_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(position = "right") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA), panel.spacing.y = unit(2, "pt"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left = element_text(angle = 90, size = 6.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y.right = element_text(face = "italic", size = 8),
        axis.title = element_blank(), legend.position = "bottom")
ggsave(fig("dotplot_grouped.svg"), p_dot, width = 5.5, height = 9.5)

# ---------------------------------------------------------------------
# 7. rho sensitivity
# ---------------------------------------------------------------------
p_rho <- ggplot(da$rho_sweep, aes(x = rho, y = p_value, group = display_code)) +
  geom_hline(yintercept = FDR_CUTOFF, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(aes(colour = display_code), linewidth = 0.7) +
  geom_point(aes(colour = display_code), size = 1.6) +
  geom_vline(xintercept = da$rho, linetype = "dotted", colour = "grey20") +
  annotate("text", x = da$rho, y = 1, label = sprintf("estimated ρ = %.2f", da$rho),
           hjust = 1.05, vjust = 1, size = 3, colour = "grey20") +
  scale_y_log10() + scale_colour_discrete(name = NULL) +
  theme_classic(base_size = 11) +
  labs(x = "assumed intra-donor correlation (ρ)", y = "raw p-value (log scale)",
       title = "Sensitivity of differential abundance to intra-donor correlation",
       subtitle = "ρ = 0: samples independent | ρ → 1: biopsies interchangeable",
       caption = paste0("ρ = 0.99 is a numerical boundary where the block-correlation matrix ",
                        "approaches singularity;\nit lies outside the range supported by the data ",
                        "(per-cluster ICC 0.51-0.88)."))
ggsave(fig("rho_sensitivity.svg"), p_rho, width = 8, height = 5.5)

save_session(STEP, d)
log_step("Done")
