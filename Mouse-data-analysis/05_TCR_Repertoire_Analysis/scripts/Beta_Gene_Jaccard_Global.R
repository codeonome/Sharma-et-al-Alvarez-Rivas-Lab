suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(svglite)
})

set.seed(123)

# ==============================================================================
# I. LOAD ANALYSIS SET
# ==============================================================================

acell <- read.csv("../objects/BETA_GENE_ACELL.csv", stringsAsFactors = FALSE) %>%
  filter(ANL01FL == "Y")

rec_cells <- acell %>% filter(COND == "REC")
rcu_cells <- acell %>% filter(COND == "RCU")
n_rec <- nrow(rec_cells)
n_rcu <- nrow(rcu_cells)

cat(sprintf("Analysis set: %d cells (REC = %d, RCU = %d)\n", nrow(acell), n_rec, n_rcu))

# ==============================================================================
# II. PRIMARY ANALYSIS — GLOBAL JACCARD (pooled across cell types)
# ==============================================================================

rec_genes <- unique(rec_cells$BETAGC)
rcu_genes <- unique(rcu_cells$BETAGC)
inter     <- intersect(rec_genes, rcu_genes)
uni       <- union(rec_genes, rcu_genes)
jaccard_obs <- length(inter) / length(uni)

table2 <- data.frame(
  Measure = c("Unique genes in REC, n", "Unique genes in RCU, n",
              "Shared genes (intersection), n", "Total unique genes (union), n",
              "Jaccard similarity coefficient"),
  Value = c(length(rec_genes), length(rcu_genes), length(inter), length(uni),
            round(jaccard_obs, 3))
)

write.csv(table2, "../objects/TABLE2_jaccard_global.csv", row.names = FALSE)

cat(sprintf(
  "\nGlobal Jaccard (primary): REC genes = %d, RCU genes = %d, shared = %d, union = %d, Jaccard = %.3f\n",
  length(rec_genes), length(rcu_genes), length(inter), length(uni), jaccard_obs
))
cat("Saved TABLE2: global Jaccard similarity\n")

# ==============================================================================
# III. SENSITIVITY ANALYSIS — SUBSAMPLING REC TO MATCH RCU SAMPLE SIZE
# ==============================================================================

n_iter <- 1000
jaccard_sub <- numeric(n_iter)

for (i in seq_len(n_iter)) {
  rec_sub_cells <- rec_cells[sample(nrow(rec_cells), n_rcu, replace = FALSE), ]
  rec_sub_genes <- unique(rec_sub_cells$BETAGC)
  inter_sub <- intersect(rec_sub_genes, rcu_genes)
  uni_sub   <- union(rec_sub_genes, rcu_genes)
  jaccard_sub[i] <- length(inter_sub) / length(uni_sub)
}

sub_summary <- data.frame(
  n_iterations = n_iter,
  median       = round(median(jaccard_sub), 3),
  mean         = round(mean(jaccard_sub), 3),
  sd           = round(sd(jaccard_sub), 3),
  ci_lower_2.5 = round(quantile(jaccard_sub, 0.025), 3),
  ci_upper_97.5 = round(quantile(jaccard_sub, 0.975), 3),
  observed     = round(jaccard_obs, 3)
)

write.csv(data.frame(iteration = seq_len(n_iter), jaccard = jaccard_sub),
          "../objects/TABLE_A2_jaccard_subsampled.csv", row.names = FALSE)
write.csv(sub_summary, "../objects/TABLE_A2_jaccard_subsampled_summary.csv", row.names = FALSE)

cat(sprintf(
  "\nSensitivity analysis (%d iterations, REC subsampled to n=%d):\n  Median = %.3f, Mean = %.3f, SD = %.3f, 95%% interval = (%.3f, %.3f)\n  Observed primary Jaccard (%.3f) %s within this interval.\n",
  n_iter, n_rcu, sub_summary$median, sub_summary$mean, sub_summary$sd,
  sub_summary$ci_lower_2.5, sub_summary$ci_upper_97.5, jaccard_obs,
  ifelse(jaccard_obs >= sub_summary$ci_lower_2.5 & jaccard_obs <= sub_summary$ci_upper_97.5,
         "falls", "does NOT fall")
))
cat("Saved TABLE_A2: subsampled Jaccard distribution + summary\n")

# ==============================================================================
# IV. FIGURE 1 — SENSITIVITY ANALYSIS HISTOGRAM
# ==============================================================================
# Layout follows the analysis plan's Figure 1: grey histogram +
# density curve overlay, solid line = observed, dashed = mean, dotted = 95% CI bounds.

df_sub <- data.frame(jaccard = jaccard_sub)

p_hist <- ggplot(df_sub, aes(x = jaccard)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "grey70", color = "white", linewidth = 0.2) +
  geom_density(color = "black", linewidth = 0.8) +
  geom_vline(xintercept = jaccard_obs, linetype = "solid", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = sub_summary$mean, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = sub_summary$ci_lower_2.5, linetype = "dotted", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = sub_summary$ci_upper_97.5, linetype = "dotted", color = "black", linewidth = 0.6) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 9)) +
  labs(
    title   = "Sensitivity Analysis of Jaccard Similarity",
    subtitle = "Distribution of subsampled Jaccard coefficients (1,000 iterations)",
    x = "Jaccard similarity coefficient", y = "Density",
    caption = sprintf("Solid = observed (%.3f); dashed = mean (%.3f); dotted = 95%% interval (%.3f, %.3f).",
                       jaccard_obs, sub_summary$mean, sub_summary$ci_lower_2.5, sub_summary$ci_upper_97.5)
  )

svglite("../figures/beta_gene_jaccard_sensitivity.svg", width = 6.5, height = 5.5)
print(p_hist)
dev.off()
cat("Saved: ../figures/beta_gene_jaccard_sensitivity.svg\n")

cat("\n=== Global Jaccard similarity analysis complete ===\n")
