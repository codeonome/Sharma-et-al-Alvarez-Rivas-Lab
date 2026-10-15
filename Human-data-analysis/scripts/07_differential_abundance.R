# =====================================================================
# Step 07 -- Composition and differential abundance, PD vs GM
#   The unit of replication is the sample or the donor, never the cell.
#   PRIMARY model: propeller (logit proportions), sample level, donor as a
#   blocking factor with the duplicateCorrelation consensus rho.
#   Supplementary: Wilcoxon (sample, donor), propeller (sample, donor),
#   leave-one-donor-out, variance components, rho sensitivity sweep.
# Input  : objects/06_tcell_final.rds
# Output : results/07_differential_abundance/*.csv, intra_donor_correlation.txt
#          objects/07_da_results.rds (inputs for step 08 figures)
# =====================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr); library(tibble) })
source("config.R"); source("R/helpers.R")
STEP <- "07_differential_abundance"
d <- step_dirs(STEP)
stopifnot(requireNamespace("speckle", quietly = TRUE), requireNamespace("limma", quietly = TRUE))
wcsv <- function(x, f) write.csv(x, file.path(d$results, f), row.names = FALSE)

tcell <- readRDS(OBJ$tcell_final)
new_labels <- unname(LABEL_MAP)
code_map   <- setNames(paste0("C", seq_along(new_labels)), new_labels)

# ---------------------------------------------------------------------
# 1. Composition: sample level, and donor level (unweighted mean of a
#    donor's biopsy proportions)
# ---------------------------------------------------------------------
comp_sample <- tcell@meta.data %>%
  count(sample, donor, condition, cluster_id, name = "n_cells") %>%
  complete(nesting(sample, donor, condition), cluster_id, fill = list(n_cells = 0)) %>%
  group_by(sample) %>%
  mutate(unit_total = sum(n_cells), proportion = n_cells / unit_total) %>%
  ungroup()
wcsv(comp_sample, "sample_x_cluster_long.csv")
wcsv(comp_sample %>% select(sample, condition, cluster_id, n_cells) %>%
       pivot_wider(names_from = cluster_id, values_from = n_cells, values_fill = 0) %>%
       arrange(condition, sample), "sample_x_cluster_counts.csv")
wcsv(comp_sample %>% select(sample, condition, cluster_id, proportion) %>%
       mutate(proportion = round(proportion, 4)) %>%
       pivot_wider(names_from = cluster_id, values_from = proportion, values_fill = 0) %>%
       arrange(condition, sample), "sample_x_cluster_proportions.csv")

comp_donor <- comp_sample %>%
  group_by(donor, condition, cluster_id) %>%
  summarise(proportion = mean(proportion), n_samples = n_distinct(sample),
            n_cells = sum(n_cells), unit_total = sum(unit_total), .groups = "drop")
wcsv(comp_donor, "donor_x_cluster_long.csv")
wcsv(comp_donor %>% select(donor, condition, cluster_id, proportion) %>%
       mutate(proportion = round(proportion, 4)) %>%
       pivot_wider(names_from = cluster_id, values_from = proportion, values_fill = 0) %>%
       arrange(condition, donor), "donor_x_cluster_proportions.csv")

# ---------------------------------------------------------------------
# 2. Tests
# ---------------------------------------------------------------------
da_list <- list()
wilcox_sample <- wilcox_da(comp_sample, "sample")
wcsv(wilcox_sample, "DA_wilcox_sample_level.csv")
da_list$wilcox_sample <- wilcox_sample %>% mutate(test_id = "wilcox_sample")

wilcox_donor <- wilcox_da(comp_donor, "donor")
wcsv(wilcox_donor, "DA_wilcox_donor_level.csv")
da_list$wilcox_donor <- wilcox_donor %>% mutate(test_id = "wilcox_donor")

prop_sample <- propeller_da(tcell$cluster_id, tcell$sample, tcell$condition, "sample")
wcsv(prop_sample, "DA_propeller_sample_level.csv")
da_list$propeller_sample <- prop_sample %>% mutate(test_id = "propeller_sample")

# Donor as the unit pools cells within a donor (cell-number weighted)
prop_donor <- propeller_da(tcell$cluster_id, tcell$donor, tcell$condition, "donor")
wcsv(prop_donor, "DA_propeller_donor_level.csv")
da_list$propeller_donor <- prop_donor %>% mutate(test_id = "propeller_donor")

# PRIMARY
blocked <- propeller_blocked(tcell$cluster_id, tcell$sample, tcell@meta.data)
wcsv(blocked$res, "DA_propeller_donor_blocked.csv")
writeLines(sprintf("consensus_intra_donor_correlation\t%.6f", blocked$rho),
           file.path(d$results, "intra_donor_correlation.txt"))
da_list$propeller_blocked_sample <- blocked$res %>% mutate(test_id = "propeller_blocked_sample")
log_step(sprintf("Primary model: consensus within-donor rho = %.3f", blocked$rho))
print(as.data.frame(blocked$res), digits = 3)

# All tests in one table, long and wide
da_all <- bind_rows(da_list) %>%
  mutate(display_code = unname(code_map[cluster_id])) %>%
  select(all_of(DA_COLS)) %>%
  arrange(test_id, p_value)
wcsv(da_all, "DA_all_tests_combined.csv")
wcsv(da_all %>% select(display_code, cluster_id, test_id, p_value, FDR) %>%
       pivot_wider(names_from = test_id, values_from = c(p_value, FDR), names_glue = "{test_id}_{.value}") %>%
       arrange(display_code), "DA_all_tests_wide.csv")

# ---------------------------------------------------------------------
# 3. Leave-one-donor-out (Wilcoxon): sample level drops each repeat-biopsy
#    donor; donor level drops every donor in turn
# ---------------------------------------------------------------------
multi_donors <- comp_sample %>% distinct(sample, donor) %>% count(donor, name = "n_samples") %>%
  filter(n_samples > 1) %>% pull(donor)
sens_list <- list()
for (dn in multi_donors) {
  sens_list[[dn]] <- wilcox_da(filter(comp_sample, donor != dn), "sample") %>%
    transmute(cluster_id, dropped_donor = dn, p_value_after_drop = p_value, FDR_after_drop = FDR)
}
loo_list <- list()
for (dn in unique(comp_donor$donor)) {
  loo_list[[dn]] <- wilcox_da(filter(comp_donor, donor != dn), "donor") %>%
    transmute(cluster_id, dropped_donor = dn, p_value_after_drop = p_value)
}
sens_all <- bind_rows(bind_rows(sens_list) %>% mutate(level = "sample"),
                      bind_rows(loo_list) %>% mutate(level = "donor", FDR_after_drop = NA_real_)) %>%
  left_join(da_all %>% filter(test_id == "wilcox_sample") %>% select(cluster_id, p_value_full = p_value),
            by = "cluster_id") %>%
  mutate(display_code = unname(code_map[cluster_id])) %>%
  arrange(level, display_code, dropped_donor)
wcsv(sens_all, "DA_leave_one_donor_out.csv")
wcsv(sens_all %>% group_by(level, display_code, cluster_id) %>%
       summarise(p_full = first(p_value_full), p_min = min(p_value_after_drop, na.rm = TRUE),
                 p_max = max(p_value_after_drop, na.rm = TRUE), .groups = "drop") %>%
       arrange(level, p_max), "DA_leave_one_donor_out_summary.csv")

# ---------------------------------------------------------------------
# 4. Variance components and rho sensitivity
#    rho rests on 3 within-donor degrees of freedom per cluster, so the
#    model is refitted across the full range of rho. rho = 0.99 is a
#    numerical boundary, retained deliberately.
# ---------------------------------------------------------------------
props  <- blocked$props
smeta  <- blocked$smeta
cluster_names <- rownames(props$TransformedProps)
RHO_GRID <- sort(unique(c(RHO_GRID_BASE, round(blocked$rho, 3))))
donor_sizes <- table(smeta$donor)

# Per-cluster within/between-donor variance on the logit scale. var_between
# pools both conditions, so these ICCs are supporting evidence only; the
# reported rho is the duplicateCorrelation estimate.
tp <- props_to_df(props$TransformedProps, smeta)
rep_donors <- names(donor_sizes)[donor_sizes > 1]
vc <- lapply(cluster_names, function(cl) {
  v <- tp[[cl]]
  ss_w <- 0; df_w <- 0
  for (dn in rep_donors) {
    x <- v[tp$donor == dn]
    ss_w <- ss_w + sum((x - mean(x))^2)
    df_w <- df_w + length(x) - 1
  }
  var_within  <- if (df_w > 0) ss_w / df_w else NA_real_
  var_between <- var(tapply(v, tp$donor, mean))
  tibble(cluster_id = cl, display_code = unname(code_map[cl]),
         sd_within = sqrt(var_within), sd_between = sqrt(var_between),
         var_within = var_within, var_between = var_between,
         ICC = var_between / (var_between + var_within))
}) %>% bind_rows() %>% arrange(display_code)
wcsv(vc, "intra_donor_variance_components.csv")

# Effective sample size k / (1 + (k - 1) rho), summed over donors
eff_n   <- function(rho) sum(donor_sizes / (1 + (donor_sizes - 1) * rho))
eff_tbl <- tibble(rho = RHO_GRID, effective_N = round(sapply(RHO_GRID, eff_n), 2), nominal_N = nrow(smeta))

fit_at_rho <- function(rho) {
  fit <- if (rho == 0) limma::lmFit(props$TransformedProps, blocked$design) else
    limma::lmFit(props$TransformedProps, blocked$design, block = smeta$donor, correlation = rho)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, blocked$contr), robust = TRUE)
  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none", confint = TRUE)
  tibble(rho = rho, cluster_id = rownames(tt), display_code = unname(code_map[rownames(tt)]),
         logFC = tt$logFC, ci_low = tt$CI.L, ci_high = tt$CI.R, t = tt$t,
         p_value = tt$P.Value, FDR = tt$adj.P.Val)
}
rho_sweep <- bind_rows(lapply(RHO_GRID, fit_at_rho)) %>%
  left_join(eff_tbl %>% select(rho, effective_N), by = "rho")
wcsv(rho_sweep, "DA_propeller_rho_sensitivity.csv")

rho_wide_p <- rho_sweep %>% select(display_code, cluster_id, rho, p_value) %>%
  pivot_wider(names_from = rho, values_from = p_value, names_prefix = "p_rho_") %>% arrange(display_code)
rho_wide_fdr <- rho_sweep %>% select(display_code, cluster_id, rho, FDR) %>%
  pivot_wider(names_from = rho, values_from = FDR, names_prefix = "FDR_rho_") %>% arrange(display_code)
wcsv(left_join(rho_wide_p, rho_wide_fdr, by = c("display_code", "cluster_id")),
     "DA_propeller_rho_sensitivity_wide.csv")
wcsv(rho_sweep %>% group_by(display_code, cluster_id) %>%
       summarise(p_min = min(p_value), p_max = max(p_value), fdr_min = min(FDR), fdr_max = max(FDR),
                 logFC_min = min(logFC), logFC_max = max(logFC), .groups = "drop") %>%
       arrange(p_max), "DA_propeller_rho_sensitivity_summary.csv")

saveRDS(list(comp_sample = comp_sample, blocked_res = blocked$res, rho = blocked$rho,
             rho_sweep = rho_sweep), OBJ$da)
save_session(STEP, d)
log_step("Done")
