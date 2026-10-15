#' ---
#' title: "TCR Dropout Analysis"
#' author: "Apoorva Sharma"
#' ---
#'
#' Task 3 — does sequencing depth drive whether a cell's TCR is recovered, and is that
#' effect independent of condition (BAS/PD/REC/RCU)?
#'
#'   Population : contamination-removed T cells (pre-alpha/beta filter), so TCR-negative
#'                ("No TCR") cells are retained — the alpha/beta-filtered objects cannot be used.
#'   Outcome    : tcr_detected = tcr_status != "No TCR" (any recovered chain).
#'   Predictor  : log10(nCount_RNA) (GEX UMI depth); condition as covariate.
#'   Plots      : (P1) recovery rate vs depth per condition, (P2) per-cell logistic fit,
#'                (P3) depth by recovery status.
#'   Model      : glm(binomial) additive [headline] + glmer (1|sample) [robustness]
#'                + depth x condition interaction LR test.
#'
#' CAVEAT: GEX UMI depth and 10x VDJ are separate libraries; nCount_RNA is a proxy for
#' overall cell quality/capture, not the VDJ read depth itself. A positive depth->recovery
#' association reflects cell quality, not a mechanical "more GEX UMIs = more TCR reads".

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(svglite)
  library(lme4)
})

set.seed(42)

# TCR status classifier — verbatim from 03_Tcell_Analysis/scripts/Tcell_analysis_reclustering.R:112-126
classify_tcr_status <- function(meta) {
  use_cols <- intersect(c("CTgene","CTaa","CTstrict","chains","c_call","tcr_chains"), colnames(meta))
  txt <- if (length(use_cols) == 0) rep("", nrow(meta)) else
    apply(meta[, use_cols, drop = FALSE], 1, function(r) { r[is.na(r)] <- ""; paste(r, collapse = ";") })
  has_a  <- grepl("TRAV|TRAJ|TRAC|\\bTRA\\b", txt, TRUE)
  has_b  <- grepl("TRBV|TRBJ|TRBC|\\bTRB\\b", txt, TRUE)
  has_gd <- grepl("TRGV|TRGJ|\\bTRG\\b|TRDV|TRDJ|\\bTRD\\b", txt, TRUE)
  status <- ifelse(has_gd & !(has_a | has_b), "γδ (TRG/TRD)",
            ifelse(has_a & has_b,              "Paired (TRA+TRB)",
            ifelse(has_a != has_b,             "Single chain (TRA or TRB)",
                                               "No TCR")))
  factor(status, levels = c("Paired (TRA+TRB)", "Single chain (TRA or TRB)", "γδ (TRG/TRD)", "No TCR"))
}

cond_colors <- c("BAS" = "grey70", "PD" = "#C0392B", "REC" = "#2C3E8C", "RCU" = "#E67E22")
cond_levels <- c("BAS", "PD", "REC", "RCU")

# ==============================================================================
# I. LOAD + BUILD POPULATION (contamination-removed T cells)
# ==============================================================================
obj <- readRDS("../../03_Tcell_Analysis/objects/tcell_only_clustered_unfilter02.rds")
cat("Loaded (pre-filter T cells):", ncol(obj), "cells\n")

remove_clusters <- c("6", "7", "10", "12", "13", "14")
clu <- if ("pca_clusters_tcells_0.3" %in% colnames(obj@meta.data))
         as.character(obj$pca_clusters_tcells_0.3) else as.character(Idents(obj))
keep <- colnames(obj)[!clu %in% remove_clusters]
obj  <- subset(obj, cells = keep)
cat("After dropping contamination clusters (", paste(remove_clusters, collapse = ","), "):",
    ncol(obj), "T cells\n")

# ==============================================================================
# II. PER-CELL VARIABLES
# ==============================================================================
md <- obj@meta.data %>%
  mutate(
    condition    = factor(sub("_.*", "", orig.ident), levels = cond_levels),
    tcr_status   = classify_tcr_status(obj@meta.data),
    tcr_detected = as.integer(classify_tcr_status(obj@meta.data) != "No TCR"),
    depth        = nCount_RNA,
    log10_depth  = log10(nCount_RNA)
  )
cat(sprintf("\nOverall TCR detection: %d / %d recovered (%.1f%%); %d 'No TCR' (%.1f%%)\n",
            sum(md$tcr_detected), nrow(md), 100 * mean(md$tcr_detected),
            sum(md$tcr_detected == 0), 100 * mean(md$tcr_detected == 0)))
cat("\nDetection by condition:\n")
print(round(100 * tapply(md$tcr_detected, md$condition, mean), 1))

# ==============================================================================
# III. LOGISTIC REGRESSION
# ==============================================================================
# Headline: additive fixed-effects GLM
m_add <- glm(tcr_detected ~ log10_depth + condition, data = md, family = binomial)
cat("\n--- GLM: tcr_detected ~ log10_depth + condition ---\n"); print(summary(m_add))

# Robustness: mixed model with per-sample random intercept (guards pseudoreplication)
m_mix <- tryCatch(
  glmer(tcr_detected ~ log10_depth + condition + (1 | orig.ident), data = md, family = binomial),
  error = function(e) { message("glmer failed: ", conditionMessage(e)); NULL })
if (!is.null(m_mix)) { cat("\n--- GLMER: + (1|orig.ident) ---\n"); print(summary(m_mix)) }

# Interaction: does the depth effect differ by condition?
m_int   <- glm(tcr_detected ~ log10_depth * condition, data = md, family = binomial)
lr_int  <- anova(m_add, m_int, test = "Chisq")
lr_drop <- drop1(m_add, test = "Chisq")   # LR tests for depth and condition main effects

# ---- coefficient table (both models) ----------------------------------------
coef_glm <- as.data.frame(summary(m_add)$coefficients)
names(coef_glm) <- c("estimate_logOR", "SE", "z", "p")
ci_glm   <- suppressMessages(confint(m_add))       # profile CI
coef_glm <- coef_glm %>%
  tibble::rownames_to_column("term") %>%
  mutate(OR = exp(estimate_logOR),
         CI_low = exp(ci_glm[term, 1]), CI_high = exp(ci_glm[term, 2]),
         model = "glm")

coef_tbl <- coef_glm
if (!is.null(m_mix)) {
  cm <- as.data.frame(summary(m_mix)$coefficients)
  names(cm) <- c("estimate_logOR", "SE", "z", "p")
  ci_m <- confint(m_mix, method = "Wald")
  coef_mix <- cm %>%
    tibble::rownames_to_column("term") %>%
    mutate(OR = exp(estimate_logOR),
           CI_low = exp(ci_m[term, 1]), CI_high = exp(ci_m[term, 2]),
           model = "glmer")
  coef_tbl <- bind_rows(coef_glm, coef_mix)
}
coef_tbl <- coef_tbl %>%
  select(model, term, estimate_logOR, OR, CI_low, CI_high, SE, z, p) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 4)))
write.csv(coef_tbl, "../objects/tcr_dropout_coefficients.csv", row.names = FALSE)
cat("\nSaved: ../objects/tcr_dropout_coefficients.csv\n")

# depth effect in interpretable units (per 2x depth) from the GLM
depth_b <- coef(m_add)["log10_depth"]
cat(sprintf("Depth OR per 10x UMI: %.2f | per 2x UMI: %.2f\n",
            exp(depth_b), exp(depth_b * log10(2))))

# ---- LR-test table -----------------------------------------------------------
tests_tbl <- data.frame(
  test       = c("depth (log10 UMI)", "condition", "depth × condition interaction"),
  statistic  = c(lr_drop["log10_depth", "LRT"], lr_drop["condition", "LRT"], lr_int$Deviance[2]),
  df         = c(1, length(cond_levels) - 1, lr_int$Df[2]),
  p          = c(lr_drop["log10_depth", "Pr(>Chi)"], lr_drop["condition", "Pr(>Chi)"], lr_int$`Pr(>Chi)`[2])
) %>% mutate(statistic = signif(statistic, 4))
write.csv(tests_tbl, "../objects/tcr_dropout_model_tests.csv", row.names = FALSE)
cat("Saved: ../objects/tcr_dropout_model_tests.csv\n"); print(tests_tbl)

# ==============================================================================
# IV. SUMMARY TABLES
# ==============================================================================
det_by_cond <- md %>%
  group_by(condition) %>%
  summarise(n_cells = n(), n_recovered = sum(tcr_detected),
            detection_rate = round(100 * mean(tcr_detected), 2),
            median_depth        = median(depth),
            median_depth_recov  = median(depth[tcr_detected == 1]),
            median_depth_drop   = median(depth[tcr_detected == 0]),
            .groups = "drop")
write.csv(det_by_cond, "../objects/tcr_dropout_detection_by_condition.csv", row.names = FALSE)
cat("Saved: ../objects/tcr_dropout_detection_by_condition.csv\n"); print(det_by_cond)

# per depth-decile recovery (global deciles), for P1
md$depth_decile <- dplyr::ntile(md$depth, 10)
rate_decile <- md %>%
  group_by(condition, depth_decile) %>%
  summarise(n = n(), median_depth = median(depth),
            recovery_fraction = mean(tcr_detected), .groups = "drop")
write.csv(rate_decile, "../objects/tcr_dropout_rate_by_depth_decile.csv", row.names = FALSE)
cat("Saved: ../objects/tcr_dropout_rate_by_depth_decile.csv\n")

# ==============================================================================
# V. PLOTS
# ==============================================================================
base_thm <- theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

# fitted probabilities from the additive GLM, over a depth grid per condition
grid <- expand.grid(
  log10_depth = seq(min(md$log10_depth), max(md$log10_depth), length.out = 200),
  condition   = factor(cond_levels, levels = cond_levels))
grid$fit   <- predict(m_add, newdata = grid, type = "response")
grid$depth <- 10 ^ grid$log10_depth

# ---- P1: recovery rate vs depth per condition (decile points + GLM curve) ----
p1 <- ggplot(rate_decile, aes(median_depth, recovery_fraction, color = condition)) +
  geom_line(data = grid, aes(depth, fit, color = condition), linewidth = 0.8) +
  geom_point(aes(size = n), alpha = 0.9) +
  scale_x_log10() +
  scale_color_manual(values = cond_colors, name = "Condition") +
  scale_size_continuous(name = "cells/bin", range = c(1, 4)) +
  labs(title = "TCR recovery rate vs sequencing depth, per condition",
       subtitle = "Points = depth-decile recovery fraction; lines = additive logistic fit",
       x = "UMI count per cell (nCount_RNA, log10)", y = "TCR detection rate") +
  base_thm
svglite("../figures/tcr_dropout_rate_vs_depth.svg", width = 7, height = 5)
print(p1); dev.off()
cat("Saved: ../figures/tcr_dropout_rate_vs_depth.svg\n")

# ---- P2: per-cell logistic fit ----------------------------------------------
jitter_layer <- if (requireNamespace("ggrastr", quietly = TRUE)) {
  ggrastr::geom_jitter_rast(height = 0.06, width = 0, size = 0.15, alpha = 0.12, raster.dpi = 150)
} else {
  geom_jitter(height = 0.06, width = 0, size = 0.15, alpha = 0.12)
}
p2 <- ggplot(md, aes(depth, tcr_detected, color = condition)) +
  jitter_layer +
  geom_line(data = grid, aes(depth, fit, color = condition), linewidth = 0.9) +
  scale_x_log10() +
  scale_color_manual(values = cond_colors, name = "Condition") +
  labs(title = "Per-cell TCR recovery vs depth (logistic fit)",
       subtitle = "Each point a cell (0 = No TCR, 1 = recovered); lines = additive logistic fit",
       x = "UMI count per cell (nCount_RNA, log10)", y = "TCR detected (0/1)") +
  base_thm
svglite("../figures/tcr_dropout_percell_logit.svg", width = 7, height = 5)
print(p2); dev.off()
cat("Saved: ../figures/tcr_dropout_percell_logit.svg\n")

# ---- P3: depth by recovery status, per condition ----------------------------
md$status_lab <- factor(ifelse(md$tcr_detected == 1, "TCR recovered", "No TCR (dropout)"),
                        levels = c("TCR recovered", "No TCR (dropout)"))
p3 <- ggplot(md, aes(condition, log10_depth, fill = status_lab)) +
  geom_violin(scale = "width", alpha = 0.6, color = NA, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.8),
               alpha = 0.9, linewidth = 0.3) +
  scale_fill_manual(values = c("TCR recovered" = "#2C7FB8", "No TCR (dropout)" = "#D95F0E"),
                    name = NULL) +
  labs(title = "Sequencing depth by TCR recovery status, per condition",
       x = "Condition", y = "UMI count per cell (log10 nCount_RNA)") +
  base_thm
svglite("../figures/tcr_dropout_depth_by_status.svg", width = 7.5, height = 5)
print(p3); dev.off()
cat("Saved: ../figures/tcr_dropout_depth_by_status.svg\n")

cat("\n=== TCR dropout analysis complete ===\n")
