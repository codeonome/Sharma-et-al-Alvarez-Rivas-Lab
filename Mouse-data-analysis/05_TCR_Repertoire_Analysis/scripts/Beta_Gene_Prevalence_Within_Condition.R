suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(svglite)
})


# Condition palette — matches project-wide convention (Module 03/04)
pal_cond <- c("REC" = "#59A14F", "RCU" = "#E15759")
pal_class <- c(
  "REC Enriched" = pal_cond[["REC"]],
  "RCU Enriched" = pal_cond[["RCU"]],
  "Shared"       = "#4E79A7",
  "Rare"         = "#BAB0AC"
)

# ==============================================================================
# I. LOAD ANALYSIS SET
# ==============================================================================

acell <- read.csv("../objects/BETA_GENE_ACELL.csv", stringsAsFactors = FALSE) %>%
  filter(ANL01FL == "Y")

n_rec <- sum(acell$COND == "REC")
n_rcu <- sum(acell$COND == "RCU")
cat(sprintf("Analysis set: %d cells (REC = %d, RCU = %d)\n", nrow(acell), n_rec, n_rcu))

# ==============================================================================
# II. PER-GENE PREVALENCE WITHIN CONDITION
# ==============================================================================

gene_counts <- acell %>%
  count(BETAGC, BETAGENE_C, COND, name = "n") %>%
  pivot_wider(names_from = COND, values_from = n, values_fill = 0) %>%
  rename(REC_n = REC, RCU_n = RCU) %>%
  mutate(
    REC_pct = round(100 * REC_n / n_rec, 1),
    RCU_pct = round(100 * RCU_n / n_rcu, 1),
    delta     = round(REC_pct - RCU_pct, 1),
    abs_delta = abs(delta)
  )

# ==============================================================================
# III. CLASSIFICATION
# ==============================================================================
# REC Enriched : delta > 0.5
# RCU Enriched : delta < -0.5
# Shared       : |delta| <= 0.5 and observed >= 3 times in either condition
# Rare         : |delta| <= 0.5 and observed < 3 times in both conditions

gene_counts <- gene_counts %>%
  mutate(
    class = case_when(
      delta >  0.5                                  ~ "REC Enriched",
      delta < -0.5                                  ~ "RCU Enriched",
      REC_n >= 3 | RCU_n >= 3                        ~ "Shared",
      TRUE                                           ~ "Rare"
    ),
    class = factor(class, levels = c("REC Enriched", "RCU Enriched", "Shared", "Rare"))
  )

table_a3 <- gene_counts %>%
  arrange(desc(abs_delta), BETAGC) %>%
  mutate(rank = row_number()) %>%
  select(rank, BETAGC, BETAGENE_C, REC_n, REC_pct, RCU_n, RCU_pct, delta, abs_delta, class)

write.csv(table_a3, "../objects/TABLE_A3_gene_prevalence_within_condition.csv", row.names = FALSE)

n_genes <- nrow(table_a3)
class_summary <- table_a3 %>%
  count(class, name = "n") %>%
  mutate(pct = round(100 * n / n_genes, 1))

cat(sprintf("\nTotal unique genes: %d\n", n_genes))
cat("Classification summary:\n")
print(class_summary)
cat("Saved TABLE_A3: gene prevalence within condition\n")

# ==============================================================================
# IV. TABLE 3 — TOP 10 REC ENRICHED, TOP 10 RCU ENRICHED, ALL SHARED
# ==============================================================================

top_rec    <- table_a3 %>% filter(class == "REC Enriched") %>% arrange(desc(delta)) %>% slice_head(n = 10)
top_rcu    <- table_a3 %>% filter(class == "RCU Enriched") %>% arrange(delta)       %>% slice_head(n = 10)
all_shared <- table_a3 %>% filter(class == "Shared")

table3 <- bind_rows(top_rec, top_rcu, all_shared) %>%
  arrange(rank)

write.csv(table3, "../objects/TABLE3_top_enriched_shared_genes.csv", row.names = FALSE)
cat(sprintf("\nSaved TABLE3: top %d REC enriched, top %d RCU enriched, %d shared genes\n",
            nrow(top_rec), nrow(top_rcu), nrow(all_shared)))

# ==============================================================================
# V. FIGURE 2 — PREVALENCE SCATTER (REC % vs RCU %), ALL GENES
# ==============================================================================

p_scatter <- ggplot(table_a3, aes(x = REC_pct, y = RCU_pct, color = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = pal_class, name = "Classification", drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Prevalence of Beta Genes in Cells Exposed to REC vs RCU Conditions",
       x = "Beta Gene Prevalence in REC (%)", y = "Beta Gene Prevalence in RCU (%)")

svglite("../figures/beta_gene_prevalence_scatter.svg", width = 6.5, height = 5.5)
print(p_scatter)
dev.off()
cat("Saved: ../figures/beta_gene_prevalence_scatter.svg\n")

# ==============================================================================
# VI. FIGURE 3 — ZOOMED VIEW (<= 3.0% ON BOTH AXES)
# ==============================================================================

p_scatter_zoom <- p_scatter +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 3)) +
  labs(title = "Prevalence of Beta Genes in Cells Exposed to REC vs RCU Conditions",
       subtitle = "Population: Beta Gene Prevalence ≤ 3.0%")

svglite("../figures/beta_gene_prevalence_scatter_zoomed.svg", width = 6.5, height = 5.5)
print(p_scatter_zoom)
dev.off()
cat("Saved: ../figures/beta_gene_prevalence_scatter_zoomed.svg\n")

# ==============================================================================
# VII. FIGURE 4 — TOP REC/RCU ENRICHED + SHARED GENES (DIVERGING BAR)
# ==============================================================================

bar_df <- table3 %>%
  mutate(gene_label = paste0(BETAGC, " (", BETAGENE_C, ")")) %>%
  arrange(delta) %>%
  mutate(gene_label = factor(gene_label, levels = gene_label))

p_bar <- ggplot(bar_df, aes(x = delta, y = gene_label, fill = class)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = pal_class, name = "Classification", drop = FALSE) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 8)) +
  labs(title = "Top 10 REC Enriched, RCU Enriched, and Shared Beta Genes",
       x = "Prevalence Difference (REC % − RCU %)", y = NULL)

svglite("../figures/beta_gene_top_enriched_shared_bar.svg", width = 8, height = 8)
print(p_bar)
dev.off()
cat("Saved: ../figures/beta_gene_top_enriched_shared_bar.svg\n")

cat("\n=== Gene prevalence within condition analysis complete ===\n")
