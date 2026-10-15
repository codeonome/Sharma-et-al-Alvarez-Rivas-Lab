suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


# ==============================================================================
# I. BUILD BETA_GENE_ACELL — standardized cell-level analysis dataset
# ==============================================================================

raw <- read.csv("../objects/Clonotype_TABLE1_per_cell_CD4_expanded_REC_RCU.csv",
                 stringsAsFactors = FALSE)

n_raw <- nrow(raw)

beta_gene_acell <- raw %>%
  transmute(
    CID      = cell_id,
    COND     = condition,
    CELLTYP  = celltype,
    BETAGENE = beta_gene_string
  ) %>%
  mutate(
    BETAGENE_C = trimws(BETAGENE),
    BETAGENE_C = ifelse(BETAGENE_C %in% c("", "NA"), NA_character_, BETAGENE_C),
    ANL01FL    = ifelse(!is.na(BETAGENE_C), "Y", "N")
  )

# BETAGC: standardized gene ID (B001, B002, ...) ranked by descending
# prevalence in the analysis set; ties broken alphabetically for reproducibility
gene_rank <- beta_gene_acell %>%
  filter(ANL01FL == "Y") %>%
  count(BETAGENE_C, name = "n_cells") %>%
  arrange(desc(n_cells), BETAGENE_C) %>%
  mutate(BETAGC = sprintf("B%03d", row_number())) %>%
  select(BETAGENE_C, BETAGC)

beta_gene_acell <- beta_gene_acell %>%
  left_join(gene_rank, by = "BETAGENE_C") %>%
  select(CID, COND, CELLTYP, BETAGENE, BETAGENE_C, BETAGC, ANL01FL)

write.csv(beta_gene_acell, "../objects/BETA_GENE_ACELL.csv", row.names = FALSE)
cat(sprintf("Saved BETA_GENE_ACELL: %d cells (raw)\n", n_raw))

# ==============================================================================
# II. ANALYSIS SET DEFINITION
# ==============================================================================

analysis_set <- beta_gene_acell %>% filter(ANL01FL == "Y")

n_analysis <- nrow(analysis_set)
n_missing  <- n_raw - n_analysis
pct_missing <- round(100 * n_missing / n_raw, 1)

analysis_set_summary <- data.frame(
  n_raw        = n_raw,
  n_missing    = n_missing,
  pct_missing  = pct_missing,
  n_analysis   = n_analysis
)

write.csv(analysis_set_summary, "../objects/Beta_Gene_AnalysisSet_Summary.csv",
          row.names = FALSE)

cat(sprintf(
  "Analysis set: %d raw cells, %d missing (%.1f%%), %d in analysis set (ANL01FL = Y)\n",
  n_raw, n_missing, pct_missing, n_analysis
))

# ==============================================================================
# III. DISTRIBUTION OF BETA GENES (Appendix Table A1)
# ==============================================================================

table_a1 <- analysis_set %>%
  count(BETAGC, BETAGENE_C, name = "n_cells") %>%
  arrange(desc(n_cells), BETAGC) %>%
  mutate(
    pct        = round(100 * n_cells / n_analysis, 1),
    cum_n      = cumsum(n_cells),
    cum_pct    = round(100 * cum_n / n_analysis, 1),
    ge5_flag   = ifelse(n_cells >= 5, "Y", "N")
  ) %>%
  select(BETAGC, BETAGENE_C, n_cells, pct, cum_n, cum_pct, ge5_flag)

write.csv(table_a1, "../objects/TABLE_A1_beta_gene_distribution.csv", row.names = FALSE)

n_unique_genes <- nrow(table_a1)
n_ge5_genes    <- sum(table_a1$ge5_flag == "Y")
n_ge5_cells    <- sum(table_a1$n_cells[table_a1$ge5_flag == "Y"])
pct_ge5_cells  <- round(100 * n_ge5_cells / n_analysis, 1)

cat(sprintf(
  "Beta gene distribution: %d unique genes (range %d-%d cells/gene); %d genes with count >= 5, covering %d cells (%.1f%%)\n",
  n_unique_genes, min(table_a1$n_cells), max(table_a1$n_cells),
  n_ge5_genes, n_ge5_cells, pct_ge5_cells
))
cat("Saved TABLE_A1: beta gene distribution\n")

# ==============================================================================
# IV. TABLE 1 — CELL TYPE BY CONDITION (analysis set)
# ==============================================================================

celltype_counts <- analysis_set %>%
  count(CELLTYP, COND, name = "n") %>%
  group_by(COND) %>%
  mutate(
    col_total = sum(n),
    pct       = round(100 * n / col_total, 1)
  ) %>%
  ungroup() %>%
  select(-col_total)

celltype_counts_all <- analysis_set %>%
  count(CELLTYP, name = "n") %>%
  mutate(
    COND = "All",
    pct  = round(100 * n / n_analysis, 1)
  )

table1_long <- bind_rows(celltype_counts, celltype_counts_all) %>%
  mutate(display = sprintf("%d (%.1f%%)", n, pct))

table1_wide <- table1_long %>%
  select(CELLTYP, COND, display) %>%
  pivot_wider(names_from = COND, values_from = display)

col_totals <- analysis_set %>%
  count(COND, name = "n") %>%
  bind_rows(data.frame(COND = "All", n = n_analysis)) %>%
  mutate(n = as.character(n)) %>%
  pivot_wider(names_from = COND, values_from = n) %>%
  mutate(CELLTYP = "Column Total")

table1_final <- bind_rows(table1_wide, col_totals)

write.csv(table1_final, "../objects/TABLE1_celltype_by_condition.csv", row.names = FALSE)
write.csv(table1_long, "../objects/TABLE1_celltype_by_condition_long.csv", row.names = FALSE)

cat("Saved TABLE1: cell type by condition (analysis set)\n")
print(table1_final)

cat("\n=== Beta gene prevalence analysis complete ===\n")
cat("All outputs written to: ../objects/\n")
