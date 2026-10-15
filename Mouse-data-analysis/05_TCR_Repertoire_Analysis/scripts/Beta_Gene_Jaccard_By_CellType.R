suppressPackageStartupMessages({
  library(dplyr)
})


# ==============================================================================
# I. LOAD ANALYSIS SET
# ==============================================================================

acell <- read.csv("../objects/BETA_GENE_ACELL.csv", stringsAsFactors = FALSE) %>%
  filter(ANL01FL == "Y")

celltypes <- sort(unique(acell$CELLTYP))
cat(sprintf("Analysis set: %d cells across %d cell types (%s)\n",
            nrow(acell), length(celltypes), paste(celltypes, collapse = ", ")))

# ==============================================================================
# II. PER-CELL-TYPE JACCARD (4.1 / 4.2)
# ==============================================================================
# Gene presence = prevalence > 0 (n >= 1) within a given CELLTYP x COND group.

table5 <- lapply(celltypes, function(ct) {
  sub <- acell %>% filter(CELLTYP == ct)
  rec_genes <- unique(sub$BETAGC[sub$COND == "REC"])
  rcu_genes <- unique(sub$BETAGC[sub$COND == "RCU"])
  inter <- intersect(rec_genes, rcu_genes)
  uni   <- union(rec_genes, rcu_genes)
  data.frame(
    CELLTYP      = ct,
    REC_n        = length(rec_genes),
    RCU_n        = length(rcu_genes),
    Intersection = length(inter),
    Union        = length(uni),
    Jaccard      = round(length(inter) / length(uni), 3)
  )
}) %>% bind_rows() %>% arrange(desc(Jaccard))

write.csv(table5, "../objects/TABLE5_jaccard_by_celltype.csv", row.names = FALSE)

cat("\nJaccard similarity by cell type (REC vs RCU gene presence):\n")
print(table5)
cat(sprintf(
  "\nHighest overlap: %s (J = %.3f). Lowest overlap: %s (J = %.3f).\n",
  table5$CELLTYP[1], table5$Jaccard[1],
  table5$CELLTYP[nrow(table5)], table5$Jaccard[nrow(table5)]
))
cat("Saved TABLE5: Jaccard similarity by cell type\n")

# ==============================================================================
# III. SHARED-GENE x CELL TYPE BREAKDOWN
# ==============================================================================
# Drill-down on the 6 genes globally classified "Shared" in
# Beta_Gene_Prevalence_Within_Condition.R (TABLE3): which cell type(s) does
# each one actually appear in, and with what REC/RCU cell counts?

table3 <- read.csv("../objects/TABLE3_top_enriched_shared_genes.csv", stringsAsFactors = FALSE)
shared_genes <- table3 %>% filter(class == "Shared") %>% pull(BETAGC)

cat(sprintf("\nShared genes (n = %d): %s\n", length(shared_genes), paste(shared_genes, collapse = ", ")))

shared_by_celltype <- acell %>%
  filter(BETAGC %in% shared_genes) %>%
  count(BETAGC, BETAGENE_C, CELLTYP, COND, name = "n") %>%
  arrange(BETAGC, CELLTYP, COND)

write.csv(shared_by_celltype, "../objects/TABLE_shared_genes_by_celltype.csv", row.names = FALSE)

n_celltypes_per_gene <- shared_by_celltype %>%
  group_by(BETAGC) %>%
  summarise(n_celltypes = n_distinct(CELLTYP), .groups = "drop")

cat("\nNumber of distinct cell types each shared gene appears in:\n")
print(n_celltypes_per_gene)
cat("Saved TABLE_shared_genes_by_celltype: shared genes x cell type x condition breakdown\n")

cat("\n=== Gene sharing within cell type analysis complete ===\n")
