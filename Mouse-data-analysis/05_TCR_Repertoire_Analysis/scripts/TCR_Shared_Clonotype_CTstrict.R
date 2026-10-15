suppressPackageStartupMessages({
  library(dplyr)
})


# ==============================================================================
# I. LOAD PER-CELL TABLE (CD4 expanded, REC + RCU)
# ==============================================================================

cd4_exp <- read.csv("../objects/Clonotype_TABLE1_per_cell_CD4_expanded_REC_RCU.csv",
                     stringsAsFactors = FALSE)

cat(sprintf("Loaded: %d cells (CD4 expanded, REC + RCU)\n", nrow(cd4_exp)))

# Require a fully paired TCR — CTstrict must carry real alpha AND beta info
# (single-chain cells have "NA;NA" on one side and can't be "everything the same")
paired <- cd4_exp %>% filter(!grepl("NA;NA", CTstrict), CTstrict != "NA_NA")

cat(sprintf("Fully paired (alpha + beta + CDR3.nt all present): %d cells\n", nrow(paired)))

# ==============================================================================
# II. IDENTIFY CTstrict CLONOTYPES PRESENT IN BOTH REC AND RCU
# ==============================================================================
# CTstrict = TRAV.TRAJ.TRAC;CDR3nt_alpha _ TRBV.TRBD.TRBJ.TRBC;CDR3nt_beta
# i.e. identical alpha genes + alpha CDR3 nucleotide + beta genes + beta CDR3
# nucleotide -- the strictest possible clonotype identity (stricter than the
# gene-only alpha_tokens/beta_tokens grouping used in TABLE2/TABLE3).

overlap <- paired %>%
  group_by(CTstrict) %>%
  summarise(
    REC_present = any(condition == "REC"),
    RCU_present = any(condition == "RCU"),
    REC_n_cells = sum(condition == "REC"),
    RCU_n_cells = sum(condition == "RCU"),
    .groups = "drop"
  ) %>%
  filter(REC_present, RCU_present)

cat(sprintf("\nCTstrict clonotypes (identical alpha+beta+CDR3.nt) shared between REC and RCU: %d\n",
            nrow(overlap)))

# ==============================================================================
# III. PULL OUT ALL CELLS BELONGING TO THESE SHARED CLONOTYPES
# ==============================================================================

shared_cells <- paired %>%
  semi_join(overlap, by = "CTstrict") %>%
  arrange(CTstrict, condition)

cat(sprintf("Cells belonging to these shared clonotypes: %d\n", nrow(shared_cells)))
cat("Condition x celltype breakdown:\n")
print(table(shared_cells$condition, shared_cells$celltype))

write.csv(shared_cells, "../objects/Clonotype_TABLE12_shared_CTstrict_REC_RCU.csv", row.names = FALSE)
cat("\nSaved TABLE12: cells sharing an identical alpha+beta+CDR3.nt clonotype between REC and RCU\n")
