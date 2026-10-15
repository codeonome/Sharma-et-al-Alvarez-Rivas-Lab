suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
})


# ==============================================================================
# CONVERGENT CLONOTYPE SUMMARY (REC vs RCU)
# ==============================================================================
# A "convergent clonotype" = one receptor AMINO-ACID identity (identical alpha
# genes + beta genes + CDR3 alpha aa + CDR3 beta aa) that is encoded by MORE THAN
# ONE distinct nucleotide sequence (CTstrict). Convergence = the same protein
# receptor arising from multiple independent V(D)J recombination events.
#
# For each convergent clonotype present in BOTH REC and RCU we report three
# resolution levels, so downstream text/figures are never ambiguous:
#   (a) total cells sharing the amino-acid identity, per condition
#   (b) number of distinct nucleotide variants encoding it, per condition
#   (c) dominant clone size = size of the single largest EXACT-nucleotide clone
#       (one CTstrict, within one sample), per condition
#
# Population: full annotated αβ T-cell object (all cell types, all samples),
# MAIT-cleaned. clone_size is defined per exact CTstrict within each orig.ident,
# matching Module 04's definition.

obj <- readRDS("../../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds")

meta <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  filter(condition %in% c("REC", "RCU"),
         !is.na(CTgene), !is.na(CTaa), !is.na(CTstrict),
         !grepl("NA;NA", CTstrict), CTstrict != "NA_NA",
         tcr_status == "Paired (TRA+TRB)")

cat(sprintf("Paired REC/RCU cells with full TCR info: %d\n", nrow(meta)))

# aa-level clonotype key = genes + CDR3 amino acid (both chains)
meta <- meta %>% mutate(aa_key = paste(CTgene, CTaa, sep = "|"))

# ---- (a) total cells per condition, and (b) # nt variants per condition -------
per_cond <- meta %>%
  group_by(aa_key, condition) %>%
  summarise(
    total_cells      = n(),
    n_nt_variants    = n_distinct(CTstrict),
    .groups = "drop"
  )

# ---- (c) dominant clone size per condition (largest single CTstrict in a sample)
dom <- meta %>%
  count(aa_key, condition, orig.ident, CTstrict, name = "clone") %>%
  group_by(aa_key, condition) %>%
  summarise(dominant_clone_size = max(clone), .groups = "drop")

# ---- assemble wide, keep only convergent (>1 nt variant somewhere) & in BOTH ---
summary_wide <- per_cond %>%
  left_join(dom, by = c("aa_key", "condition")) %>%
  pivot_wider(
    names_from  = condition,
    values_from = c(total_cells, n_nt_variants, dominant_clone_size),
    values_fill = 0
  ) %>%
  # present in both conditions
  filter(total_cells_REC > 0, total_cells_RCU > 0) %>%
  # convergent: encoded by >1 nucleotide variant overall
  filter((n_nt_variants_REC + n_nt_variants_RCU) > 1)

# ---- split aa_key back into readable gene + CDR3aa columns --------------------
summary_out <- summary_wide %>%
  separate(aa_key, into = c("CTgene", "CTaa"), sep = "\\|", remove = TRUE) %>%
  mutate(
    alpha_gene    = sub("_.*", "", CTgene),
    beta_gene     = sub(".*_", "", CTgene),
    cdr3_alpha_aa = sub("_.*", "", CTaa),
    cdr3_beta_aa  = sub(".*_", "", CTaa),
    total_nt_variants = n_nt_variants_REC + n_nt_variants_RCU
  ) %>%
  transmute(
    alpha_gene, beta_gene, cdr3_alpha_aa, cdr3_beta_aa,
    REC_total_cells      = total_cells_REC,
    RCU_total_cells      = total_cells_RCU,
    REC_nt_variants      = n_nt_variants_REC,
    RCU_nt_variants      = n_nt_variants_RCU,
    REC_dominant_clone   = dominant_clone_size_REC,
    RCU_dominant_clone   = dominant_clone_size_RCU
  ) %>%
  arrange(desc(REC_total_cells + RCU_total_cells))

write.csv(summary_out, "../objects/Clonotype_TABLE13_convergent_REC_RCU.csv", row.names = FALSE)

cat(sprintf("\nConvergent clonotypes present in BOTH REC and RCU (>1 nt variant): %d\n", nrow(summary_out)))
cat("\nTop rows:\n")
print(head(summary_out, 10))
cat("\nSaved TABLE13: convergent clonotype summary (REC vs RCU, 3 resolution levels)\n")
