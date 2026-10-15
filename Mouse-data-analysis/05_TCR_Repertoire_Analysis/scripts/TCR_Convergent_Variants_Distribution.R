suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
})


# ==============================================================================
# DISTRIBUTION OF THE 3 SHARED NUCLEOTIDE VARIANTS OF THE DOMINANT CONVERGENT
# CLONOTYPE (TRAV6N-7.TRAJ13 / TRBV5.TRBD1.TRBJ1-4, CDR3aa CALGANSGTYQRF /
# CASSQEGQGTNERLFF) — across the FULL T-cell population (all cell types, CD4+CD8).
#
# These are the 3 exact-nucleotide (CTstrict) variants that were found shared
# between REC and RCU in TABLE12. All 3 encode the identical receptor protein;
# they differ by 1-2 synonymous bases in the CDR3 nt seq. Stable IDs V1/V2/V3 are
# assigned HERE (not stored elsewhere) by their defining sequence.
# ==============================================================================

variants <- tibble::tribble(
  ~variant_id, ~CTstrict,
  "V1", "TRAV6N-7.TRAJ13.TRAC;TGTGCTCTGGGGGCAAATTCTGGGACTTACCAGAGGTTT_TRBV5.TRBD1.TRBJ1-4.TRBC1;TGTGCCAGCAGCCAAGAGGGACAGGGGACCAACGAAAGATTATTTTTC",
  "V2", "TRAV6N-7.TRAJ13.TRAC;TGTGCTCTGGGTGCAAATTCTGGGACTTACCAGAGGTTT_TRBV5.TRBD1.TRBJ1-4.TRBC1;TGTGCCAGCAGCCAAGAAGGACAGGGGACCAACGAAAGATTATTTTTC",
  "V3", "TRAV6N-7.TRAJ13.TRAC;TGTGCTCTGGGTGCAAATTCTGGGACTTACCAGAGGTTT_TRBV5.TRBD1.TRBJ1-4.TRBC1;TGTGCCAGCAGCCAAGAGGGACAGGGGACCAACGAAAGATTATTTTTC"
)

obj <- readRDS("../../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds")

meta <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  inner_join(variants, by = "CTstrict") %>%
  mutate(
    alpha_cdr3_nt = sub(".*;", "", sub("_.*", "", CTstrict)),
    beta_cdr3_nt  = sub(".*;", "", CTstrict)
  ) %>%
  select(variant_id, cell_id, condition, sample_id = orig.ident,
         celltype = celltype_merged, clone_size = clone_size_in_sample,
         expansion_status, tcr_status, alpha_cdr3_nt, beta_cdr3_nt,
         CTstrict, CTaa) %>%
  arrange(variant_id, condition, sample_id, celltype)

cat(sprintf("Total cells across the 3 variants (full T-cell population): %d\n", nrow(meta)))

# ---- complete per-cell table -------------------------------------------------
write.csv(meta, "../objects/Clonotype_TABLE14_convergent_variants_allcells.csv", row.names = FALSE)
cat("Saved TABLE14: per-cell table for V1/V2/V3 across all T cells\n\n")

# ---- summary: variant x condition x celltype counts --------------------------
summary_tab <- meta %>%
  count(variant_id, condition, celltype, name = "n_cells") %>%
  arrange(variant_id, condition, desc(n_cells))
write.csv(summary_tab, "../objects/Clonotype_TABLE14b_convergent_variants_summary.csv", row.names = FALSE)

cat("=== Cells per variant x condition (all cell types pooled) ===\n")
print(meta %>% count(variant_id, condition, name = "n_cells") %>%
        pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0))

cat("\n=== Cells per variant x cell type (all conditions pooled) ===\n")
print(meta %>% count(variant_id, celltype, name = "n_cells") %>%
        pivot_wider(names_from = celltype, values_from = n_cells, values_fill = 0))

cat("\n=== CD4 vs CD8 split per variant ===\n")
cd8 <- c("Teff CD8-like","Teff/Tem-like CD8","Trm/Exh-like CD8")
cd4 <- c("Tfh/Tph-like","Th17-like","Th2-like","Tregs-like")
print(meta %>%
        mutate(lineage = case_when(celltype %in% cd8 ~ "CD8",
                                   celltype %in% cd4 ~ "CD4",
                                   celltype == "Trm type mixed CD8+ CD4" ~ "mixed",
                                   celltype == "Prolif T" ~ "Prolif",
                                   TRUE ~ "other")) %>%
        count(variant_id, lineage, name = "n_cells") %>%
        pivot_wider(names_from = lineage, values_from = n_cells, values_fill = 0))

cat("\nSaved TABLE14b: variant x condition x celltype summary\n")
cat("\n=== Convergent variant distribution complete ===\n")
