suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})


obj <- readRDS("../../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds")

# Derive CD_broad (not stored in v2 object). 3-class: cluster 6
# ("Trm type mixed CD8+ CD4") is a confirmed CD4/CD8-mixed population (Module 03
# section 4.8c) → its own "mixed" class (excluded from CD4/CD8-specific analyses,
# retained in whole-population clonotype discovery).
cd8_types_all <- c("Teff CD8-like", "Trm/Exh-like CD8", "Teff/Tem-like CD8")
cd4_types_all <- c("Tfh/Tph-like", "Th17-like", "Th2-like", "Tregs-like")
obj$CD_broad <- dplyr::case_when(
  obj$celltype_merged %in% cd8_types_all           ~ "CD8",
  obj$celltype_merged %in% cd4_types_all           ~ "CD4",
  obj$celltype_merged == "Trm type mixed CD8+ CD4" ~ "mixed",
  TRUE                                             ~ NA_character_
)

# NOTE on "shared clonotypes": cells across conditions come from different mice.
# A clonotype appearing in both REC and RCU reflects convergent or highly
# abundant memory T cells — NOT longitudinal tracking of the same clone.

# ==============================================================================
# HELPERS: parse TCR chain gene tokens from CTgene
# ==============================================================================

# Extract TRA gene tokens (e.g. "TRAV6N-7;TRAJ13;TRAC")
extract_tra_tokens <- function(ctgene) {
  ctgene <- as.character(ctgene)
  ctgene[ctgene %in% c("", "NA")] <- NA
  alpha_raw <- str_split_fixed(ctgene, "_", 2)[, 1]
  tokens <- str_extract_all(alpha_raw, "TRA[A-Z0-9\\-]+")
  vapply(tokens, function(x)
    if (length(x) == 0) NA_character_ else paste(unique(x), collapse = ";"),
    character(1))
}

# Extract TRB gene tokens (e.g. "TRBV5;TRBD1;TRBJ2-4;TRBC1")
extract_trb_tokens <- function(ctgene) {
  ctgene <- as.character(ctgene)
  ctgene[ctgene %in% c("", "NA")] <- NA
  beta_raw <- str_split_fixed(ctgene, "_", 2)[, 2]
  tokens <- str_extract_all(beta_raw, "TRB[A-Z0-9\\-]+")
  vapply(tokens, function(x)
    if (length(x) == 0) NA_character_ else paste(unique(x), collapse = ";"),
    character(1))
}

# Extract raw beta gene string (for V–J–C level grouping, less stringent)
extract_beta_string <- function(ctgene) {
  ctgene <- as.character(ctgene)
  ctgene[ctgene %in% c("", "NA")] <- NA
  parts <- str_split_fixed(ctgene, "_", 2)
  beta  <- ifelse(parts[, 2] == "", NA, parts[, 2])
  beta  <- trimws(beta)
  beta  <- gsub("[^A-Za-z0-9+.-]", "", beta)
  beta[beta == ""] <- NA
  beta
}

# ==============================================================================
# I. BUILD PER-CELL TABLE — CD4 Expanded, REC + RCU
# ==============================================================================

md <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    alpha_tokens     = extract_tra_tokens(CTgene),
    beta_tokens      = extract_trb_tokens(CTgene),
    beta_gene_string = extract_beta_string(CTgene)
  )

cd4_exp <- md %>%
  filter(
    CD_broad          == "CD4",
    expansion_status  == "Expanded",
    condition         %in% c("REC", "RCU"),
    !is.na(beta_gene_string)
  ) %>%
  transmute(
    cell_id,
    condition,
    sample_id        = orig.ident,
    celltype         = celltype_merged,
    clone_size       = clone_size_in_sample,
    tcr_status,
    CTgene,
    CTstrict,
    CTaa,
    alpha_tokens,
    beta_tokens,
    beta_gene_string
  )

cat(sprintf("CD4 expanded cells (REC+RCU, β-parsed): %d\n", nrow(cd4_exp)))
cat("\nCells per condition × celltype:\n")
print(table(cd4_exp$condition, cd4_exp$celltype))

write.csv(cd4_exp, "../objects/Clonotype_TABLE1_per_cell_CD4_expanded_REC_RCU.csv",
          row.names = FALSE)
cat("\nSaved TABLE1: per-cell CD4 expanded REC/RCU\n")

# ==============================================================================
# II. CLONE OVERLAP: which clonotypes are shared between REC and RCU?
# ==============================================================================

clone_overlap <- cd4_exp %>%
  group_by(alpha_tokens, beta_tokens) %>%
  summarise(
    REC_present   = any(condition == "REC"),
    RCU_present   = any(condition == "RCU"),
    REC_n_cells   = sum(condition == "REC"),
    RCU_n_cells   = sum(condition == "RCU"),
    REC_celltypes = paste(sort(unique(celltype[condition == "REC"])), collapse = ";"),
    RCU_celltypes = paste(sort(unique(celltype[condition == "RCU"])), collapse = ";"),
    .groups       = "drop"
  ) %>%
  mutate(
    shared     = REC_present & RCU_present,
    rec_only   = REC_present & !RCU_present,
    rcu_only   = RCU_present & !REC_present
  ) %>%
  arrange(desc(shared), desc(REC_n_cells + RCU_n_cells))

cat(sprintf("\nClone overlap (by full αβ pair):\n  Shared REC+RCU: %d\n  REC only: %d\n  RCU only: %d\n",
            sum(clone_overlap$shared),
            sum(clone_overlap$rec_only),
            sum(clone_overlap$rcu_only)))

write.csv(clone_overlap, "../objects/Clonotype_TABLE2_clone_overlap_REC_RCU.csv",
          row.names = FALSE)
cat("Saved TABLE2: clone overlap REC vs RCU\n")

# Cells belonging to shared clonotypes
shared_clones <- clone_overlap %>% filter(shared)
shared_cells  <- cd4_exp %>%
  semi_join(shared_clones, by = c("alpha_tokens", "beta_tokens"))

cat(sprintf("\nCells from clones shared between REC and RCU: %d\n", nrow(shared_cells)))
cat("Celltype breakdown of shared clone cells:\n")
print(table(shared_cells$condition, shared_cells$celltype))

write.csv(shared_cells, "../objects/Clonotype_TABLE3_cells_in_shared_clones.csv",
          row.names = FALSE)
cat("Saved TABLE3: cells from shared clones\n")

# ==============================================================================
# III. BETA GENE COUNTS AND PROPORTIONS (condition × celltype)
# ==============================================================================

beta_counts <- cd4_exp %>%
  count(condition, celltype, beta_gene_string, name = "n_cells") %>%
  arrange(celltype, condition, desc(n_cells))

beta_props <- beta_counts %>%
  group_by(condition, celltype) %>%
  mutate(
    total_in_ct     = sum(n_cells),
    prop_within_ct  = n_cells / total_in_ct
  ) %>%
  ungroup()

write.csv(beta_counts, "../objects/Clonotype_TABLE4_beta_counts_condition_celltype.csv",
          row.names = FALSE)
write.csv(beta_props,  "../objects/Clonotype_TABLE5_beta_props_condition_celltype.csv",
          row.names = FALSE)
cat("Saved TABLE4: beta counts (condition × celltype)\n")
cat("Saved TABLE5: beta proportions (condition × celltype)\n")

# Wide format: which beta genes are shared/condition-specific
beta_wide <- beta_counts %>%
  select(condition, celltype, beta_gene_string, n_cells) %>%
  pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0) %>%
  mutate(
    shared   = (REC > 0 & RCU > 0),
    rec_only = (REC > 0 & RCU == 0),
    rcu_only = (RCU > 0 & REC == 0)
  )

write.csv(beta_wide, "../objects/Clonotype_TABLE6_beta_wide_REC_RCU_by_celltype.csv",
          row.names = FALSE)
cat("Saved TABLE6: beta genes wide (REC vs RCU, shared flag)\n")

# ==============================================================================
# IV. PER-SAMPLE BREAKDOWN (sample_id × condition × celltype)
# ==============================================================================

beta_counts_sample <- cd4_exp %>%
  count(sample_id, condition, celltype, beta_gene_string, name = "n_cells") %>%
  arrange(celltype, condition, sample_id, desc(n_cells))

beta_props_sample <- beta_counts_sample %>%
  group_by(sample_id, condition, celltype) %>%
  mutate(
    total_in_sample_ct     = sum(n_cells),
    prop_within_sample_ct  = n_cells / total_in_sample_ct
  ) %>%
  ungroup()

write.csv(beta_counts_sample, "../objects/Clonotype_TABLE7_beta_counts_per_sample.csv",
          row.names = FALSE)
write.csv(beta_props_sample,  "../objects/Clonotype_TABLE8_beta_props_per_sample.csv",
          row.names = FALSE)
cat("Saved TABLE7: beta counts per sample\n")
cat("Saved TABLE8: beta proportions per sample\n")

# ==============================================================================
# V. SUMMARY TABLES
# ==============================================================================

# Cells per sample × condition × celltype
cells_summary <- cd4_exp %>%
  count(sample_id, condition, celltype, name = "n_cells") %>%
  arrange(celltype, condition, sample_id)

write.csv(cells_summary, "../objects/Clonotype_TABLE9_cells_per_sample_condition_celltype.csv",
          row.names = FALSE)
cat("Saved TABLE9: cell counts per sample × condition × celltype\n")

# Beta gene diversity per sample × condition × celltype
diversity_summary <- cd4_exp %>%
  group_by(sample_id, condition, celltype) %>%
  summarise(
    n_cells      = n(),
    n_beta_genes = n_distinct(beta_gene_string),
    .groups      = "drop"
  ) %>%
  arrange(celltype, condition, sample_id)

write.csv(diversity_summary, "../objects/Clonotype_TABLE10_beta_diversity_per_sample.csv",
          row.names = FALSE)
cat("Saved TABLE10: beta gene diversity per sample\n")

# ==============================================================================
# VI. FULL ALPHA+BETA CHAIN TABLE (CD4 expanded, PD/REC/RCU, clone_size >= 5)
# ==============================================================================
# Adapted from module 05 section 7 — all three post-BAS conditions

cd4_all_conds <- md %>%
  filter(
    CD_broad         == "CD4",
    condition        %in% c("PD", "REC", "RCU"),
    clone_size_in_sample >= 5,
    !is.na(CTgene)
  ) %>%
  group_by(condition, orig.ident, celltype_merged,
           CTstrict, CTaa, alpha_tokens, beta_tokens, clone_size_in_sample) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(condition, desc(clone_size_in_sample))

write.csv(cd4_all_conds,
          "../objects/Clonotype_TABLE11_CD4_expanded_PD_REC_RCU_full_chains.csv",
          row.names = FALSE)
cat("Saved TABLE11: full αβ chain table (CD4 expanded, PD/REC/RCU, clone_size ≥ 5)\n")

cat("\n=== Module 05 — clonotype tracking complete ===\n")
cat("All tables written to: ../objects/\n")
