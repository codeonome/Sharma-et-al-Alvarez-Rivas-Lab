# =====================================================================
# Step 06 -- T-cell cleanup and final annotation (C1-C9)
#   1. remove T/plasma doublets (> IG_DOUBLET_CUTOFF % of counts from Ig genes)
#   2. remove res-0.4 cluster 9 (ribosomal-gene-dominated, no lineage markers;
#      validation in step 05)
#   3. merge res-0.4 cluster 10 into cluster 0 (over-clustering)
#   4. assign final labels and display codes C1-C9
# UMAP / PCA / Harmony are NOT recomputed: coordinates stay those of step 04.
# Input  : objects/04_tcell_reclustered.rds
# Output : objects/06_tcell_final.rds  (20,802 cells; cluster_id, cluster_code)
#          results/06_tcell_cleanup_annotation/*.csv
# =====================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr); library(tibble) })
source("config.R"); source("R/helpers.R")
STEP <- "06_tcell_cleanup_annotation"
d <- step_dirs(STEP)

tcell <- readRDS(OBJ$tcell_cluster)
stopifnot(TCELL_RES_COL %in% colnames(tcell@meta.data), TCELL_UMAP %in% Reductions(tcell),
          all(c(GM_LABEL, PD_LABEL) %in% unique(tcell$condition)))

# Work from an explicit column, not the active identity
tcell$res04 <- factor(as.character(tcell@meta.data[[TCELL_RES_COL]]), levels = EXPECTED_CLUSTERS)
stopifnot(!any(is.na(tcell$res04)))

# Donor = sample without the repeat-biopsy suffix (PD164b -> PD164)
tcell$donor <- sub(SAMPLE_SUFFIX_RE, "", tcell$sample)
donor_key <- tcell@meta.data %>% count(donor, sample, condition, name = "n_cells") %>%
  arrange(condition, donor, sample)
write.csv(donor_key, file.path(d$results, "sample_to_donor_map.csv"), row.names = FALSE)
log_step("Samples: ", n_distinct(tcell$sample), " | donors: ", n_distinct(tcell$donor))

n_start      <- ncol(tcell)
counts_start <- table(tcell$res04)

# ---------------------------------------------------------------------
# 1. Ig-high T/plasma doublets
# ---------------------------------------------------------------------
tcell <- PercentageFeatureSet(tcell, features = ig_genes_in(tcell), col.name = "percent.ig")
ig_high <- tcell$percent.ig > IG_DOUBLET_CUTOFF
ig_doublet_tbl <- tcell@meta.data[ig_high, ] %>%
  count(res04, condition, name = "n_doublets") %>%
  complete(res04, condition, fill = list(n_doublets = 0))
write.csv(ig_doublet_tbl, file.path(d$results, "ig_doublets_by_cluster.csv"), row.names = FALSE)
doublet_cells <- colnames(tcell)[ig_high]
tcell <- subset(tcell, cells = setdiff(colnames(tcell), doublet_cells))
log_step("Removed ", length(doublet_cells), " Ig-high doublets -> ", ncol(tcell), " cells")

# ---------------------------------------------------------------------
# 2. Drop cluster 9; 3. merge cluster 10 into 0
# ---------------------------------------------------------------------
n_dropped <- sum(tcell$res04 %in% DROP_CLUSTERS)
tcell <- subset(tcell, cells = colnames(tcell)[!tcell$res04 %in% DROP_CLUSTERS])
log_step("Removed cluster ", DROP_CLUSTERS, ": ", n_dropped, " cells -> ", ncol(tcell), " cells")

res_chr <- as.character(tcell$res04)
res_chr <- ifelse(res_chr %in% names(MERGE_MAP), unname(MERGE_MAP[res_chr]), res_chr)
tcell$cluster_res04_merged <- droplevels(factor(res_chr, levels = EXPECTED_CLUSTERS))

# NK / gamma-delta contamination check (hypothesis rejected; kept as a record)
nk_genes <- intersect(c("KLRF1", "NCAM1", "TRDC", "TRGC1", "TRGC2", "CD3D", "CD3E", "TRAC"), rownames(tcell))
nk_check <- fetch_expr(tcell, c(nk_genes, "cluster_res04_merged")) %>%
  group_by(cluster = cluster_res04_merged) %>%
  summarise(across(all_of(nk_genes), ~ round(mean(.x > 0), 3)), .groups = "drop")
write.csv(nk_check, file.path(d$results, "nk_check_by_cluster.csv"), row.names = FALSE)

# Before / after per cluster
counts_end <- table(tcell$cluster_res04_merged)
change_tbl <- tibble(cluster  = names(counts_start),
                     n_before = as.integer(counts_start),
                     n_after  = as.integer(counts_end[names(counts_start)])) %>%
  mutate(n_after = replace_na(n_after, 0L),
         delta   = n_after - n_before,
         note    = case_when(
           cluster %in% DROP_CLUSTERS     ~ "REMOVED - low-complexity ribosomal artifact",
           cluster %in% names(MERGE_MAP)  ~ paste("MERGED into cluster", MERGE_MAP[cluster]),
           cluster %in% unname(MERGE_MAP) ~ "gained merged cluster; lost Ig-high doublets",
           TRUE                           ~ "lost Ig-high doublets only"))
write.csv(change_tbl, file.path(d$results, "tcell_cleanup_changelog.csv"), row.names = FALSE)
log_step("Cells: ", n_start, " -> ", ncol(tcell))

# ---------------------------------------------------------------------
# 4. Final labels and display codes
# ---------------------------------------------------------------------
new_labels <- unname(LABEL_MAP)
stopifnot(setequal(names(LABEL_MAP), levels(tcell$cluster_res04_merged)))
# unname(): a character-indexed named vector would carry names into factor()
tcell$cluster_id <- factor(unname(LABEL_MAP[as.character(tcell$cluster_res04_merged)]), levels = new_labels)
code_map <- setNames(paste0("C", seq_along(new_labels)), new_labels)
tcell$cluster_code <- factor(unname(code_map[as.character(tcell$cluster_id)]), levels = unname(code_map))
stopifnot(!any(is.na(tcell$cluster_id)), !any(is.na(tcell$cluster_code)))

renaming_map <- tibble(
  res04_original = names(LABEL_MAP),
  merged_from    = ifelse(names(LABEL_MAP) %in% unname(MERGE_MAP),
                          paste(names(LABEL_MAP), "+", paste(names(MERGE_MAP), collapse = "+")),
                          names(LABEL_MAP)),
  display_code   = unname(code_map[new_labels]),
  final_label    = new_labels,
  n_cells        = as.integer(table(tcell$cluster_res04_merged)[names(LABEL_MAP)]))
write.csv(renaming_map, file.path(d$results, "cluster_renaming_map.csv"), row.names = FALSE)
print(as.data.frame(renaming_map))

saveRDS(tcell, OBJ$tcell_final)
save_session(STEP, d)
log_step("Done: ", OBJ$tcell_final, " (", ncol(tcell), " cells)")
