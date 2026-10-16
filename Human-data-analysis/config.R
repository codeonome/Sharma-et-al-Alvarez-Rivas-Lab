# =====================================================================
# Pipeline configuration -- every path and analysis parameter lives here.
# Sourced by every script in scripts/. Run scripts from the pipeline/
# directory:   Rscript scripts/02_global_clustering.R
#
# Paths can be overridden with environment variables (see README):
#   PIPELINE_OUT    output root                      (default: output)
#   RAW_10X_DIR     Cell Ranger outputs for step 00
# =====================================================================

# --- Paths ------------------------------------------------------------
OUT_ROOT    <- Sys.getenv("PIPELINE_OUT", "output")
RAW_10X_DIR <- Sys.getenv("RAW_10X_DIR", "data/cellranger")

OBJ_DIR <- file.path(OUT_ROOT, "objects")
step_dirs <- function(step) {
  d <- list(results = file.path(OUT_ROOT, "results", step),
            figures = file.path(OUT_ROOT, "figures", step),
            logs    = file.path(OUT_ROOT, "logs"))
  for (x in c(d, OBJ_DIR)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  d
}

# Seurat objects passed between steps
OBJ <- list(
  integrated     = file.path(OBJ_DIR, "00_integrated.rds"),
  global         = file.path(OBJ_DIR, "02_global_clustered.rds"),
  immune         = file.path(OBJ_DIR, "03_immune_subset.rds"),
  tcell          = file.path(OBJ_DIR, "03_tcell_subset.rds"),
  tcell_cluster  = file.path(OBJ_DIR, "04_tcell_reclustered.rds"),
  tcell_final    = file.path(OBJ_DIR, "06_tcell_final.rds"),
  da             = file.path(OBJ_DIR, "07_da_results.rds")
)

# --- Study design -----------------------------------------------------
GM_LABEL <- "Control"            # value of meta.data$condition for healthy gingiva
PD_LABEL <- "PD"                 # value for periodontitis
SAMPLE_SUFFIX_RE <- "[a-z]+$"    # repeat biopsy suffix: PD164b -> donor PD164

# --- Step 00: QC + integration (Tyler Chinsky) ------------------------
QC_MIN_CELLS_PER_GENE <- 3
QC_MIN_FEATURES       <- 200
QC_MAX_FEATURES       <- 6000
QC_MAX_PERCENT_MT     <- 20
QC_MIN_CELLS_SAMPLE   <- 50
N_VARIABLE_FEATURES   <- 2000
INT_NPCS              <- 50
INT_UMAP_DIMS_UNINT   <- 1:30
INT_CLUSTER_DIMS      <- 1:50
INT_CLUSTER_RES       <- 0.5

# --- Step 02: global clustering ---------------------------------------
GLOBAL_REDUCTION <- "integrated.harmony"
GLOBAL_UMAP      <- "umap.harmony"
GLOBAL_DIMS      <- 1:30
GLOBAL_RES_SWEEP <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5)
GLOBAL_RES_FINAL <- 0.1

# --- Step 03: immune compartment --------------------------------------
# Global res-0.1 clusters: 1 T/NK, 3 plasma, 4 myeloid, 7 mast, 8 B cells
IMMUNE_CLUSTERS  <- c("1", "3", "4", "7", "8")
IMMUNE_NPCS      <- 50
IMMUNE_DIMS      <- 1:20
IMMUNE_RES_SWEEP <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
IMMUNE_RES_FINAL <- 0.2
# Immune res-0.2 clusters: 0 CD4 T, 1 CD8 T
TCELL_CLUSTERS   <- c("0", "1")

# --- Step 04: T-cell reclustering -------------------------------------
TCELL_NPCS      <- 30
TCELL_DIMS      <- 1:15
TCELL_RES_SWEEP <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
TCELL_RES_FINAL <- 0.4

# --- Seeds (values used for the published run) -------------------------
SEED_GLOBAL     <- 1234   # set.seed() at the start of steps 00 and 04
SEED_PCA_UMAP   <- 42     # seed.use for RunPCA / RunUMAP
SEED_HARMONY_T  <- 1234   # set.seed() immediately before T-cell Harmony
SEED_UMAP_IMM   <- 42     # set.seed() before the immune Harmony UMAP
SEED_CLUSTERS   <- 0      # FindClusters random.seed (Seurat default)

# --- Step 06: T-cell cleanup and annotation ----------------------------
TCELL_RES_COL     <- paste0("tcell_res.", TCELL_RES_FINAL)
TCELL_UMAP        <- "umap.tcell"
EXPECTED_CLUSTERS <- as.character(0:10)
DROP_CLUSTERS     <- "9"                   # ribosomal-gene-dominated, no lineage markers
MERGE_MAP         <- c("10" = "0")         # names = merged cluster, values = target
IG_DOUBLET_CUTOFF <- 10                    # % of counts from Ig genes
IG_PATTERN        <- "^IGH[ADEGM]|^IGHG[1-4]|^IGHGP|^IGK[CV]|^IGL[CV]"
IG_EXCLUDE        <- "^IGHMBP|^IGLON"      # non-immunoglobulin genes the pattern catches

# Final labels (res-0.4 cluster -> annotation). Order = display order C1-C9.
LABEL_MAP <- c(
  "0" = "CD8 TEM (GZMK/GZMH)",
  "1" = "CD4 naive/CM",
  "2" = "CD8 TEM (ZEB2-hi)",
  "3" = "Th17",
  "4" = "Treg",
  "5" = "CD4 resting/CM (TSHZ2+)",
  "6" = "CD8 cytotoxic (GZMB/GNLY)",
  "7" = "MAIT",
  "8" = "CD4 Tph-like (CXCL13+)"
)
PALETTE <- c(
  "CD8 TEM (GZMK/GZMH)"       = "#8B5E83",
  "CD4 naive/CM"              = "#A99B8E",
  "CD8 TEM (ZEB2-hi)"         = "#E38A3B",
  "Th17"                      = "#E8C84A",
  "Treg"                      = "#C0392B",
  "CD4 resting/CM (TSHZ2+)"   = "#EFA3AC",
  "CD8 cytotoxic (GZMB/GNLY)" = "#3E7A4E",
  "MAIT"                      = "#3B5CA8",
  "CD4 Tph-like (CXCL13+)"    = "#6FA8B8"
)
COND_COLOURS <- setNames(c("#1E7A32", "#A5231C"), c(GM_LABEL, PD_LABEL))

# --- Step 07: differential abundance ----------------------------------
FDR_CUTOFF     <- 0.05
PROP_TRANSFORM <- "logit"
RHO_GRID_BASE  <- c(0, 0.25, 0.5, 0.85, 0.95, 0.99)   # 0.99 kept deliberately (see README)

# --- Step 09: pseudobulk differential expression ------------------------
DE_ASSAY                <- "RNA"
DE_MIN_CELLS_PER_SAMPLE <- 20      # samples below this, within a cluster, are dropped
DE_MIN_COUNT            <- 10      # filterByExpr
DE_MIN_TOTAL_COUNT      <- 15      # filterByExpr
DE_FDR_CUTOFF           <- 0.05
DE_LFC_CUTOFF           <- 1
DE_N_LABEL              <- 20
DE_HIGHLIGHT_GENES      <- c("TOX2", "BCL6", "CXCR5", "CXCL13", "PDCD1")
# Hand-curated Tph signature of res-0.4 cluster 8 (= C9), flagged in the pooled table
DE_C9_SIGNATURE         <- c("CXCL13", "TOX2", "PTPN13", "TOX", "PDE7B", "ICOS", "MAF")
