# Module 01 — Global Analysis

## What this step does
Annotates the four major populations (T, B, myeloid, stromal) in the already QC'd, CCA-integrated scRNA-seq + TCR object, and draws the global figures. QC, integration and TCR integration happened upstream (see *Input provenance*).

## Input
| File | Description |
|---|---|
| `objects/rivas_3111_integrated_scrna_tcr.rds` | 51,572 cells from 8 samples (BAS_1/2, PD_1/2, REC_1/2, RCU_1/2), CCA-integrated, clustered (`cca_clusters_0.05`: 8 clusters), TCR metadata attached. Not in git — download from Zenodo. |

### Input provenance
The script that built this object is not part of the repository. The parameters below were read from the object's own Seurat command history (`obj@commands`):

| Step | Parameters |
|---|---|
| TCR | scRepertoire `combineExpression(cloneCall = "strict", chain = "both", proportion = TRUE, filterNA = FALSE)` |
| Normalization | `NormalizeData`, LogNormalize, scale factor 10,000 |
| Variable genes | `FindVariableFeatures`, vst, 2,000 genes |
| Scaling / PCA | `ScaleData` (no regression); `RunPCA`, 50 PCs |
| Integration | CCA across the 8 sample layers → `integrated.cca` (50 dims) |
| Graph / clusters | `FindNeighbors` on `integrated.cca` dims 1–30 (k = 20); Louvain; resolution 0.05 used for annotation |
| UMAP | `RunUMAP` on `integrated.cca` dims 1–30 (30 neighbours, min.dist 0.3, cosine) → `umap.cca` |

Cell Ranger version, reference genome and QC cutoffs are not recorded in the object; retained cells have ≥ ~500 UMIs and < 19% mitochondrial reads.

## Key Analysis Steps
1. Set default assay to RNA; inspect cluster structure
2. `FindAllMarkers` (Wilcoxon, logFC ≥ 0.25, min.pct ≥ 0.1) to identify top 100 markers per cluster
3. Manual cluster → cell type annotation based on canonical markers
4. Condition labels (BAS / PD / REC / RCU) derived from `orig.ident` prefix

## Cluster Annotation
| Clusters | Cell Type |
|---|---|
| 0, 2, 5, 6 | Myeloid |
| 1, 4 | T |
| 3 | B |
| 7 | Stromal |

## Key Parameters
- CCA resolution: 0.05
- FindAllMarkers: `logfc.threshold = 0.25`, `min.pct = 0.1`, `test.use = "wilcox"`
- UMAP: `reduction = "umap.cca"`

## Outputs
| File | Description |
|---|---|
| `objects/global_clustered_annotated.rds` | Annotated Seurat object → **passed to Modules 02 and 03** |
| `objects/cca0.05_top100_markers.csv` | Top 100 DE markers per cluster |
| `figures/Fig4_A_B_combined.svg` | **Fig. S2A–B**: cluster UMAP (A) + marker dotplot by cell type (B) |
| `figures/UMAP_major_types_by_condition.svg` | Annotated UMAP split by BAS / PD / REC / RCU |
| `figures/pie_celltypes_by_condition.svg` | Cell type proportions per condition (pooled cells of both samples) |
| `figures/UMAP_pies_by_condition_combined.svg` | **Fig. 4B**: split UMAP (top) + pies (bottom) |

## Color Palettes
- **Cluster UMAP (panel A):** `cluster_cols` — NRC-inspired 8-color palette
  - 0 `#E64B35` · 1 `#4DBBD5` · 2 `#00A087` · 3 `#3C5488` · 4 `#F39B7F` · 5 `#8491B4` · 6 `#91D1C2` · 7 `#7E6148`
- **Annotated UMAP + pie charts:** `type_cols` — Okabe-Ito (colorblind-safe)
  - T `#0072B2` · B `#E69F00` · Myeloid `#009E73` · Stromal `#CC79A7`

## Figure layout: Fig4_A_B_combined.svg (paper Fig. S2A–B; file name kept from an earlier draft)
- Dimensions: 14 × 12 in, 600 dpi, SVG only
- **Panel A:** Cluster UMAP (square aspect ratio), clusters 0–7
- **Panel B:** DotPlot with cell-type group labels on right side via `facet_grid`
  - Order (top → bottom): Stromal → B Cells → T Cells → Myeloid
- Width ratio: `plot_layout(widths = c(1, 0.60))`
