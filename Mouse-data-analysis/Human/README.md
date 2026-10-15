# Human — reanalysis of GSE164241 (input to Module 06)

Human gingival T cells from the public scRNA-seq dataset GSE164241 (Williams et al., 2021; healthy gingival mucosa vs periodontitis), used in Module 06 for the cross-species and B-helper comparisons. This folder holds the **inputs** Module 06 reads and the script that produced the final human object. The full human pipeline that produces these files is in [`../../Human-data-analysis`](../../Human-data-analysis) (steps 00–09, verified against the published results).

## Files
| File | Role | In git |
|---|---|---|
| `tcell_final_seurat_object.rds` | Final annotated human T-cell object: 20,802 cells, 21 samples from 18 donors (Control 9,113; PD 11,689), nine clusters C1–C9 in `cluster_code`, UMAP `umap.tcell`. **Input to Module 06.** | no — Zenodo |
| `tcell_subset_markers.csv` | `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)` on the human T-cell clusters (raw cluster IDs 0–8 ↔ C1–C9). Input to `06/01_marker_overlap_specificity.R`. | yes |
| `human_mouse_Tph_chronic_activation_gene_table.csv` | Curated 81-row human–mouse gene table (categories, specificity flags). Input to `06/02_crossspecies_module_score.R`. | yes |

## How the final object was made
By [`../../Human-data-analysis`](../../Human-data-analysis). Steps 00–04 produce the T-cell object clustered at resolution 0.4 (11 clusters, 21,899 cells). Step 06 then:
1. removes T cell–plasma cell doublets (> 10% of counts from immunoglobulin genes);
2. drops cluster 9 (low-complexity, ribosomal-gene-high artifact) and merges cluster 10 into cluster 0;
3. assigns the final labels and display codes:

| Code | Label |
|---|---|
| C1 | CD8 TEM (GZMK/GZMH) |
| C2 | CD4 naive/CM |
| C3 | CD8 TEM (ZEB2-hi) |
| C4 | Th17 |
| C5 | Treg |
| C6 | CD4 resting/CM (TSHZ2+) |
| C7 | CD8 cytotoxic (GZMB/GNLY) |
| C8 | MAIT |
| C9 | CD4 Tph-like (CXCL13+) |

Step 07 runs the condition differential-abundance analysis (propeller on logit proportions, donor-blocked).

Rerunning the human pipeline is **not** needed to reproduce Modules 01–06: `run_all.sh` uses `tcell_final_seurat_object.rds` directly.
