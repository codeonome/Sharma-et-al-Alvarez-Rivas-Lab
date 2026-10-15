# Module 03 — T Cell Analysis

## What this module does
Subsets T cells from the global object, removes contaminating non-T populations, classifies TCR chain pairing, and performs αβ T cell recluster with full annotation. Produces all T cell figures for publication.

Core scripts, run in sequence (`../run_all.sh 03` runs the whole module):
1. `Tcell_analysis_reclustering.R` — Stages 1–3 (contamination removal + TCR classification)
2. `Tcell_ab_analysis_reclustering.R` — Stage 4+ (recluster, annotation, refinement, figures). `start_from_clustered <- FALSE` (default) reclusters from `tcell_TCRab_only03.rds`; set it to `TRUE` to reload the saved `tcell_TCRab_clustered04a.rds` checkpoint instead.
3. `Tcell_remove_MAIT.R` — physically drops the 27 `celltype_merged == NA` MAIT/innate-like cells flagged in step 2, producing the object all downstream modules actually use
4. `Tcell_doublet_check.R` — scDblFinder doublet rate per final cluster (`objects/cl6_doublet_check.csv`), the evidence cited for the cluster-6 annotation

Supporting diagnostics (`Tcell_Trm_CD8_identity_check.R`, `Tcell_Trm_CD8_feature_umap.R`, `Tph_marker_violin.R`) and figure-assembly scripts (`assemble_*.R`) follow; see the sections below.

---

## Input
| File | Description |
|---|---|
| `../../01_Global_Analysis/objects/global_clustered_annotated.rds` | Output of Module 01 |

---

## Script 1 — `Tcell_analysis_reclustering.R` (Stages 1–3)

### Steps
1. Subset T cells (`major_celltype == "T"`)
2. Normalize → HVG (n=3000) → Scale (regress `percent.mt`) → PCA (50 PCs)
3. Cluster at res 0.3 with 20 PCs → `FindAllMarkers` → identify contamination
4. **Cluster 9 TCR pairing check** — tabulate before filtering (Th2 vs ILC2 clarification)
5. Remove contamination: clusters 6 (NK), 7 (PMN), 10 (cPMN), 12 (germ B), 13 (act B), 14 (pro B)
6. Classify each cell: Paired TRA+TRB | Single chain | γδ | No TCR
7. Keep αβ TCR cells (Paired + Single chain) → save

### Key parameters
| Parameter | Value |
|---|---|
| PCA dims | 1:20 |
| Clustering resolution | 0.3 |
| αβ TCR filter | Paired + Single chain |

### Outputs
| File | Description |
|---|---|
| `objects/tcell_only_unclustered01.RDS` | T cell subset before clustering |
| `objects/tcell_only_clustered_unfilter02.rds` | All T cell clusters (pre-filter) |
| `objects/tcell_top500_markers_initial.csv` | Top 500 markers per initial cluster |
| `objects/cluster9_TCR_pairing_check.csv` | TCR pairing breakdown for cluster 9 |
| `objects/tcr_overall_filtered.csv` | Overall TCR status summary |
| `objects/tcell_TCRab_only03.rds` | **αβ TCR T cells → input to Script 2** |
| `figures/elbow_plot_tcell.svg` | PCA elbow plot |
| `figures/TCR_status_by_condition.svg` | TCR pairing proportions per condition |
| `figures/cell_counts_per_cluster_per_condition.svg` | Cell counts per cluster |

---

## Script 2 — `Tcell_ab_analysis_reclustering.R` (Stage 4+)

Set `start_from_clustered <- TRUE` to skip reclustering and load from checkpoint (fast re-run).

### Steps
1. **Stage 4a** — Recluster αβ T cells (HVG with TCR genes excluded, 8 PCs, res 0.8)
2. **USER CONFIGURATION** — review numbered UMAP + heatmap, set cluster labels
3. **Stage 4b** — Cluster refinement:
   - Cluster 7: MAIT GMM (Zbtb16 + Trav1) → remove innate-like cells, confirm Th17-like
   - Cluster 2: sub-clustering + Tbx21/Batf test → retained as Tfh/Tph-like (identity inconclusive)
4. **Stage 4c** — Final figures: UMAP, pies, heatmap, violins, stacked violin, identity heatmap

### Key parameters
| Parameter | Value |
|---|---|
| HVG | 3000 (2,887 after TCR-gene exclusion) |
| PCA dims | 1:8 (of 50 computed) |
| Clustering resolution | 0.8 |
| UMAP n.neighbors | 30 |
| UMAP min.dist | 0.15 |

### Cell type annotations (current)
| Cluster | Label |
|---|---|
| 0 | Th2-like |
| 1 | Teff CD8-like |
| 2 | Tfh/Tph-like |
| 3 | Tregs-like |
| 4 | Teff CD8-like |
| 5 | Naive T-like |
| 6 | Trm type mixed CD8+ CD4 (CD4/CD8-mixed; resolved in §4.8c) |
| 7 | Th17-like (after MAIT removal) |
| 8 | Trm/Exh-like CD8 |
| 9 | Teff/Tem-like CD8 |
| 10 | Prolif T |

### Outputs
| File | Description |
|---|---|
| `objects/tcell_TCRab_clustered04a.rds` | Clustered checkpoint |
| `objects/tcell_TCRab_annotated_final04.rds` | Final annotated object (9,251 cells; 27 MAIT/innate-like cells flagged `celltype_merged = NA`, not yet removed) |
| `objects/tcell_TCRab_annotated_final04_clean.rds` | Same object with the 27 NA cells physically removed (9,224 cells) — produced by `Tcell_remove_MAIT.R`; used by all downstream modules |
| `objects/cl7_innate_MAIT.rds` | MAIT cells identified in cluster 7 |
| `objects/cl2_Tph_subcluster.rds` | Cluster 2 sub-clustering result |
| `objects/cl2_subcluster_Tbx21_Batf_means.csv` | Per-sub-cluster Tbx21/Batf means |
| `objects/MAIT_contamination_per_condition.csv` | MAIT fraction per condition |
| `figures/UMAP_TCRab_annotated.svg` | Main annotated UMAP |
| `figures/UMAP_TCRab_byCondition.svg` | UMAP split by condition |
| `figures/Pie_TCRab_byCondition.svg` | Cell type composition per condition |
| `figures/heatmap_TCRab.svg` | Marker heatmap (avg exp, z-scored) |
| `figures/violin_markers.svg` | 27-gene violin across cell types |
| `figures/violin_markers_stacked.svg` | Stacked violin (ENC style) |
| `figures/heatmap_Tph_identity.svg` | Tfh/Tph-like identity gene heatmap |
| `figures/Supp_MixedTrm_annotation.svg` | Supplementary — justification of the cluster-6 `Trm type mixed CD8+ CD4` annotation (§4.8c) |
| `objects/cl6_mixed_lineage_composition.csv` | Cluster-6 CD4/CD8 composition, raw + kNN-imputed (§4.8c) |
| `objects/cl6_doublet_check.csv` | scDblFinder doublet rate per cell type (cluster 6 = 1.2%, lowest) |
| `figures/cl2_Tph_featureplot.svg` | Cluster 2 marker feature plots |
| `figures/cl2_subcluster_resolution_sweep.svg` | Sub-clustering resolution comparison |

---

## Diagnostic — "Trm-like CD8" (cluster 6) identity check (`Tcell_Trm_CD8_identity_check.R`)

Prompted by Module 05 flagging cluster 6 (`RNA_snn_res.0.8` = 6, provisional `celltype_merged` = "Trm-like CD8", 781 cells) as only ~50% CD8-marker-positive. This script is the **exploratory record** — it changes no labels itself and reads `objects/tcell_TCRab_annotated_final04_clean.rds`. The conclusion below was **acted on**: the canonical re-annotation and the supplementary figure now live in `Tcell_ab_analysis_reclustering.R` **section 4.8c** (see *Resolution* below).

### What the evidence shows

1. **The label is only partially supported.** In the cluster's own FindAllMarkers list, `Cd8a` ranks **146** (log2FC 0.96; expressed in only 41% of cells vs 33% elsewhere), `Cd8b1` ranks **330** (48% vs 46% — essentially flat), and `Cd4` is **not a marker at all**. The genuine top markers are lineage-neutral: immediate-early/activation (`Fos`, `Egr1`, `Tnf`, `Hspa1b`), tissue-resident/IEL (`Itga1` rank 46, `Id2` rank 80, `Itgae`, `Cd69`, `Cxcr6`, `Jaml`, `Adgrg1`, `Cdh1`), and — notably — **γδ TCR genes** (`Trgc2` rank 2, `Trgv2` rank 13). So "Trm/IEL" is well-supported; "CD8" rested on weakly-enriched `Cd8a`.

2. **It is a genuine ~50/50 CD4/CD8 mixture, not doublets.** Per-cell: 49% Cd8⁺/Cd4⁻, 26% Cd4⁺/Cd8⁻, only **1% double-positive** (rules out doublets), 24% double-negative (dropout). Dropout-robust CD4-vs-CD8 module scores give 53% CD4 / 45% CD8 / 1% ambiguous, and the module-score scatter (`figures/Trm_CD8_module_scores.svg`) shows **two cleanly separated arms** — distinct CD4 and CD8 cells, not a continuum. The CD4 fraction is **not Treg** (1% Foxp3⁺).

3. **Unconventional / innate-like features.** 47% of the cluster transcribes a γδ TCR gene (mostly `Trgv2` 41%; vs 21% object-wide) despite all cells having αβ TCRs, and the activation/immediate-early score is elevated (0.68 vs 0.49). Single-chain TCR fraction (33%) is **not** elevated vs the object (35%), so the ambiguity is not just poor TCR recovery.

4. **Sub-clustering does not separate the lineages.** Sub-clustering cluster 6 (silhouette-best res 0.2 → 2 sub-clusters) splits on **activation state** (sub-cluster 1 = `Atf3`, `Cdkn1a`, `Ccl4`, `Dusp1`, `Nr4a1`, `Fosb`, `Zfp36` — stress/IEG), not on CD4/CD8 — both sub-clusters stay lineage-mixed. Lineage is not the dominant transcriptional axis here, which is exactly why a CD4-Trm and a CD8-Trm population merged into one cluster at res 0.8 / 8 PCs.

### Resolution (applied)

**✅ RESOLVED — cluster 6 re-annotated `Trm type mixed CD8+ CD4`.** Acting on the diagnostic above, the cluster is now labeled honestly as a lineage-mixed Trm/IEL-like population inside the reclustering pipeline (`Tcell_ab_analysis_reclustering.R`, **section 4.8c**, parallel to the MAIT-GMM and Tph refinement blocks), and the label is propagated to every downstream object and figure.

- **Doublets ruled out** (`Tcell_doublet_check.R`): scDblFinder (`samples = orig.ident`) puts cluster 6 at **1.2% doublets (9 of 781) — the lowest of any T-cell cluster** (overall 4.4%); see `objects/cl6_doublet_check.csv`. This is a genuine biological mixture, not a technical artifact.
- **Composition (CD8-majority mixed):** ~**60% CD8 / ~40% CD4**, concordant across methods — raw marker-positive resolvable cells 66/34; kNN label transfer (majority vote; local PCA, k=30, over the 187 marker-negative dropout cells; 8 Cd4⁺Cd8⁺ cells excluded) 61/39. Reported as robust ~60/40 rather than a single fragile integer split.
- **Two functionally distinct arms, not a "CD8-helper" population:** the CD8 arm is cytotoxic (`Nkg7`/`Ifng`/`Gzmb`/`Prf1`-high, helper-low) and the CD4 arm is a helper (`Cd40lg`/`Icos`/`Maf`-high) — verified as part of a "is cluster 6 CD8 T-helper?" check.
- **Why "mixed" rather than a per-cell split:** ~24% of cells carry no lineage marker (dropout) and cannot be split reliably; sub-clustering separates on activation state, not lineage. The mixed label makes cluster 6 **neither CD4 nor CD8**, so it correctly drops out of all lineage-specific (CD4/CD8) downstream analyses.
- **Supplementary figure:** `figures/Supp_MixedTrm_annotation.svg` (panels A composition, B two-arm module-score scatter, C lineage/Trm-IEL/innate program dotplot, D cytotoxic-CD8-vs-helper-CD4 contrast) justifies the annotation; evidence table `objects/cl6_mixed_lineage_composition.csv`.

Modules 04–06 and the additional analyses (07–09) were re-run after the relabel; the CD4-only beta-gene/clonotype core (Module 05 TABLE1–11) is unaffected because the mixed cluster was never part of the CD4-expanded population.

**Companion feature plots (`Tcell_Trm_CD8_feature_umap.R`):** plots the CD8/CD4 module scores and the raw per-cell lineage category on (A) the main T-cell UMAP and (B) cluster 6's own local re-embedding, to see whether the double-negative "dropout" cells colocalize with the CD8 or CD4 arm. Diagnostic only; reads the `_clean` object, changes nothing. Outputs: `figures/Trm_CD8_featureplot_mainUMAP.svg`, `figures/Trm_CD8_featureplot_localUMAP.svg`.

### Outputs
| File | Description |
|---|---|
| `objects/Trm_CD8_marker_means_by_cluster.csv` | Marker mean + % expressing, every cell type |
| `objects/Trm_CD8_top100_marker_audit.csv` | Cluster-6 FindAllMarkers list, ranked + lineage/Trm/γδ/IEG tags |
| `objects/Trm_CD8_lineage_composition.csv` | Crude >0 CD4/CD8 cross-tab (49/26/1/24%) |
| `objects/Trm_CD8_module_scores.csv` | Per-cell CD8 & CD4 module scores + call (53% CD4 / 45% CD8) |
| `objects/Trm_CD8_gd_activation_diagnostics.csv` | γδ-gene %, activation score, single-chain % vs whole object |
| `objects/Trm_CD8_subcluster_markers.csv` | Cd4/Cd8/Trm profile per sub-cluster |
| `objects/Trm_CD8_subcluster_top_DE.csv` | Genes distinguishing sub-clusters (activation, not lineage) |
| `figures/Trm_CD8_identity_heatmap.svg` | Marker profile by cell type (z-scored) — cluster 6 intermediate |
| `figures/Trm_CD8_identity_dotplot.svg` | Marker dot plot by cell type |
| `figures/Trm_CD8_lineage_composition.svg` | Cluster 6 CD4/CD8 composition bar |
| `figures/Trm_CD8_module_scores.svg` | Per-cell CD8 vs CD4 module-score scatter (two arms) |
| `figures/Trm_CD8_subcluster_umap.svg` | Sub-cluster UMAP + Cd4/Cd8a/Cd8b1 feature UMAPs |

---

## Supplementary figure — Tfh/Tph-like marker violins (`Tph_marker_violin.R`)

Grouped violin panel supporting the Cluster 2 **Tfh/Tph-like** annotation. Reads `objects/tcell_TCRab_annotated_final04_clean.rds` and plots four curated marker groups — B-cell chemoattraction / core Tph (`Cxcl13`, `Icos`, `Il21`, `Cd40lg`, `Pdcd1`), exhaustion-associated TFs (`Tox`, `Batf`, `Prdm1`, `Bcl6`), Tfh-exclusion markers (`Cxcr5`, `Cxcr6`, `Sell`, `Ccr7`), and costimulation/survival (`Tnfrsf4`, `Tnfsf8`, `Cd200`) — across all annotated cell types. Figure only, no object modified.

| File | Description |
|---|---|
| `figures/Tph_marker_violin_combined.svg` | Tfh/Tph-like marker violin panel (4 gene groups) |

---

## Figure-assembly scripts (SVG inlining)

These scripts do **no analysis** — they combine already-generated panel SVGs into publication figures by true SVG inlining (each element stays a real, editable node; per-panel ID prefixing avoids conflicts). Run only after the upstream figures exist.

| Script | Combines | Output |
|---|---|---|
| `assemble_figure_panels.R` | Global UMAP+pies (Module 01) + TCR status + annotated UMAP + by-condition UMAP + pies | `figures/Figure_Tcell_combined.svg` |
| `assemble_panel_D.R` | Annotated UMAP + by-condition UMAPs | `figures/panel_D.svg` |
| `assemble_extended_fig.R` | Marker heatmap + stacked violin | `figures/Figure_extended_heatmap_violin.svg` |
