# Module 05 — TCR repertoire analysis

TCRβ V–J gene usage, β-gene prevalence and repertoire overlap between recovery (REC) and recurrence (RCU), TCRs shared between REC and RCU (exact and convergent), and a TCR-recovery QC. Produces Figure 6B–F and Figure S2D.


## Inputs
| File | From |
|---|---|
| `../04_Tcell_Expansion_Analysis/objects/tcell_TCRab_expansion_annotated05.rds` | Module 04 (9,224 αβ T cells with `expansion_status`) |
| `../03_Tcell_Analysis/objects/tcell_only_clustered_unfilter02.rds` | Module 03 (all T cells incl. no-TCR; dropout analysis only) |

**CD4 definition.** `CD_broad` is derived on the fly from `celltype_merged`: CD4 = Tfh/Tph-like, Th17-like, Th2-like, Tregs-like; CD8 = Teff CD8-like, Trm/Exh-like CD8, Teff/Tem-like CD8; "Trm type mixed CD8+ CD4" = mixed (excluded from CD4/CD8-specific analyses). Naive T-like and Prolif T are not assigned a lineage.

## Scripts (run order — `../run_all.sh 05`)
| # | Script | Purpose | Main outputs |
|---|---|---|---|
| 1 | `TCR_clonotype_tracking.R` | Per-cell αβ chain table of expanded CD4 cells in REC + RCU; clonotype overlap and β-gene tables | `Clonotype_TABLE1`–`TABLE11` |
| 2 | `Beta_Gene_Prevalence_Analysis.R` | Standardized β-gene dataset and analysis set (reads TABLE1) | `BETA_GENE_ACELL.csv`, `TABLE_A1`, `TABLE1_celltype_by_condition*` |
| 3 | `Beta_Gene_Prevalence_Within_Condition.R` | Prevalence per condition, REC − RCU difference, enrichment class | `TABLE_A3`, `TABLE3`, `figures/beta_gene_*` |
| 4 | `Beta_Gene_Jaccard_Global.R` | REC/RCU Jaccard of β-gene sets + subsampling sensitivity | `TABLE2_jaccard_global.csv`, `TABLE_A2*`, `figures/beta_gene_jaccard_sensitivity.svg` |
| 5 | `Beta_Gene_Jaccard_By_CellType.R` | Jaccard within each CD4 cell type; shared genes by cell type | `TABLE5_jaccard_by_celltype.csv`, `TABLE_shared_genes_by_celltype.csv` |
| 6 | `TCR_VDJ_chord.R` | TCRβ V–J chord diagrams | `figures/CD4_combined_chord.svg`, `ALLcells_combined.svg`, `Expanded_celltype_combined.svg` |
| 7 | `TCR_Shared_Clonotype_CTstrict.R` | Exact-nucleotide clonotypes shared REC↔RCU (CD4 expanded, paired) | `Clonotype_TABLE12` |
| 8 | `TCR_Convergent_Clonotype_Summary.R` | Convergent (same amino acid, ≥ 2 nucleotide variants) clonotypes present in REC and RCU | `Clonotype_TABLE13` |
| 9 | `TCR_Convergent_Variants_Distribution.R` | All cells carrying the dominant clonotype's 3 shared nucleotide variants (V1–V3) | `Clonotype_TABLE14`, `TABLE14b`, `TABLE17b` |
| 10 | `TCR_Convergent_Variant_Heatmaps.R` | Per-cell marker heatmaps of V1–V3 and V1 sequence table | `figures/convergent_V1_heatmap.svg`, `convergent_variant_heatmaps.svg`, `Clonotype_TABLE16`, `TABLE17` |
| 11 | `TCR_Dropout_Analysis.R` | TCR recovery vs sequencing depth (logistic GLM + mixed model) | `figures/tcr_dropout_*.svg`, `objects/tcr_dropout_*.csv` |

Scripts 2–5 read TABLE1 from script 1; scripts 10 reads TABLE14 from script 9. All tables are written to `objects/`.

## Figure map
| Paper panel | Output |
|---|---|
| Fig. 6B | `figures/CD4_combined_chord.svg` (top 12 TRBV × top 12 TRBJ per panel, ribbons normalized within condition) |
| Fig. 6C | `figures/beta_gene_top_enriched_shared_bar.svg` |
| Fig. 6D | `objects/TABLE5_jaccard_by_celltype.csv` + global row from `TABLE2_jaccard_global.csv` |
| Fig. 6E | `objects/Clonotype_TABLE17_V1_TCR_sequence.csv` (V1 TCR sequences) |
| Fig. 6F | `figures/convergent_V1_heatmap.svg` (V1 marker heatmap; values in `objects/Clonotype_TABLE16_convergent_variant_heatmap_matrix.csv`) |
| Fig. S2D | `figures/tcr_dropout_depth_by_status.svg`, `tcr_dropout_percell_logit.svg`, `tcr_dropout_rate_vs_depth.svg` |

`convergent_variant_heatmaps.svg` (V1–V3 side by side), `ALLcells_combined.svg`, `Expanded_celltype_combined.svg` and the scatter plots are produced alongside but not used in the paper.

## Key results

**β-gene analysis set** (expanded CD4 T cells, REC + RCU): 678 cells, 33 (4.9%) without a β chain, **645 analysed** (REC 442, RCU 203); 125 distinct TCRβ gene combinations.

**Prevalence classes** (REC − RCU difference): REC-enriched 39, RCU-enriched 40, shared 6, rare 40. Enriched = difference > 0.5 percentage points; shared = |difference| ≤ 0.5 and ≥ 3 cells in either condition.

**Jaccard similarity** of REC vs RCU β-gene sets: global **0.200** (90 REC, 60 RCU, 25 shared, 125 union). Subsampling REC to n = 203 (1,000×): median 0.194, 95% interval 0.162–0.220, so the overlap is not an artifact of REC's larger size.

| CD4 cell type | REC | RCU | Shared | Union | Jaccard |
|---|---|---|---|---|---|
| Tfh/Tph-like | 63 | 23 | 10 | 76 | 0.132 |
| Tregs-like | 34 | 23 | 5 | 52 | 0.096 |
| Th17-like | 17 | 21 | 2 | 36 | 0.056 |
| Th2-like | 11 | 10 | 1 | 20 | 0.050 |

**Shared and convergent clonotypes.** "Shared" depends on how strictly "same TCR" is defined and on which cells are counted, so the stored tables are not directly comparable:

| Table | Identity | Population | Shared REC↔RCU |
|---|---|---|---|
| TABLE2 | α + β gene tokens | CD4 expanded, incl. single-chain cells | 9 clonotypes (69 cells); 138 REC-only, 77 RCU-only |
| TABLE12 | exact `CTstrict` (genes + CDR3 nucleotide) | CD4 expanded, paired | 3 clonotypes (16 cells) |
| TABLE13 | genes + CDR3 amino acid, ≥ 2 nucleotide variants | all T cells | 12 clonotypes |

The 3 exact shared clonotypes are nucleotide variants of one convergent TCR: α `TRAV6N-7.TRAJ13` / CDR3 `CALGANSGTYQRF`, β `TRBV5.TRBD1.TRBJ1-4` / CDR3 `CASSQEGQGTNERLFF` (29 REC + 25 RCU cells at the amino-acid level, from 8 and 9 nucleotide variants; largest single clone 6 cells in each condition — different nucleotide sequences). Its three shared variants (V1–V3) total 47 cells across all four conditions (TABLE14). **V1**, the most expanded, is carried by 12 cells (7 REC, 4 RCU, 1 PD; 8 Tfh/Tph-like, 3 mixed-Trm cluster, 1 Th2-like). Because each condition is a different cohort of mice, shared TCRs are public/convergent TCRs, not tracked clones.

**TCR dropout.** Of 14,530 contamination-removed T cells, 63.7% have a recovered TCR (BAS 39.6%, PD 73.7%, REC 62.2%, RCU 56.8%). Recovery rises with GEX UMI depth independent of condition: OR 49 per 10-fold UMI (95% CI 41–59; GLM), 46 (38–56) in a mixed model with sample as random intercept. Condition remains associated after adjusting for depth (likelihood-ratio χ²₃ = 211, p = 1.5 × 10⁻⁴⁵). Non-recovery is therefore largely technical, and BAS (lowest depth and recovery) needs the most caution.

## Cluster label vs single-cell identity (note for reviewers)

`celltype_merged` is a cluster-level label, so individual cells can disagree with their cluster's majority identity. Lineage claims are made at the cluster level and quantified; where a claim rests on specific rare cells (e.g. the convergent-clonotype cells in the mixed-Trm cluster), those cells were checked directly with a lineage-marker check. The mixed-Trm cluster itself is documented in Module 03 (§4.8c of `Tcell_ab_analysis_reclustering.R`; doublets ruled out by `Tcell_doublet_check.R`, scDblFinder 1.2%, lowest of any T-cell cluster).
