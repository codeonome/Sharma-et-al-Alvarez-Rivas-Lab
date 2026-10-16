# Sharma et al. — Alvarez-Rivas Lab

Analysis code for:

> Sharma A, Olmedo V, Rojas L, *et al.*, Alvarez Rivas C. **Clonally Convergent Tph/Tfh-like CD4⁺ T Cells are Shared Across Recovery and Recurrence in Experimental Periodontitis.** *[journal, year — to be added]*

The repository has two parts:

| Folder | Contents |
|---|---|
| [`Mouse-data-analysis/`](Mouse-data-analysis/) | **Sensitized recurrent periodontitis (SeRPD) mouse model.** Paired 10x 5′ scRNA-seq and TCR-seq of gingival CD45⁺ cells at baseline (BAS), periodontitis (PD), recovery (REC) and recurrence (RCU). Covers global, myeloid and T-cell clustering; clonal expansion; the TCR repertoire and convergent clonotypes; and the cross-species Tph/Tfh comparison. |
| [`Human-data-analysis/`](Human-data-analysis/) | **Reanalysis of human gingiva scRNA-seq** (GSE164241; 21 biopsies, 18 donors, healthy vs periodontitis). Covers T-cell subclustering (C1–C9), differential abundance and pseudobulk differential expression. The pipeline is step-wise and verified to reproduce the published results. |

Each folder's README explains how to run it, describes its computational environment, and maps figure panels to scripts and output files.

## Figures and tables → code
Each row links to the script that produces the panel and to the table of values behind it.
- Human tables are the published results in `Human-data-analysis/published_results/`; the pipeline regenerates them exactly.
- Mouse tables are written by each module to its `objects/` folder.

Panels drawn directly from a Seurat object have no separate table. The objects are deposited at Zenodo (see below).

| Item | Content | Code | Data table(s) |
|---|---|---|---|
| Fig. 1A–B | Imaging mass cytometry, human gingiva | not in this repository (QuPath/Fiji) | — |
| Fig. 1C | Human T-cell UMAP; density in GM and PD | [08_figures.R](Human-data-analysis/scripts/08_figures.R) | [cluster_renaming_map.csv](Human-data-analysis/published_results/04_Tcell_Final/results/cluster_renaming_map.csv) |
| Fig. 1D | Grouped marker dot plot (CD8 · MAIT · CD4) | [08_figures.R](Human-data-analysis/scripts/08_figures.R) | computed from the final T-cell object |
| Fig. 1E | Cluster composition in GM and PD | [08_figures.R](Human-data-analysis/scripts/08_figures.R) | [donut_pooled_percentages.csv](Human-data-analysis/published_results/04_Tcell_Final/results/donut_pooled_percentages.csv), [sample_x_cluster_counts.csv](Human-data-analysis/published_results/04_Tcell_Final/results/sample_x_cluster_counts.csv) |
| Fig. 1F | Pseudobulk DE, PD vs GM (all T cells) | [09_pseudobulk_de.R](Human-data-analysis/scripts/09_pseudobulk_de.R) | [DE_PD_vs_GM_allTcells.csv](Human-data-analysis/published_results/04_Tcell_Final/results/DE_PD_vs_GM_allTcells.csv) |
| Fig. 1G | Differential abundance (odds ratio, FDR) | [07_differential_abundance.R](Human-data-analysis/scripts/07_differential_abundance.R), [08_figures.R](Human-data-analysis/scripts/08_figures.R) | [DA_bar_blocked_values.csv](Human-data-analysis/published_results/04_Tcell_Final/results/DA_bar_blocked_values.csv), [DA_propeller_donor_blocked.csv](Human-data-analysis/published_results/04_Tcell_Final/results/DA_propeller_donor_blocked.csv) |
| Fig. 2, Fig. 3 | Flow cytometry, micro-CT, qPCR, TRAP, adoptive transfer | not in this repository (FlowJo/GraphPad Prism) | — |
| Fig. 4A | Experimental design (schematic) | — | — |
| Fig. 4B | Mouse immune landscape by phase | [Global_analyis_postprocessing.R](Mouse-data-analysis/01_Global_Analysis/scripts/Global_analyis_postprocessing.R) | computed from the global object |
| Fig. 4C–E | T-cell subclusters; by phase; composition | [Tcell_ab_analysis_reclustering.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_ab_analysis_reclustering.R) | [markers_TCRab_res0.8_top500.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/markers_TCRab_res0.8_top500.csv) |
| Fig. 5A–E | Clonal expansion | [Tcell_expansion_analysis.R](Mouse-data-analysis/04_Tcell_Expansion_Analysis/scripts/Tcell_expansion_analysis.R) | [expanded_pool_composition.csv](Mouse-data-analysis/04_Tcell_Expansion_Analysis/objects/expanded_pool_composition.csv), [expanded_pool_denominators.csv](Mouse-data-analysis/04_Tcell_Expansion_Analysis/objects/expanded_pool_denominators.csv), [clonotype_by_sample_summary.csv](Mouse-data-analysis/04_Tcell_Expansion_Analysis/objects/clonotype_by_sample_summary.csv), `expansion_{BAS,PD,REC,RCU}.csv` |
| Fig. 6A | TCR structure (schematic) | — | — |
| Fig. 6B | TCRβ V–J chord diagrams, CD4 T cells | [TCR_VDJ_chord.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/TCR_VDJ_chord.R) | [Clonotype_TABLE1_per_cell_CD4_expanded_REC_RCU.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/Clonotype_TABLE1_per_cell_CD4_expanded_REC_RCU.csv) |
| Fig. 6C | β-chain prevalence, REC vs RCU | [Beta_Gene_Prevalence_Within_Condition.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/Beta_Gene_Prevalence_Within_Condition.R) | [TABLE3_top_enriched_shared_genes.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/TABLE3_top_enriched_shared_genes.csv), [TABLE_A3_gene_prevalence_within_condition.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/TABLE_A3_gene_prevalence_within_condition.csv) |
| Fig. 6D | Jaccard index per CD4 cluster | [Beta_Gene_Jaccard_By_CellType.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/Beta_Gene_Jaccard_By_CellType.R), [Beta_Gene_Jaccard_Global.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/Beta_Gene_Jaccard_Global.R) | [TABLE5_jaccard_by_celltype.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/TABLE5_jaccard_by_celltype.csv), [TABLE2_jaccard_global.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/TABLE2_jaccard_global.csv) |
| Fig. 6E | Convergent clonotype V1: TCR sequences | [TCR_Convergent_Variant_Heatmaps.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/TCR_Convergent_Variant_Heatmaps.R) | [Clonotype_TABLE17_V1_TCR_sequence.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/Clonotype_TABLE17_V1_TCR_sequence.csv) |
| Fig. 6F | Convergent clonotype V1: marker heatmap | [TCR_Convergent_Variant_Heatmaps.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/TCR_Convergent_Variant_Heatmaps.R) | [Clonotype_TABLE16_convergent_variant_heatmap_matrix.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/Clonotype_TABLE16_convergent_variant_heatmap_matrix.csv) |
| Fig. 7A | B-helper module genes, human and mouse clusters | [04_bhelper_featureumap.R](Mouse-data-analysis/06_CrossSpecies_Tph/scripts/04_bhelper_featureumap.R) | [bhelper_gene_by_cluster_expression.csv](Mouse-data-analysis/06_CrossSpecies_Tph/objects/bhelper_gene_by_cluster_expression.csv) |
| Fig. 7B–C | B-helper score on UMAP; score per cluster | [04_bhelper_featureumap.R](Mouse-data-analysis/06_CrossSpecies_Tph/scripts/04_bhelper_featureumap.R) | [bhelper_score_by_cluster.csv](Mouse-data-analysis/06_CrossSpecies_Tph/objects/bhelper_score_by_cluster.csv), [bhelper_target_vs_cluster_stats.csv](Mouse-data-analysis/06_CrossSpecies_Tph/objects/bhelper_target_vs_cluster_stats.csv) |
| Supp. Fig. 1 | Flow cytometry gating and frequencies | not in this repository (FlowJo/GraphPad Prism) | — |
| Supp. Fig. 2A–B | Global clusters; marker dot plot | [Global_analyis_postprocessing.R](Mouse-data-analysis/01_Global_Analysis/scripts/Global_analyis_postprocessing.R) | [cca0.05_top100_markers.csv](Mouse-data-analysis/01_Global_Analysis/objects/cca0.05_top100_markers.csv) |
| Supp. Fig. 2C | TCR chain pairing per phase | [Tcell_analysis_reclustering.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_analysis_reclustering.R) | [tcr_overall_filtered.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/tcr_overall_filtered.csv) |
| Supp. Fig. 2D | TCR recovery vs sequencing depth | [TCR_Dropout_Analysis.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/TCR_Dropout_Analysis.R) | [tcr_dropout_detection_by_condition.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/tcr_dropout_detection_by_condition.csv), [tcr_dropout_rate_by_depth_decile.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/tcr_dropout_rate_by_depth_decile.csv), [tcr_dropout_model_tests.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/tcr_dropout_model_tests.csv) |
| Supp. Fig. 3A–C | Myeloid subclusters | [Myeloid_cells_analysis_refactored.R](Mouse-data-analysis/02_Myeloid_Analysis/scripts/Myeloid_cells_analysis_refactored.R) | [cca0.3_top500_markers_myeloid01.csv](Mouse-data-analysis/02_Myeloid_Analysis/objects/cca0.3_top500_markers_myeloid01.csv) |
| Supp. Fig. 4A–B | T-cell subtype heatmap and violins | [Tcell_ab_analysis_reclustering.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_ab_analysis_reclustering.R), [assemble_extended_fig.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/assemble_extended_fig.R) | [markers_TCRab_res0.8_top500.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/markers_TCRab_res0.8_top500.csv) |
| Supp. Fig. 4C | Tfh/Tph-like cluster markers | [Tcell_ab_analysis_reclustering.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_ab_analysis_reclustering.R) | [cl2_Tph_marker_summary.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/cl2_Tph_marker_summary.csv) |
| Supp. Fig. 4D | Trm-like CD4/CD8 mixed-lineage cluster | [Tcell_ab_analysis_reclustering.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_ab_analysis_reclustering.R), [Tcell_Trm_CD8_identity_check.R](Mouse-data-analysis/03_Tcell_Analysis/scripts/Tcell_Trm_CD8_identity_check.R) | [cl6_mixed_lineage_composition.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/cl6_mixed_lineage_composition.csv), [Trm_CD8_lineage_composition.csv](Mouse-data-analysis/03_Tcell_Analysis/objects/Trm_CD8_lineage_composition.csv) |
| Supp. Table 1 | Human differential abundance, all clusters | [07_differential_abundance.R](Human-data-analysis/scripts/07_differential_abundance.R) | [Propeller_DA_results.xlsx](Human-data-analysis/published_results/04_Tcell_Final/results/Propeller_DA_results.xlsx) (from [DA_propeller_donor_blocked.csv](Human-data-analysis/published_results/04_Tcell_Final/results/DA_propeller_donor_blocked.csv), [DA_all_tests_wide.csv](Human-data-analysis/published_results/04_Tcell_Final/results/DA_all_tests_wide.csv)) |
| Supp. Table 2 | Jaccard index per CD4 cluster | [Beta_Gene_Jaccard_By_CellType.R](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/Beta_Gene_Jaccard_By_CellType.R) | [TABLE5_jaccard_by_celltype.csv](Mouse-data-analysis/05_TCR_Repertoire_Analysis/objects/TABLE5_jaccard_by_celltype.csv) |
| Supp. Tables 3–4 | Convergent clonotype variants | [05_TCR_Repertoire_Analysis/scripts/](Mouse-data-analysis/05_TCR_Repertoire_Analysis/scripts/) | *[to be confirmed]* |

## Data and code availability
- **Mouse scRNA-seq and TCR-seq:** GEO [GSE350875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350875)
- **Human scRNA-seq:** GEO [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241) (Williams et al., *Cell* 2021)
- **Processed Seurat objects:** too large for git; deposited at Zenodo *[DOI — to be added]*. Each README lists which objects go where.
- **Software versions:** `Mouse-data-analysis/renv.lock` and `Human-data-analysis/environment/`

## Contact
Apoorva Sharma , (apoorva.nyu@gmail.com)
Dr. Carla Alvarez Rivas, Oral Biology Department, University at Buffalo (carlaalv@buffalo.edu)
