# Human gingival T cells in periodontitis — GSE164241 reanalysis

Code for the human single-cell analyses in **Sharma et al., *Clonally Convergent Tph/Tfh-like CD4⁺ T Cells are Shared Across Recovery and Recurrence in Experimental Periodontitis*** (Figure 1C–G, Supplementary Table 1, and the human panels of Figure 7).

The data are public scRNA-seq of human gingiva ([GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241); Williams et al., 2021): 21 biopsies from 18 donors, 13 healthy (GM) and 8 periodontitis (PD). Two donors gave repeat biopsies: PD164 (3) and PD161 (2). The pipeline goes from the integrated object to nine T-cell clusters (C1–C9). It then tests their differential abundance between PD and GM and runs pseudobulk differential expression.

## Figures and tables → code
The **published values** behind every panel are in [`published_results/04_Tcell_Final/results/`](published_results/04_Tcell_Final/results/); a run writes the same files under `output/`.

<details>
<summary><b>Figure 1, Supp. Table 1 and Fig. 7 → script and output files</b></summary>

| Manuscript item | Content | Script | Output (`output/…`) |
|---|---|---|---|
| **Fig. 1C** | T-cell UMAP, all subjects; density in GM and PD | [`08_figures.R`](scripts/08_figures.R) | `figures/08_figures/umap_tcell_final.svg`, `umap_density_by_condition.svg` |
| **Fig. 1D** | Grouped marker dot plot (CD8 · MAIT · CD4) | [`08_figures.R`](scripts/08_figures.R) | `figures/08_figures/dotplot_grouped.svg` |
| **Fig. 1E** | Cluster composition, GM and PD (pooled cells) | [`08_figures.R`](scripts/08_figures.R) | `figures/08_figures/donut_GM.svg`, `donut_PD.svg`; `results/08_figures/donut_pooled_percentages.csv` |
| **Fig. 1F** | Pseudobulk DE, PD vs GM, all T cells; TOX2, BCL6, CXCR5 and CXCL13 highlighted | [`09_pseudobulk_de.R`](scripts/09_pseudobulk_de.R) | `figures/09_pseudobulk_de/DE_volcano_PD_vs_GM_allTcells.svg`; `results/09_pseudobulk_de/DE_PD_vs_GM_allTcells.csv` |
| **Fig. 1G** | Differential abundance: odds ratio (95% CI) and FDR, donor-blocked propeller | [`07_differential_abundance.R`](scripts/07_differential_abundance.R), [`08_figures.R`](scripts/08_figures.R) | `figures/08_figures/DA_bar_blocked.svg`; `results/07_differential_abundance/DA_propeller_donor_blocked.csv` |
| **Supp. Table 1** | DA for all clusters: donor-blocked (primary), sample-level and donor-level models | [`07_differential_abundance.R`](scripts/07_differential_abundance.R) | `results/07_differential_abundance/DA_propeller_donor_blocked.csv`, `DA_all_tests_wide.csv` (assembled as `published_results/04_Tcell_Final/results/Propeller_DA_results.xlsx`) |
| **Fig. 7A–C** (human) | B-helper module in human T-cell clusters | [`06_CrossSpecies_Tph`](../Mouse-data-analysis/06_CrossSpecies_Tph/scripts/04_bhelper_featureumap.R) | inputs: `objects/06_tcell_final.rds`, `results/04_tcell_reclustering/tcell_subset_markers.csv` |

</details>

Figure 1A–B (imaging mass cytometry) is not produced by code in this repository.

### Cluster names
The manuscript uses functional names. The code and the tables use marker-based labels (`cluster_id`). Both refer to the same clusters, identified by the code C1–C9.

<details>
<summary><b>C1–C9: manuscript names vs code labels</b></summary>

| Code | Manuscript (Fig. 1) | Label in code and tables | Lineage |
|---|---|---|---|
| C1 | effector/cytotoxic-like | CD8 TEM (GZMK/GZMH) | CD8 |
| C2 | naive/central memory-like | CD4 naive/CM | CD4 |
| C3 | resting/effector-memory-like | CD8 TEM (ZEB2-hi) | CD8 |
| C4 | Th17/Treg-like | Th17 | CD4 |
| C5 | Treg-like | Treg | CD4 |
| C6 | activated/memory-like | CD4 resting/CM (TSHZ2+) | CD4 |
| C7 | terminally differentiated effector-like | CD8 cytotoxic (GZMB/GNLY) | CD8 |
| C8 | MAIT | MAIT | MAIT |
| C9 | Tfh/Tph-like | CD4 Tph-like (CXCL13+) | CD4 |

</details>

## Steps
| Step | Script | What it does | Feeds | Env. |
<details>
<summary><b>Steps 00–09 (Env. 1 = R 4.2.0 / Seurat 5.1.0 / harmony 0.1.0; Env. 2 = R 4.5.1 / limma / edgeR / speckle)</b></summary>

|---|---|---|---|---|
| 00 | [`00_integration.R`](scripts/00_integration.R) | Per-sample QC (200–6,000 genes, < 20% mitochondrial), merge, CCA and Harmony integration, clustering | — | 1 |
| 01 | [`01_qc_integrated_object.R`](scripts/01_qc_integrated_object.R) | Sanity checks on the integrated object | — | 1 |
| 02 | [`02_global_clustering.R`](scripts/02_global_clustering.R) | Resolution sweep on the Harmony embedding; res 0.1 = 12 clusters | — | 1 |
| 03 | [`03_immune_subset.R`](scripts/03_immune_subset.R) | Immune clusters re-integrated and reclustered (res 0.2, 13 clusters); T cells = immune clusters 0 and 1 | — | 1 |
| 04 | [`04_tcell_reclustering.R`](scripts/04_tcell_reclustering.R) | T cells re-integrated and reclustered; res 0.4 = 11 clusters; markers | Fig. 7 (markers) | 1 |
| 05 | [`05_tcell_validation.R`](scripts/05_tcell_validation.R) | Validation of clusters 3 (Treg/Th17), 8 (Tph/Tfh) and 9 (QC) | — | 1 |
| 06 | [`06_tcell_cleanup_annotation.R`](scripts/06_tcell_cleanup_annotation.R) | Removes Ig-high doublets and cluster 9, merges 10 into 0, labels C1–C9 | Fig. 1, 7 | 2 |
| 07 | [`07_differential_abundance.R`](scripts/07_differential_abundance.R) | Composition; Wilcoxon and propeller tests; **primary: donor-blocked propeller**; leave-one-donor-out; ρ sensitivity | Fig. 1G, Supp. Table 1 | 2 |
| 08 | [`08_figures.R`](scripts/08_figures.R) | Figure 1 panels C, D, E and G, plus supporting plots | Fig. 1C–E, G | 2 |
| 09 | [`09_pseudobulk_de.R`](scripts/09_pseudobulk_de.R) | Pseudobulk DE, PD vs GM: all T cells and each cluster C1–C9 | Fig. 1F | 2 |

</details>

Cell numbers: 110,166 integrated cells → 48,346 immune cells → 21,899 T cells → 20,802 T cells in C1–C9.

Step 00 was written and run by Tyler Chinsky. The deposited integrated object is the verified starting point for steps 01–09.

## Methods summary
- **Integration:**
  - Normalisation, 2,000 variable genes and PCA.
  - Harmony with sample as the batch; each biopsy is its own batch.
  - Clustering on the Harmony embedding (30 dims), global res 0.1.
- **Immune and T-cell subclustering:**
  - The immune clusters, and then the T cells, are each re-normalised, re-integrated with Harmony and reclustered: immune with 20 dims at res 0.2, T cells with 15 dims at res 0.4.
- **Cleanup:**
  - Cells with more than 10% of counts from immunoglobulin genes are removed as T/plasma doublets.
  - Cluster 9 is removed: it is dominated by ribosomal genes and has no lineage markers.
  - Cluster 10 is merged into cluster 0.
- **Differential abundance:**
  - Model: `speckle::propeller` on logit-transformed per-sample proportions, fitted at sample level with donor as a blocking factor (`limma::duplicateCorrelation`, consensus ρ = 0.726), with robust empirical Bayes.
  - Sample-level and donor-level models and Wilcoxon tests are reported for concordance.
- **Differential expression:**
  - Counts are summed per sample (samples with fewer than 20 cells are dropped), then `filterByExpr` and TMM.
  - Model: `edgeR::voomLmFit` with sample quality weights and donor blocking, then `eBayes(robust = TRUE)`.
  - 10,796 genes are tested in the all-T-cell analysis.

## Folder layout
```
Human-data-analysis/
├── config.R            all paths and parameters (the only file to edit)
├── R/helpers.R         shared functions (plotting, markers, DA models)
├── scripts/00-09_*.R   one script per step
├── run/                run scripts for steps 00-05 (Environment 1) and 06-09 (Environment 2)
├── environment/        R and package versions for both environments, renv.lock
├── verify/             comparison of outputs with the published results, and its reports
├── published_results/  tables and figures as published (the verification reference)
└── output/             created by a run (git-ignored): objects/, results/<step>/, figures/<step>/, logs/
```

## How to run
Run every script from this folder; each one reads `config.R`. Set `PIPELINE_OUT` to change the output location. The analysis uses **two R environments**, recorded in [`environment/`](environment/ENVIRONMENT.md), each with an `renv.lock` to recreate it.

**1. Steps 00–05: Environment 1** (R 4.2.0, Seurat 5.1.0, harmony 0.1.0; ~100 GB RAM):
```bash
mkdir -p output/objects
cp <deposited>/integrated_seurat_object.rds output/objects/00_integrated.rds
bash run/run_steps00-05.sh                                   # steps 01-05
# or from raw Cell Ranger output:
RAW_10X_DIR=/path/to/cellranger STEPS="00 01 02 03 04 05" bash run/run_steps00-05.sh
```
harmony 0.1.0 is what reproduces the published clusters exactly; harmony ≥ 1.0 gives a different embedding. On an HPC cluster, submit the same script with your scheduler (e.g. `sbatch`).

**2. Steps 06–09: Environment 2** (R 4.5.1, limma 3.66.0, edgeR 4.8.2, speckle 1.10.0):
```bash
bash run/run_steps06-09.sh        # reads output/objects/04_tcell_reclustered.rds from steps 00-05
```

## Reproducibility check
`verify/verify_outputs.R` compares a run's outputs with [`published_results/`](published_results/). Tables are compared byte for byte. For objects it compares cells, clusters at every resolution (ARI) and embeddings. Object checks also need the deposited `.rds` files.
```bash
Rscript verify/verify_outputs.R published_results 06 07 08 09
```
On 2026-10-15 the pipeline was run from the deposited integrated object. **All 91 checks were identical:**
- Steps 02–05 (Environment 1): 47/47 ([report](verify/report_steps02-05.csv)).
- Steps 06–09 (Environment 2): 44/44 ([report](verify/report_steps06-09.csv)).

## Analysis notes
- **Repeat biopsies** (PD161/PD161b, PD164/PD164b/PD164c) are separate samples and separate integration batches throughout. Donor is a blocking factor in DA and DE.
- **Cluster 9 removal does not drive the C9 result.** The donor-blocked C9 p-value is 0.0071 as published, 0.0055 with cluster 9 kept, 0.0059 with it merged into C2, and 0.0048 with no cleanup at all.
- **ρ = 0.99 is kept in the ρ-sensitivity sweep** (`DA_propeller_rho_sensitivity*.csv`). It is a numerical boundary and the one grid point where the result does not hold.
- **Immunoglobulin genes are flagged (`is_ig`) and left out of the volcano plots,** because ambient Ig is donor-driven in this dataset.

## Authors
- Step 00 (QC and integration): Tyler Chinsky
- Steps 01–09: Apoorva Sharma
