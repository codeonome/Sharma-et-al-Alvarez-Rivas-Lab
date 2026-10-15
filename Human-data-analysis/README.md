# T-cell analysis of human gingiva in periodontitis — reproducible pipeline

This is a reanalysis of **GSE164241**: single-cell RNA-seq of human gingiva, 21 biopsies from 18 donors (13 healthy, GM; 8 periodontitis, PD). It starts from the integrated object and ends with nine T-cell populations (C1–C9), their differential abundance between PD and GM, and pseudobulk differential expression.

This folder holds the clean, step-wise analysis code (`scripts/`, `R/`, `config.R`) and, in [`published_results/`](published_results/), the tables and figures reported in the manuscript. **Every published table and cluster assignment is regenerated exactly by this pipeline** (see [Verification](#verification)).

## Manuscript result
> **Human periodontitis is characterized by gingival enrichment of antigen-experienced and tissue-adapted CD4⁺ T cells**
>
> To examine the transcriptomic landscape of antigen-experienced CD4⁺ T cells in GM and PD, we analyzed a publicly available scRNA-seq dataset comprising 13 GM and 8 PD samples. In this, 2 of 18 donors contributed multiple biopsies (donor PD164, three samples; donor PD161, two samples) (GSE164241; Williams et al., 2021). We focused on the T cell compartment by performing sequential unsupervised clustering of immune cells (upon eliminating innate and other contaminants) and identifying T cells by the expression of canonical T cell markers.
>
> Nine transcriptionally distinct clusters were identified among all T cells (C1–C9), displaying similar dominance in GM and PD.
>
> Five clusters were identified as CD4 lineage (C2, C4, C5, C6 and C9), three as CD8 lineage (C1, C3 and C7) and one as MAIT (C8).

### Figure panels → pipeline outputs
| Panel | Content | Figure (`output/figures/`) | Data behind it (`output/results/`) |
|---|---|---|---|
| C | T-cell UMAP, all subjects; density in GM and PD | `08_figures/umap_tcell_final.svg`, `08_figures/umap_density_by_condition.svg` | `objects/06_tcell_final.rds` (`umap.tcell`, `cluster_code`) |
| D | Grouped marker dot plot (CD8 · MAIT · CD4) | `08_figures/dotplot_grouped.svg` | `objects/06_tcell_final.rds` |
| E | Cluster composition, GM and PD (pooled cells) | `08_figures/donut_GM.svg`, `donut_PD.svg` | `08_figures/donut_pooled_percentages.csv` |
| F | Pseudobulk DE, PD vs GM, all T cells; TOX2, BCL6, CXCR5, CXCL13 highlighted | `09_pseudobulk_de/DE_volcano_PD_vs_GM_allTcells.svg` | `09_pseudobulk_de/DE_PD_vs_GM_allTcells.csv` |
| G | Differential abundance: odds ratio (95% CI), FDR; donor-blocked propeller | `08_figures/DA_bar_blocked.svg` | `07_differential_abundance/DA_propeller_donor_blocked.csv`, `08_figures/DA_bar_blocked_values.csv` |

Key numbers, regenerated exactly by the pipeline:

| | C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 |
|---|---|---|---|---|---|---|---|---|---|
| % of T cells, GM | 31.9 | 25.3 | 8.8 | 7.6 | 5.6 | 5.9 | 6.1 | 5.9 | 2.8 |
| % of T cells, PD | 32.7 | 20.3 | 12.2 | 7.0 | 7.4 | 7.0 | 5.1 | 3.0 | 5.3 |
| Odds ratio PD vs GM | 0.83 | 0.79 | 1.23 | 1.42 | 1.27 | 1.18 | 0.83 | 0.99 | **2.26** |
| FDR (donor-blocked propeller) | 0.54 | 0.54 | 0.54 | 0.54 | 0.54 | 0.59 | 0.54 | 0.98 | **0.064** (p = 0.0071) |

The within-donor correlation is ρ = 0.726 (DA) and ρ = 0.29 (pooled DE). Pooled DE gives 616 genes with |log2FC| > 1 at FDR < 0.05 (247 up in PD, 369 up in GM). TOX2, BCL6, CXCR5 and CXCL13 are all up in PD (FDR 0.001–0.03).

## Steps

| Step | Script | What it does | Runs on |
|---|---|---|---|
| 00 | `00_integration.R` | Per-sample QC, merge, CCA + Harmony integration, annotation (10x, scType, CloudAzimuth). *Author: Tyler Chinsky* | CCR |
| 01 | `01_qc_integrated_object.R` | Sanity checks on the integrated object | CCR |
| 02 | `02_global_clustering.R` | Resolution sweep on the Harmony embedding; res 0.1 = 12 clusters | CCR |
| 03 | `03_immune_subset.R` | Immune clusters → re-integration → 13 immune clusters → T cells (immune clusters 0, 1) | CCR |
| 04 | `04_tcell_reclustering.R` | T cells re-integrated and reclustered; res 0.4 = 11 clusters | CCR |
| 05 | `05_tcell_validation.R` | Validation of clusters 3 (Treg/Th17), 8 (Tph/Tfh), 9 (QC) | CCR |
| 06 | `06_tcell_cleanup_annotation.R` | Removes Ig-high doublets and cluster 9, merges 10 into 0, labels C1–C9 | local |
| 07 | `07_differential_abundance.R` | Composition; Wilcoxon and propeller tests; **primary: donor-blocked propeller**; leave-one-donor-out; ρ sensitivity | local |
| 08 | `08_figures.R` | UMAPs, density, donuts, DA figure, proportions, dot plot, ρ sensitivity | local |
| 09 | `09_pseudobulk_de.R` | Pseudobulk DE, PD vs GM: all T cells and each cluster C1–C9 | local |

Cell numbers: 110,166 integrated cells → 48,346 immune cells → 21,899 T cells → 20,802 T cells in C1–C9.

| Code | Population | Code | Population |
|---|---|---|---|
| C1 | CD8 TEM (GZMK/GZMH) | C6 | CD4 resting/CM (TSHZ2+) |
| C2 | CD4 naive/CM | C7 | CD8 cytotoxic (GZMB/GNLY) |
| C3 | CD8 TEM (ZEB2-hi) | C8 | MAIT |
| C4 | Th17 | C9 | CD4 Tph-like (CXCL13+) |
| C5 | Treg | | |

## Layout
```
pipeline/
├── config.R            all paths and parameters (the only file to edit)
├── R/helpers.R         shared functions (plotting, markers, DA models)
├── scripts/00-09_*.R   one script per step
├── run/                SLURM job for steps 00-05 (CCR), shell script for 06-09
├── verify/             comparison of outputs against the published results
├── environment/        R/package versions for both environments + renv.lock
├── published_results/  published tables + figures (modules 01–05), the verification reference
└── output/             (created by the run; git-ignored)
    ├── objects/        Seurat objects passed between steps
    ├── results/<step>/ tables
    ├── figures/<step>/ figures
    └── logs/           per-step logs and sessionInfo
```

## How to run
Every script is run from this folder and reads `config.R`. The output location can be changed with `PIPELINE_OUT`.

The analysis spans **two environments**, exactly as the published run did:

**1. Steps 00–05 on UB CCR** (`module load foss/2021b r/4.2.0`: R 4.2.0, Seurat 5.1.0, harmony 0.1.0)
```bash
# Start from the deposited integrated object (recommended):
mkdir -p output/objects
cp <deposited>/integrated_seurat_object.rds output/objects/00_integrated.rds
sbatch run/run_ccr_steps00-05.slurm                       # steps 01-05

# or from raw Cell Ranger output (step 00 needs internet for scType/CloudAzimuth):
RAW_10X_DIR=/path/to/cellranger STEPS="00 01 02 03 04 05" sbatch run/run_ccr_steps00-05.slurm
```
`presto` (1.0.0) must be installed in the CCR user library. Seurat's `FindAllMarkers` uses it automatically.

**2. Steps 06–09 locally** (R 4.5.1, limma 3.66.0, edgeR 4.8.2, speckle 1.10.0)
```bash
# copy output/objects/04_tcell_reclustered.rds from CCR, then
bash run/run_local_steps06-09.sh
```

**3. Verify**
```bash
Rscript verify/verify_outputs.R published_results            # all steps
Rscript verify/verify_outputs.R published_results 06 07 08 09
# Object checks also need the deposited .rds files under published_results/<module>/objects/
```

## Verification
The final pipeline was run from the deposited integrated object (md5 `a3926732393aa6a5e7387cf13565016e`) and every output was compared with the published files on 2026-10-15:

| Steps | Where | Checks | Result |
|---|---|---|---|
| 01 | CCR (job 26499400) | 9 QC checks | all pass |
| 02–05 | CCR (job 26499400, verified by job 26499401) | 12 tables; 4 objects: cells, every clustering resolution, PCA/Harmony/UMAP | **47 / 47 byte-identical** — [`verify/report_ccr_steps02-05.csv`](verify/report_ccr_steps02-05.csv) |
| 06–09 | local | 39 tables (cleanup, composition, all DA tests, all DE tables); final object: cells, labels, UMAP | **44 / 44 byte-identical** — [`verify/report_local_steps06-09.csv`](verify/report_local_steps06-09.csv) |

In total, **91 / 91 checks are identical**: same cells, the same 21,899 T cells, the same cluster at every resolution (ARI = 1), identical embeddings, and every published table byte for byte.

The two steps that are not re-verified are:
- **Step 00**, which needs raw data and network services.
- **The figures.** Figures depend on the graphics device and font setup; the tables and objects behind every figure are verified. Figures are written as SVG (the DA bar chart also as PDF). Four diagnostic plots are PNG only: the global UMAPs at res 0.1 and 0.2 (step 02), the T-marker violins (step 03) and the voom mean–variance plot (step 09).

## Changes from the original scripts
Every step whose logic changed the data or the analysis is kept, with identical results.

**Structure and code quality** (no effect on results):
- **Paths and parameters:** every path and parameter is in `config.R`; every script writes to its own output folder.
- **Metadata column names:**
  - Global clusters are now `global_res.<r>` instead of `RNA_snn_res.<r>`.
  - The original immune object relabelled inherited global columns as `immune_res.0.1/1/1.2/1.5`; here only true immune-level columns carry the `immune_res.` prefix.
- **Seeds:** the seeds the published run used are now passed explicitly.
- **Script 07 split into steps 06–09:** cleanup, DA, figures and DE.

**Code that produced published results but was not in a committed script, now included:**
- **Pseudobulk DE model:** `edgeR::voomLmFit` with sample quality weights and donor blocking; samples with fewer than 20 cells dropped. It produced all published DE tables, including the per-cluster tables.
- **T-cell marker panels:** `tcell_subset_marker_avg_expr.csv` and the panels, made interactively.

**Removed because not used and not retained in the published results:**
- **Dead code:** a duplicate `readRDS`, an undefined `saveRDS`, an unused subset and a duplicate object save.
- **Unused outputs:** a clustree plot and an Ig-filtered copy of the T-cell marker table.
- **Step 00 extras:** exploratory per-cell DE by cell type. SingleR was in the script but was never run.

## Analysis notes
- **Repeat biopsies** (GM184a, PD161/PD161b, PD164/PD164b/PD164c) are separate samples, and separate integration batches, throughout. Donor is modelled as a blocking factor in DA and DE.
- **Cluster 9 removal (step 06):** removed as ribosomal-gene-dominated with no lineage markers. The primary DA result is not sensitive to this choice. With donor-blocked propeller, C9 has p = 0.0071 as published, p = 0.0055 with cluster 9 kept as its own cluster, p = 0.0059 with it merged into C2, and p = 0.0048 with no cleanup at all.
- **ρ = 0.99 is kept in the ρ sweep (step 07).** It is a numerical boundary and the one grid point where the result does not hold. It is reported rather than removed.
- **Ig genes:** `is_ig` genes are excluded from the volcano plots, because ambient Ig is donor-driven in this dataset.

## Environment
See [`environment/ENVIRONMENT.md`](environment/ENVIRONMENT.md). For each environment it records the OS, BLAS, `sessionInfo()`, every package version and an `renv.lock`.

## Authors
- Step 00 (QC and integration): **Tyler Chinsky**
- Steps 01–09: Apoorva
