# Computational environment

The pipeline runs in two environments, exactly as the published analysis did. Both were captured on 2026-10-15, from the final verified run, with `env_capture.R`, scanning the pipeline scripts, and both were verified to regenerate the published results.

| Steps | Environment | Folder |
|---|---|---|
| 00–05 (integration, clustering, subsetting) | UB CCR, `module load foss/2021b r/4.2.0` (ccrsoft/2023.01) | `ccr_steps00-05/` |
| 06–09 (cleanup, DA, figures, DE) | local workstation, R 4.5.1 (macOS, Apple M4 Pro) | `local_steps06-09/` |

## Files in each folder
| File | Contents |
|---|---|
| `system.txt` | OS, CPU, R, BLAS/LAPACK, RNG kind, loaded modules, library paths, packages called |
| `sessionInfo.txt` | `sessionInfo()` with every package the scripts use attached |
| `packages_used.csv` | packages the scripts call plus all their recursive dependencies (version, library, build) |
| `packages_installed.csv` | every package installed in that environment |
| `renv.lock` | lockfile for the packages used |
| `<step>_sessionInfo.txt` | `sessionInfo()` written by each step of the verified run |

The CCR folder also has `module_list.txt`, `module_show_r-4.2.0.txt` and `flexiblas_config.txt`.

## Key versions
| | CCR (00–05) | Local (06–09) |
|---|---|---|
| OS | Ubuntu 24.04.5, x86_64 | macOS 15.7.4, arm64 |
| R | **4.2.0** | **4.5.1** |
| BLAS / LAPACK | FlexiBLAS 3.0.4 → **OpenBLAS 0.3.18** | R reference BLAS/LAPACK |
| Bioconductor | 3.16 | 3.22 |
| Seurat / SeuratObject | **5.1.0 / 5.0.2** | 5.4.0 / 5.3.0 |
| harmony | **0.1.0** | not used |
| uwot / RcppAnnoy / irlba | 0.2.2 / 0.0.22 / 2.3.5 | not used |
| Matrix | 1.6-5 | 1.7-4 |
| presto (used implicitly by `FindAllMarkers`) | 1.0.0 (user library) | not used |
| limma / edgeR / speckle | not used | **3.66.0 / 4.8.2 / 1.10.0** |
| ggplot2 / dplyr / tidyr | 3.5.0 / 1.0.8 / 1.2.0 | 4.0.1 / 1.1.4 / 1.3.2 |

## Why two environments
- **Harmony version decides the clustering.** On the local machine (harmony 1.2.4) the immune and T-cell embeddings diverge from the published ones. In the CCR stack (harmony 0.1.0) they are identical.
- **The CCR stack is on read-only cvmfs** and has not changed since the published run.
- **limma, edgeR and speckle are not in the CCR stack,** so DA and DE run locally.

## Recreating
**CCR:**
```bash
module load foss/2021b r/4.2.0
Rscript -e 'install.packages("presto", repos = NULL, type = "source")'   # or remotes::install_github("immunogenomics/presto")
```
The original session was an OnDemand RStudio server with 16 cores. FlexiBLAS's default backend is OpenBLAS.

**Local:**
```r
# R 4.5.1
renv::restore(lockfile = "environment/local_steps06-09/renv.lock")
```

## Step 00 environment
Step 00 was run by Tyler Chinsky (CCR, 2026-08-28) with the same CCR stack. The deposited integrated object is the verified starting point for steps 01–09.
