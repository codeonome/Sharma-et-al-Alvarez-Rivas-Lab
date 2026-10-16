# Computational environment

The pipeline uses two R environments. Both were recorded from the final run, which was verified to regenerate the published results exactly.

| Steps | Environment | Record |
|---|---|---|
| 00–05 (integration, clustering, subsetting) | **Environment 1:** R 4.2.0, Seurat 5.1.0, harmony 0.1.0 (Linux x86_64, OpenBLAS) | [`env1_steps00-05/`](env1_steps00-05/) |
| 06–09 (cleanup, DA, figures, DE) | **Environment 2:** R 4.5.1, limma 3.66.0, edgeR 4.8.2, speckle 1.10.0 | [`env2_steps06-09/`](env2_steps06-09/) |

<details>
<summary><b>Key package versions</b></summary>

| | Environment 1 (00–05) | Environment 2 (06–09) |
|---|---|---|
| R | **4.2.0** | **4.5.1** |
| Bioconductor | 3.16 | 3.22 |
| BLAS / LAPACK | OpenBLAS 0.3.18 (via FlexiBLAS 3.0.4) | R reference BLAS/LAPACK |
| Seurat / SeuratObject | **5.1.0 / 5.0.2** | 5.4.0 / 5.3.0 |
| harmony | **0.1.0** | not used |
| uwot / RcppAnnoy / irlba | 0.2.2 / 0.0.22 / 2.3.5 | not used |
| Matrix | 1.6-5 | 1.7-4 |
| presto (used automatically by `FindAllMarkers`) | 1.0.0 | not used |
| limma / edgeR / speckle | not used | **3.66.0 / 4.8.2 / 1.10.0** |
| ggplot2 / dplyr / tidyr | 3.5.0 / 1.0.8 / 1.2.0 | 4.0.1 / 1.1.4 / 1.3.2 |

</details>

## Recreating the environments
**Environment 1** (R 4.2.0):
```r
renv::restore(lockfile = "environment/env1_steps00-05/renv.lock")
remotes::install_github("immunogenomics/presto", ref = "a24772a135c7895a8183b007376050556c60a05b")
```
**Environment 2** (R 4.5.1):
```r
renv::restore(lockfile = "environment/env2_steps06-09/renv.lock")
```

## Notes
- **harmony 0.1.0 is what reproduces the published clusters.** Harmony ≥ 1.0 is a rewrite and gives a different embedding, so the immune and T-cell clusters change. The cell-type conclusions hold, but cluster membership is not identical.
- **presto is an implicit dependency.** Seurat's `FindAllMarkers` uses it automatically when it is installed, so install it in Environment 1 to reproduce the marker tables exactly.
- **Exact reproduction was verified in the recorded environments.** On a different platform or BLAS, numerical results may differ in the last digits.
- **Environment 1 was originally an EasyBuild `foss/2021b` module stack.** The loaded modules are recorded in `env1_steps00-05/module_list.txt`.
- **Step 00 (Tyler Chinsky) used Environment 1.** The deposited integrated object is the verified starting point for steps 01–09.

<details>
<summary><b>Files in each environment folder</b></summary>

| File | Contents |
|---|---|
| `system.txt` | OS, CPU, R, BLAS/LAPACK, RNG kind, library paths, the packages the scripts call |
| `sessionInfo.txt` | `sessionInfo()` with every package the scripts use attached |
| `<step>_sessionInfo.txt` | `sessionInfo()` written by each step of the verified run |
| `packages_used.csv` | the packages the scripts call plus all their dependencies, with version and build |
| `packages_installed.csv` | every package installed in that environment |
| `renv.lock` | lockfile for the packages used |

`env1_steps00-05/` also has `module_list.txt`, `module_show_r-4.2.0.txt` and `flexiblas_config.txt`.

</details>
