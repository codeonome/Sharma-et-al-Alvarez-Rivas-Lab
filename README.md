# Sharma et al. — Alvarez-Rivas Lab

Analysis code for Sharma et al., *Clonally Convergent Tph/Tfh-like CD4+ T Cells are Shared Across Recovery and Recurrence in Experimental Periodontitis*: single-cell RNA and TCR sequencing of gingival T cells in periodontitis, in mouse and human.

| Folder | Contents |
|---|---|
| [`Mouse-data-analysis/`](Mouse-data-analysis/) | Sensitized Recurrent Periodontitis (SeRPD) mouse model. 10x 5′ scRNA-seq + paired V(D)J at baseline, periodontitis, recovery and recurrence. Global, myeloid and T-cell analysis; clonal expansion; TCR repertoire; cross-species Tph comparison. |
| [`Human-data-analysis/`](Human-data-analysis/) | Reanalysis of public human gingiva scRNA-seq (GSE164241; 21 biopsies, 18 donors, healthy vs periodontitis). Covers T-cell subclustering (C1–C9), differential abundance and pseudobulk differential expression. Includes the published results and a verified, step-wise reproducible pipeline. |

Each folder has its own README describing:
- the steps and how to run them
- the computational environment (`renv.lock` / `environment/`)
- the mapping from manuscript figures to code

## Data
- **Raw sequencing data:**
  - Mouse: [GSE350875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350875)
  - Human: [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241) (Williams et al., 2021)
- **Seurat objects:** not stored in git because of their size; deposited separately *[Zenodo DOI — to be added]*.

## Cross-species link
`Human-data-analysis` produces the final human T-cell object and the curated human–mouse Tph gene table, which `Mouse-data-analysis/06_CrossSpecies_Tph` uses to compare mouse Tfh/Tph-like cells with human Tph-like cells (C9).
