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
| Item | Content | Code |
|---|---|---|
| Fig. 1A–B | Imaging mass cytometry, human gingiva | not in this repository (QuPath/Fiji) |
| Fig. 1C–G | Human T-cell clusters, composition, pseudobulk DE, differential abundance | `Human-data-analysis/` (steps 07–09) |
| Fig. 2, Fig. 3 | Flow cytometry, micro-CT, qPCR, TRAP, adoptive transfer | not in this repository (FlowJo/GraphPad Prism) |
| Fig. 4B–E | Mouse immune landscape; T-cell subclusters by phase | `Mouse-data-analysis/` modules 01, 03 |
| Fig. 5A–E | Clonal expansion | `Mouse-data-analysis/` module 04 |
| Fig. 6B–F | TCRβ V–J usage, repertoire overlap, convergent clonotype | `Mouse-data-analysis/` module 05 |
| Fig. 7A–C | B-helper module, human vs mouse | `Mouse-data-analysis/` module 06 (uses outputs of `Human-data-analysis/`) |
| Supp. Fig. 1 | Flow cytometry gating and frequencies | not in this repository |
| Supp. Fig. 2–4 | Global clusters and TCR recovery; myeloid compartment; T-cell annotation | `Mouse-data-analysis/` modules 01, 02, 03, 05 |
| Supp. Table 1 | Human differential abundance, all clusters | `Human-data-analysis/` step 07 |
| Supp. Tables 2–4 | Jaccard per CD4 cluster; convergent clonotype variants | `Mouse-data-analysis/` module 05 |

Figures 4A and 6A are schematics.

## Data and code availability
- **Mouse scRNA-seq and TCR-seq:** GEO [GSE350875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350875)
- **Human scRNA-seq:** GEO [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241) (Williams et al., *Cell* 2021)
- **Processed Seurat objects:** too large for git; deposited at Zenodo *[DOI — to be added]*. Each README lists which objects go where.
- **Software versions:** `Mouse-data-analysis/renv.lock` and `Human-data-analysis/environment/`

## Contact
Dr. Carla Alvarez Rivas, Oral Biology Department, University at Buffalo (carlaalv@buffalo.edu)
