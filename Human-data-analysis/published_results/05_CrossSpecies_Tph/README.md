# Module 05 — Cross-species Tph handoff

A curated human/mouse gene table used to ask whether **C9 "CD4 Tph-like (CXCL13+)"** from
[Module 04](../04_Tcell_Final/) and the `Tfh/Tph-like` cluster from the companion mouse study
represent a conserved transcriptional program.

**No analysis runs in this module.** It holds the curated input; the analysis lives in the mouse
repository, which reads this table and the Module 04 Seurat object.

---

## Contents

`human_mouse_Tph_chronic_activation_gene_table.csv` — 81 hand-curated human/mouse gene pairs.

| Column | |
|---|---|
| `human_gene` / `mouse_gene` | ortholog pair |
| `category` | see below |
| `specificity_flag` | `high` (34), `medium` (39), `low` (8) |
| `note` | why the gene is in the set |

| Category | n | |
|---|---|---|
| `Tph_mouse_specific` | 26 | from the mouse Tph marker list only |
| `Tph_human_specific` | 22 | from the human C9 marker list only |
| `chronic_activation` | 19 | exhaustion / chronic-activation program |
| `cross_species_anchor` | 13 | direct overlap between both lists |
| `Tph_human_specific_and_chronic_activation` | 1 | both |

---

## Where the analysis lives

```
Alvarez-lab/Clean_code_Version01/Publication_Ready_Project-v2/07_CrossSpecies_Tph/
├── scripts/01_marker_overlap_specificity.R    hypergeometric marker overlap
├── scripts/02_crossspecies_module_score.R     reciprocal module scoring  ← reads this table
├── scripts/03_metaneighbor.R                  MetaNeighbor AUROC replicability
└── scripts/04_bhelper_featureumap.R           B-helper feature UMAPs
```

That module expects a `Human/` folder holding
`tcell_final_seurat_object.rds`, `tcell_subset_markers.csv`, `07_tcell_final_v4.R` and this gene
table — all of which live here in Modules 04 and 05.

Its target clusters are human `cluster_code == "C9"` (873 cells) and mouse
`celltype_merged == "Tfh/Tph-like"` (1,045 cells).

---

## Result

**The three tests do not support a conserved Tph cell type.** MetaNeighbor — purpose-built for
this question — scores C9 ↔ mouse Tfh/Tph-like at **AUROC 0.66**, not a reciprocal top hit, and
0.43 (below chance) under the stringent one-vs-best mode. The same run recovers Treg↔Treg at 0.86
and two CD8↔CD8 pairs at 0.82 / 0.77, so the method demonstrably works on this data. C9's nearest
mouse relative is Th2-like (0.76); the mouse Tph cluster's nearest human relative is C4/Th17
(0.73).

The one analysis that looks positive — the 13-gene `cross_species_anchor` set — **is circular by
construction**, since those genes were selected *because* they overlap. It cannot be used as
evidence. The `specificity_flag` column exists to make that circularity auditable.

See the mouse module's own README for the full argument.

---

## Note on placement

This table sits in the human repository because it was curated alongside the human C9 annotation,
but it is consumed entirely by the mouse repository. If the two are ever published as one
deposit, this module folds into the mouse module 07.

**Related:** [Module 04 — Final T-cell analysis](../04_Tcell_Final/)
