# Module 04 — T Cell Clonal Expansion Analysis

## What this module does

Takes the fully annotated αβ T cell object from Module 03 and adds clonal expansion status to each cell using TCR metadata from scRepertoire. Produces summary tables and publication figures characterising where expanded clones sit across cell types and timepoints.

Input: `../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds` (9,224 cells — MAIT/innate-like cells physically removed; see Module 03 `Tcell_remove_MAIT.R`)

---

## Script — `Tcell_expansion_analysis.R`

### Key parameters

| Parameter | Value | Notes |
|---|---|---|
| Clonotype identity | `CTstrict` | Paired α/β V+J gene + CDR3 nucleotide — most specific |
| Clone size unit | per `orig.ident` | Prevents cross-sample inflation |
| Expansion threshold | ≥ 2 cells | Singletons = Non-Expanded |
| Aggregation | sample-averaged | % per sample → mean across samples per condition |

### Steps

1. **STAGE 1** — Load object, compute `clone_size_in_sample` (CTstrict × orig.ident counts), assign `expansion_status`
2. **STAGE 2** — Export summary CSVs
3. **STAGE 3** — Produce all figures (SVG)
4. **STAGE 4** — Save annotated object

---

## Outputs

### Objects

| File | Description |
|---|---|
| `objects/tcell_TCRab_expansion_annotated05.rds` | Annotated object with expansion_status added to metadata |
| `objects/cells_expansion_table.csv` | One row per cell: orig.ident, condition, CTstrict, clone_size_in_sample, expansion_status |
| `objects/clonotype_by_sample_summary.csv` | One row per sample × clonotype: clone size per sample |
| `objects/clonotype_global_summary.csv` | One row per clonotype: global frequency, n samples, n conditions |
| `objects/expansion_BAS/PD/REC/RCU.csv` | % Expanded within each cell type per condition (denominator = cells of that type) |
| `objects/expanded_pool_composition.csv` | Cluster composition of the **expanded pool** per condition — `pct_pooled` (report this) + `pct_sampavg` (sensitivity check). Columns sum to 100% |
| `objects/expanded_pool_denominators.csv` | n expanded cells per condition (BAS 31, PD 750, REC 748, RCU 403) |

### ⚠ Three denominators — do not mix them

The same cells give very different percentages depending on the denominator:

| Output | Denominator | REC Tfh/Tph-like |
|---|---|---|
| `expanded_pool_composition.csv` | all **expanded** cells in the condition | **42.2%** (pooled) |
| alluvial figures | all T cells per sample | **8.41%** (sample-avg) |
| `expansion_REC.csv` | all **Tfh/Tph-like** cells per sample | **44.9%** (sample-avg) |

**The alluvial band heights answer "what % of all T cells is this?" — they do not answer
"what % of the expanded pool is this cell type?"** For REC Tfh/Tph those are 8.41% and
42.2% respectively. The figure gives you the first; `expanded_pool_composition.csv` gives
you the second.

The manuscript's "42% Th1/Tfh-like at recovery" is the pooled expanded-pool composition.
State the denominator *and* the aggregation wherever any of these numbers appear.

### Figures

| File | Description |
|---|---|
| `figures/expansion_UMAP_all.svg` | UMAP coloured by expansion status (all conditions combined) |
| `figures/expansion_UMAP_byCondition.svg` | UMAP faceted by condition (4 panels) |
| `figures/expansion_bar_pct.svg` | % Expanded vs Non-Expanded per condition (stacked, sample-averaged) with % labels |
| `figures/expansion_bar_counts.svg` | Mean cell counts per condition (stacked) |
| `figures/expansion_alluvial_expanded.svg` | Expanded: cell type composition (% of all T cells) |
| `figures/expansion_alluvial_expanded_counts.svg` | Expanded: cell type composition (mean cell count) |
| `figures/expansion_alluvial_nonexpanded.svg` | Non-Expanded: cell type composition (% of all T cells) |
| `figures/expansion_alluvial_nonexpanded_counts.svg` | Non-Expanded: cell type composition (mean cell count) |
| `figures/expansion_pct_heatmap.svg` | % Expanded within each cell type × condition (tile heatmap) |
| `figures/expansion_combined.svg` | Combined publication figure (A–F, 18×15 in) |

### Combined figure panel layout

| Panel | Content | Size in combined figure |
|---|---|---|
| A | % Expanded bar chart (sample-averaged, with % labels) | 6 × 5.3 in |
| B | Mean cell count bar chart (no legend) | 3 × 5.3 in |
| C | UMAP all conditions (no axis lines) | 9 × 5.3 in |
| D | UMAP faceted by condition (no axis lines) | 18 × 4.5 in |
| E | Alluvial — Expanded clones, % of T cells (no legend) | 9 × 5.3 in |
| F | Alluvial — Expanded clones, mean cell count | 9 × 5.3 in |
