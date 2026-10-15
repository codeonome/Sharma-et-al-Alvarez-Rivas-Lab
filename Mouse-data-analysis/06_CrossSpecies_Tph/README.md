# Module 06 — Cross-species Tph comparison

Do the mouse `Tfh/Tph-like` cluster and the human `C9 "CD4 Tph-like (CXCL13+)"` cluster,
derived independently, represent a conserved transcriptional program?

**Short answer: no — not as a conserved cell type.** Three independent tests agree:
marker-list overlap, reciprocal module scoring, and MetaNeighbor. The one analysis
that appears positive (the 13-gene anchor set) is circular by construction and
cannot be used as evidence.

The decisive number is MetaNeighbor, which is purpose-built for this question and
demonstrably works on this data — it recovers Treg↔Treg at AUROC 0.86 and two
CD8↔CD8 pairs at 0.82 / 0.77. In that same run, **C9 ↔ mouse Tfh/Tph-like scores
0.66 and is not a reciprocal top hit**; under the stringent one-vs-best mode it is
0.43, below chance. C9's actual nearest mouse relative is Th2-like (0.76), and the
mouse Tph cluster's nearest human relative is C4/Th17 (0.73).

**In the paper:** Figure 7 (B-helper module heatmaps, feature UMAPs and violins) comes
from script 04. Scripts 01–03 are the formal cross-species tests behind the statement that
the two clusters are not the same cell type.

---

## How to run

`./run_all.sh 06` from the repository root, or:

```bash
cd 06_CrossSpecies_Tph/scripts
Rscript 01_marker_overlap_specificity.R     # ~10 min (recomputes both marker tables)
Rscript 02_crossspecies_module_score.R      # ~8 min
Rscript 03_metaneighbor.R                   # ~12 min
Rscript 04_bhelper_featureumap.R            # ~4 min  (independent of 01-03)
```

Script 01 caches `../objects/markers_matched_*.csv` and reuses them if present;
delete those two files to force a recompute. Script 02 is independent. **Script 03
requires script 01** to have run, for `objects/ortholog_map_1to1.csv`.

`MetaNeighbor` (Bioconductor, v1.29.0) is recorded in the root `renv.lock` and installed by `renv::restore()`.

## Inputs

| Role | Path |
|---|---|
| Curated gene table (81 ortholog pairs) | `../Human/human_mouse_Tph_chronic_activation_gene_table.csv` |
| Human object (20,802 cells) | `../Human/tcell_final_seurat_object.rds` |
| Human markers | `../Human/tcell_subset_markers.csv` |
| Mouse object (9,224 cells) | `../03_Tcell_Analysis/objects/tcell_TCRab_annotated_final04_clean.rds` |
| Mouse markers | `../03_Tcell_Analysis/objects/markers_TCRab_res0.8_top500.csv` |

Target clusters: human `cluster_code == "C9"` (873 cells; = raw cluster **8**, before
`Human-data-analysis` step 06 renames 0–8 to C1–C9) and mouse `celltype_merged == "Tfh/Tph-like"`
(1,045 cells; = `RNA_snn_res.0.8` cluster 2).

---

## Result 1 — Marker overlap does not single out the Tph pair

The hypergeometric answer depends almost entirely on the background universe N,
so three are reported. Two marker-table versions are also reported, because the
published tables used different `min.pct` (human 0.25, mouse 0.1) and the mouse
one was truncated to the top 500 per cluster.

C9 × mouse Tfh/Tph-like:

| version | universe | N | k | expected | fold | p |
|---|---|---|---|---|---|---|
| published | permissive | 12,155 | 35 | 11.2 | 3.13× | 2.0e-09 |
| published | **testable** | 1,159 | 30 | 28.8 | **1.04×** | 0.43 |
| published | marker-space | 680 | 35 | 43.1 | 0.81× | 0.96 |
| matched | permissive | 12,155 | 47 | 10.6 | 4.45× | 3.3e-18 |
| matched | **testable** | 1,159 | 47 | 49.3 | **0.95×** | 0.69 |
| matched | marker-space | 787 | 47 | 54.6 | 0.86× | 0.93 |

**Use the `testable` universe.** It is the only one where both positive controls
are recovered in both versions — C2 CD4-naive → mouse Naive T-like, and C4 Th17 →
mouse Th17-like. The permissive universe fails the Th17 control in both versions;
marker-space fails it in the matched version. Full record in
`objects/reciprocal_best_hits.csv`.

In that validated universe the Tph overlap is **exactly chance** (1.04× and 0.95×).
The 3–4× seen in the permissive universe is an artifact of testing against 12,155
genes when subset markers can only be drawn from ~1,000 — there, every one of the
99 cluster pairs scores 2–7×.

Neither Tph cluster picks the other as its best match in any configuration
(`objects/reciprocal_best_hits.csv`).

## Result 2 — Reciprocal module scoring agrees

`AddModuleScore`, `full` variant. "verdict" = clusters where the target scores
significantly higher (Dunn BH < 0.05) out of all contrasts; delta = median Cliff's delta.

| species | gene set | role | verdict | delta | peaks at |
|---|---|---|---|---|---|
| human | ANCHOR | shared (**circular**) | 8/8 | 0.72 | C9 |
| human | TPH_HUMAN | within-species control | 8/8 | 0.60 | C9 ✓ |
| human | **TPH_MOUSE** | **cross-species probe** | **2/8** | **−0.16** | **C3** ✗ |
| human | CHRONIC | orthogonal | 2/8 | −0.06 | C8 (MAIT) |
| mouse | ANCHOR | shared (**circular**) | 9/9 | 0.67 | Tfh/Tph-like |
| mouse | TPH_MOUSE | within-species control | 9/9 | 0.81 | Tfh/Tph-like ✓ |
| mouse | **TPH_HUMAN** | **cross-species probe** | **5/9** | **0.06** | **Tregs-like** ✗ |
| mouse | CHRONIC | orthogonal | 8/9 | 0.54 | Tfh/Tph-like |

Both within-species positive controls pass, so the pipeline works. **Both
cross-species probes fail**: the mouse Tph signature peaks at human C3
(CD8 TEM ZEB2-hi) with a *negative* effect size at C9, and the human Tph signature
peaks at mouse Tregs-like. UCell reproduces every one of these (Spearman ρ 0.85–0.96
against AddModuleScore, `objects/score_method_concordance.csv`), and the human
donor-level check reproduces them too (`objects/donor_level_human.csv`).

### The anchor set is circular — do not cite it as evidence

The 13 anchor genes were *defined* as the genes appearing in both clusters' marker
lists. Scoring them back onto those same clusters is guaranteed to succeed; the
8/8 and 9/9 with delta ≈ 0.7 is a pipeline sanity check, not a conservation result.

`objects/anchor_gene_specificity.csv` shows why it also fails on its own terms:
only **2 of 13** anchors (ARAP2, DOCK8) are exclusive to C9 among the nine human
clusters. RORA and ETV6 each mark 5 of 9, STAT3 and LDLRAD4 mark 4. The set is
largely a shared T-cell activation program, not a Tph identity program.

## Result 3 — Chronic activation is a mouse-Tph feature, not a shared one

Mouse Tfh/Tph-like is the top-scoring cluster for the chronic-activation set
(8/9, delta 0.54). Human C9 is not (2/8, delta −0.06; the human peak is C8/MAIT).

Within human C9, Control vs PD (`objects/chronic_by_condition_in_tph.csv`):

| level | n | p | verdict |
|---|---|---|---|
| cell-level | 257 vs 616 | 7.9e-05 | significant, but pseudo-replicated |
| donor-mean | 7 vs 4 | 0.23 | **not significant** |
| donor-mean, PD164 dropped | 7 vs 3 | 0.38 | **not significant** |

The cell-level result does not survive collapsing to donors. PD164 alone supplies
343 of the 616 PD cells in C9 (56%), so the donor-level result is the one to quote.
Mouse across BAS/PD/REC/RCU is Kruskal p = 3.9e-17 (but only 2 samples per condition —
see the limitation below).

---

## Result 4 — MetaNeighbor: the arbiter, and it agrees

MetaNeighbor (Crow et al. 2018) asks, over the whole transcriptome on 1:1 orthologs,
how well cells of a cluster in one dataset are identified by their nearest neighbours
in the other. AUROC 0.5 = chance, >0.8 = good, >0.9 = strong.

**The method works on this data.** Best mouse hit for each human cluster:

| human | best mouse match | AUROC | sensible? |
|---|---|---|---|
| C5 Treg | Tregs-like | **0.861** | ✓ Treg ↔ Treg |
| C4 Th17 | Tregs-like (Th17-like 0.808, 2nd) | 0.829 | ~ see note |
| C1 CD8 TEM GZMK | Teff CD8-like | 0.823 | ✓ CD8 ↔ CD8 |
| C3 CD8 TEM ZEB2 | Teff/Tem-like CD8 | 0.773 | ✓ CD8 ↔ CD8 |
| **C9 CD4 Tph-like** | **Th2-like** | **0.761** | ✗ **not Tph** |
| C2 CD4 naive/CM | Naive T-like | 0.739 | ✓ naive ↔ naive |
| C6 CD4 resting TSHZ2+ | Naive T-like | 0.734 | ✓ resting ↔ naive |
| C7 CD8 cytotoxic | Teff CD8-like | 0.700 | ✓ CD8 ↔ CD8 |
| C8 MAIT | Prolif T | 0.657 | ✗ no counterpart exists (mouse MAIT were removed in Module 03) |

Six of nine human clusters find a biologically correct counterpart, and `topHits`
calls five reciprocal top hits. That is a validated instrument.

**In that context, the Tph pairing fails:**

```
C9 <-> mouse Tfh/Tph-like     one_vs_rest AUROC = 0.661   rank 5 of 10 (human side)
                              one_vs_best AUROC = 0.428   below chance
                              reciprocal top hit: NO
C9      -> nearest mouse relative:  Th2-like       (0.761)
mouse Tfh/Tph-like -> nearest human relative: C4 / Th17 (0.727), then C9 (0.661)
```

C9's row contains four higher values than the Tph cell (Th2-like 0.76, Th17-like 0.72,
Naive T-like 0.69, Tregs-like 0.69). Under `one_vs_best` — each cluster against its
single closest competitor, the mode Crow et al. recommend for deciding whether a cell
type replicates — the pair scores 0.428, i.e. worse than a coin flip.

*Note on C4:* it pairs marginally better with mouse Tregs-like (0.829) than with
Th17-like (0.808). That is consistent with the annotation history recorded in
the human cluster-3 validation (`../../Human-data-analysis/scripts/05_tcell_validation.R`) — C4 was first read as a Th17-like effector Treg
because TNFRSF18/GITR, TNFRSF4/OX40 and IL2RA/CD25 appear in its markers.

## B-helper module feature UMAPs (script 04)

A separate, simpler deliverable: one feature UMAP per species for a curated
B-helper module score, each beside a cluster-labelled reference UMAP of the same
embedding.

**This is a within-species readout, not a cross-species test.** The two lists share
11 genes but differ in 4 each, so the two scores are **not numerically comparable
between species** — each panel says only "the B-helper program marks this cluster in
this dataset". The divergence is forced by biology: CXCL13 is expressed in 0.0% of
mouse T cells, so it cannot be in the mouse list. For the cross-species question see
Results 1–4 above.

```
HUMAN (15)  PDCD1 ICOS MAF BATF IL21 CD200 SH2D1A TOX2 TNFRSF4 CD4
            CXCL13 IL6ST PASK AFF3 CD40LG
MOUSE (15)  Pdcd1 Icos Maf Batf Il21 Cd200 Sh2d1a Tox2 Tnfrsf4 Cd4
            Cd40lg Bcl6 Il21r Cd84 Ebi3
```

`CXCR3` was excluded from both lists. It was the only gene running opposite to the
target cluster in either species — 2.6% of human C9 cells vs 6.9% of human T cells
overall, and 10.5% of mouse Tfh/Tph-like vs 17.7% overall — a Th1-axis marker that
was pulling both scores down. After its removal **every remaining gene is enriched
in the target cluster in human (15/15), and 14/15 in mouse** (`Sh2d1a` is the one
exception, 18.9% in target vs 22.9% overall). See `objects/bhelper_gene_detection.csv`.

**The score peaks in the target cluster in both species:**

```
HUMAN  mean score            MOUSE  mean score
C9  CD4 Tph-like    0.316    Tfh/Tph-like        0.172
C4  Th17            0.175    Tregs-like          0.137
C5  Treg            0.089    Th2-like            0.041
C6  CD4 resting     0.036    Th17-like          -0.037
C2  CD4 naive/CM   -0.032    Trm/Exh-like CD8   -0.141
C3  CD8 TEM ZEB2   -0.102    Teff CD8-like      -0.206
C1  CD8 TEM GZMK   -0.118    Trm type mixed     -0.206
C8  MAIT           -0.127    Prolif T           -0.240
C7  CD8 cytotoxic  -0.180    Teff/Tem-like CD8  -0.240
                             Naive T-like       -0.322
```

Scoring is `AddModuleScore(ctrl = 100, nbin = 24, seed = 1234)` on the `data` layer,
the same parameters used in script 02. Figure colour is clipped to the 1st–99th
percentile with the lower ~40% of the ramp held neutral grey, so the hotspot reads
without a handful of extreme cells compressing the scale.

### Statistics — read the effect size, not the p-value

Kruskal–Wallis across clusters underflows double precision in both species
(`p < 2.2e-308`). That number is not evidence: the clusters were *defined* from the
same expression data the score is computed on, so a difference is guaranteed, and
cells within a donor are not independent so any cell-level p is inflated. The
defensible test collapses to one mean per donor/sample:

| species | replicates | target − rest | paired Wilcoxon |
|---|---|---|---|
| human | 18 donors | +0.274 | **p = 1.5e-05** |
| mouse | 8 samples | +0.282 | **p = 7.8e-03** |

Both hold. `objects/bhelper_replicate_level_test.csv`.

Per-cluster contrasts (Dunn BH + Cliff's delta) are in
`objects/bhelper_target_vs_cluster_stats.csv` and marked on the violin figures:

```
HUMAN  C9 vs         delta  effect       MOUSE  Tfh/Tph-like vs     delta  effect
C7 CD8 cytotoxic     0.924  large        Naive T-like             0.906  large
C8 MAIT              0.860  large        Teff/Tem-like CD8        0.783  large
C1 CD8 TEM GZMK      0.847  large        Prolif T                 0.761  large
C3 CD8 TEM ZEB2      0.836  large        Teff CD8-like            0.761  large
C2 CD4 naive/CM      0.726  large        Trm type mixed CD8+ CD4  0.703  large
C6 CD4 resting       0.581  large        Trm/Exh-like CD8         0.640  large
C5 Treg              0.477  large        Th17-like                0.449  medium
C4 Th17              0.296  small        Th2-like                 0.282  small
                                         Tregs-like               0.067  negligible (p = 0.47, ns)
```

**Two caveats this exposes.** Human C9 separates from every other cluster, but its
margin over C4/Th17 is only a *small* effect. In mouse the score does **not**
separate Tfh/Tph-like from Tregs-like at all — δ 0.067, the single non-significant
contrast in either species. The mouse B-helper peak should therefore be described as
shared between Tfh/Tph-like and Tregs-like, not as specific to Tfh/Tph-like.

### Which genes drive the score

`figures/bhelper_gene_heatmap_{human,mouse}.svg` — the 15 genes × clusters, fill =
mean expression z-scored within each gene, number = % of cells expressing, clusters
in module-score order with the target outlined. This is what shows whether the score
is carried by a few genes or spread across the set.

- **human**: CXCL13, PASK, AFF3 and CD200 are the genuinely C9-specific drivers.
  ICOS, BATF, MAF, PDCD1 and TNFRSF4 are shared with C4/C5 and so contribute little
  discrimination; SH2D1A actually peaks in C1 and TNFRSF4 in C5.
- **mouse**: Il21, Cd200, Cd4, Pdcd1, Bcl6 and Cd40lg drive Tfh/Tph-like, but
  Tnfrsf4, Ebi3, Tox2, Icos and Maf are all *higher in Tregs-like* — which is the
  mechanistic reason those two clusters are not separable on this score.

**The mouse UMAP is flipped vertically for display** (`flip_y = TRUE` negates UMAP_2).
This is a presentation-only reorientation; it changes no data, no score and no
coordinate distance. The human UMAP is unflipped.

The feature UMAPs are single-panel — point size 0.60 (human) and 0.85 (mouse), and
no cluster reference panel. Cluster identity is carried by the violin and gene
heatmap figures, which use the same cluster ordering.

| Output | Contents |
|---|---|
| `figures/bhelper_featureumap_human.svg` | score on `umap.tcell`, single panel |
| `figures/bhelper_featureumap_mouse.svg` | score on `umap.pca_TCRab`, single panel (**flipped vertically**) |
| `figures/bhelper_score_violin_human.svg` | score distribution per cluster, ordered by median |
| `figures/bhelper_score_violin_mouse.svg` | score distribution per cluster, ordered by median |
| `figures/bhelper_gene_heatmap_human.svg` | 15 genes × clusters, z-scored mean expression |
| `figures/bhelper_gene_heatmap_mouse.svg` | 15 genes × clusters, z-scored mean expression |
| `figures/bhelper_gene_heatmap_{human,mouse}_nolabels.svg` | same heatmaps, fill only — no % text in the tiles |
| `objects/bhelper_target_vs_cluster_stats.csv` | Dunn BH p + Cliff's delta, target vs each cluster |
| `objects/bhelper_replicate_level_test.csv` | donor/sample-level paired Wilcoxon |
| `objects/bhelper_gene_by_cluster_expression.csv` | mean expr, z-score, % expressing per gene per cluster |
| `objects/bhelper_score_by_cluster.csv` | n cells, mean/median score, rank, Kruskal–Wallis p per species |
| `objects/bhelper_gene_detection.csv` | per gene, % cells expressing overall vs target cluster |
| `objects/bhelper_scores_percell_{human,mouse}.csv.gz` | per-cell scores |

---

## Why the two clusters do not correspond

They share the label, not the markers behind it (% cells expressing):

```
                  human C9   mouse Tfh/Tph
CXCL13 / Cxcl13     32.3%        0.0%      diverge
TOX2   / Tox2       45.8%        3.3%
PDCD1  / Pdcd1       9.2%       61.2%      diverge the other way
TNFSF11/ Tnfsf11     1.9%       76.3%
CXCR6  / Cxcr6       1.5%       78.6%
ICOS   / Icos       65.6%       84.8%      agree
BATF   / Batf       44.0%       54.0%      agree
MAF    / Maf        50.9%       34.7%      agree
IL21   / Il21        4.8%        5.0%      agree
CXCR5  / Cxcr5       6.6%        0.2%      agree (both low = both non-Tfh)
```

Canonical human Tph is PD-1-high / CXCR5-negative / CXCL13-positive. Human C9 has
CXCL13 but not PD-1 — a caveat already recorded in the human cluster-8 validation (`../../Human-data-analysis/scripts/05_tcell_validation.R`).
Mouse cluster 2 has Pdcd1 but zero Cxcl13, and Module 03 records its identity as
"inconclusive". Each was labelled from a different half of the definition.

Part of this is genuine species divergence: mouse T cells essentially do not express
Cxcl13, so the defining human gene cannot appear in the mouse list however well the
cell types match. What the two do share — ICOS, BATF, MAF, JUNB, STAT3 — is generic
activation machinery common to many T-cell subsets.

## What can be claimed

- **Not supportable:** "human C9 and the mouse Tfh/Tph-like cluster are the same
  cell type" / "a conserved Tph program". Three independent tests say no — including
  MetaNeighbor, which is purpose-built for the question and demonstrably works here —
  and the one supporting analysis is circular.
- **Supportable:** both are B-helper-competent, CXCR5-low activated CD4 populations
  that share ICOS/BATF/MAF/IL21 expression, reached by different routes — human
  CXCL13-driven, mouse PD-1/RANKL-driven.
- **Probably the stronger story:** `Tnfsf11` (RANKL) in 76.3% of mouse Tph cells vs
  1.9% in human C9 — a mouse-specific, directly osteoimmune-relevant signal.

MetaNeighbor has now been run (Result 4) and settles it. CCA/anchor integration and
label transfer were deliberately not run: integration is designed to find shared
structure, so co-clustering would be close to circular, and MetaNeighbor already
answers the question without that problem.

## Limitations

1. The anchor set is circular (above) and cannot serve as evidence of conservation.
2. Marker-list overlap is a fragile instrument here — best-hit assignments move with
   `min.pct` and universe choice. It is reported as a specificity check, not a headline.
3. Human C9 is 873 cells across 18 donors, and only 11 donors contribute ≥10 C9 cells.
   The PD arm has 5 donors, one of which (PD164) supplies 56% of PD C9 cells.
4. Mouse has 8 samples, 2 per condition. Sample-level means are written
   (`objects/sample_level_means_mouse.csv`) but **no sample-level test is run** —
   n = 2 per group is not testable. All mouse condition p-values are cell-level and
   therefore pseudo-replicated.
5. The `strict` variant leaves ANCHOR and CHRONIC untouched — no gene in either
   carries a `low` specificity flag. It only perturbs TPH_HUMAN (−3 genes) and
   TPH_MOUSE (−5). An unchanged anchor score is not independent confirmation.
6. MetaNeighbor AUROCs are moderate across the whole matrix (max 0.861). That is
   normal for cross-species comparisons and is why the analysis is read as *relative*
   ranking plus reciprocal-top-hit status, not against an absolute 0.9 cutoff. The
   validated correct pairings (0.74–0.86) set the scale against which the Tph pair
   (0.66) is judged.
7. Human C8 (MAIT) has no mouse counterpart by construction — Module 03 removed the
   27 mouse MAIT/innate-like cells. Its 0.657 pairing with mouse Prolif T is a
   no-match artifact, not a finding.
8. Human clustering is Harmony-integrated across donors; mouse clustering is
   unintegrated PCA. Human keeps TRBC1/TRBC2/TRDC in its HVG set (and one V segment,
   TRBV28) whereas mouse excludes TCR V/J genes. Neither difference is large enough
   to explain the result, but both belong in a methods paragraph.

## Outputs

**`objects/`**

| File | Contents |
|---|---|
| `markers_matched_{human,mouse}.csv` | Both marker sets recomputed at identical parameters |
| `overlap_matrix_<version>_<universe>.csv` | 6 full matrices: k, expected, fold, p, p_BH, p_Bonferroni |
| `reciprocal_best_hits.csv` | Best match per cluster per configuration + positive-control audit |
| `universe_sizes.csv` | N for each version × universe |
| `anchor_gene_specificity.csv` | Per-anchor logFC/pct in both targets + how many clusters each marks |
| `gene_sets_parsed.csv`, `gene_set_coverage.csv` | Parsed table; genes requested vs scored per set |
| `module_scores_percell_{human,mouse}.csv.gz` | Per-cell scores, all sets × variants × methods |
| `kruskal_results.csv`, `dunn_results.csv` | Across-cluster tests (Dunn with both BH and Bonferroni) |
| `target_vs_rest_summary.csv`, `target_vs_rest_headline.csv` | Per-contrast direction, Cliff's delta, adjusted p |
| `score_peak_cluster.csv` | Which cluster each gene set peaks in |
| `score_method_concordance.csv` | AddModuleScore vs UCell Spearman ρ |
| `donor_level_human.csv` | Donor-mean paired Wilcoxon, C9 vs rest |
| `sample_level_means_mouse.csv` | Mouse per-sample means (descriptive only) |
| `chronic_by_condition_in_tph.csv`, `chronic_condition_dunn_mouse_*.csv` | Condition tests within the Tph cluster |
| `ortholog_map_1to1.csv` | The 12,155 strict 1:1 ortholog pairs (written by 01, read by 03) |
| `metaneighbor_auroc_one_vs_{rest,best}.csv` | Full cross-species AUROC matrices |
| `metaneighbor_top_hits.csv`, `metaneighbor_tph_summary.csv` | Reciprocal top hits; the Tph verdict beside its positive controls |
| `metaneighbor_variable_genes.txt` | Genes variable in both studies, used by MetaNeighbor |

**`figures/`** — `overlap_fold_heatmap.svg`, `violin_scores_by_cluster_{human,mouse}.svg`,
`featureplot_anchor_score.svg`, `heatmap_mean_score_by_cluster.svg`,
`dotplot_anchor_genes.svg`, `full_vs_strict_scatter.svg`, `chronic_by_condition_in_tph.svg`,
`metaneighbor_auroc_heatmap.svg`, `metaneighbor_one_vs_best.svg`

**`results/`** — `sessionInfo_01.txt`, `sessionInfo_02.txt`, `sessionInfo_03.txt`
