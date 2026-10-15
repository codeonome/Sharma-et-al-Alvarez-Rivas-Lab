#!/usr/bin/env bash
# Run the full analysis pipeline, or selected modules, in dependency order.
#
#   ./run_all.sh                 # all modules, 01-06
#   ./run_all.sh 03 04           # only the listed modules (their inputs must already exist)
#
# Each script runs from its own scripts/ directory (all paths inside the scripts are
# relative to it). If the renv project library has been restored (renv::restore()
# creates renv/library/), it is activated for every script; otherwise the scripts use
# your default R library. Logs are written to <module>/results/logs/.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"

PIPELINE=(
  "01|01_Global_Analysis|Global_analyis_postprocessing.R"
  "02|02_Myeloid_Analysis|Myeloid_cells_analysis_refactored.R"
  "03|03_Tcell_Analysis|Tcell_analysis_reclustering.R"
  "03|03_Tcell_Analysis|Tcell_ab_analysis_reclustering.R"
  "03|03_Tcell_Analysis|Tcell_remove_MAIT.R"
  "03|03_Tcell_Analysis|Tcell_doublet_check.R"
  "03|03_Tcell_Analysis|Tcell_Trm_CD8_identity_check.R"
  "03|03_Tcell_Analysis|Tcell_Trm_CD8_feature_umap.R"
  "03|03_Tcell_Analysis|Tph_marker_violin.R"
  "03|03_Tcell_Analysis|assemble_figure_panels.R"
  "03|03_Tcell_Analysis|assemble_panel_D.R"
  "03|03_Tcell_Analysis|assemble_extended_fig.R"
  "04|04_Tcell_Expansion_Analysis|Tcell_expansion_analysis.R"
  "05|05_TCR_Repertoire_Analysis|TCR_clonotype_tracking.R"
  "05|05_TCR_Repertoire_Analysis|Beta_Gene_Prevalence_Analysis.R"
  "05|05_TCR_Repertoire_Analysis|Beta_Gene_Prevalence_Within_Condition.R"
  "05|05_TCR_Repertoire_Analysis|Beta_Gene_Jaccard_Global.R"
  "05|05_TCR_Repertoire_Analysis|Beta_Gene_Jaccard_By_CellType.R"
  "05|05_TCR_Repertoire_Analysis|TCR_VDJ_chord.R"
  "05|05_TCR_Repertoire_Analysis|TCR_Shared_Clonotype_CTstrict.R"
  "05|05_TCR_Repertoire_Analysis|TCR_Convergent_Clonotype_Summary.R"
  "05|05_TCR_Repertoire_Analysis|TCR_Convergent_Variants_Distribution.R"
  "05|05_TCR_Repertoire_Analysis|TCR_Convergent_Variant_Heatmaps.R"
  "05|05_TCR_Repertoire_Analysis|TCR_Dropout_Analysis.R"
  "06|06_CrossSpecies_Tph|01_marker_overlap_specificity.R"
  "06|06_CrossSpecies_Tph|02_crossspecies_module_score.R"
  "06|06_CrossSpecies_Tph|03_metaneighbor.R"
  "06|06_CrossSpecies_Tph|04_bhelper_featureumap.R"
)

want() {  # is module $1 selected? (no arguments = all modules)
  [ -z "${SELECTED:-}" ] && return 0
  for m in $SELECTED; do [ "$m" = "$1" ] && return 0; done
  return 1
}
SELECTED="$*"

for entry in "${PIPELINE[@]}"; do
  IFS='|' read -r num dir script <<< "$entry"
  want "$num" || continue
  mkdir -p "$ROOT/$dir/objects" "$ROOT/$dir/figures" "$ROOT/$dir/results/logs"
  log="$ROOT/$dir/results/logs/${script%.R}.log"
  echo "[$(date '+%H:%M:%S')] $dir/scripts/$script"
  (
    cd "$ROOT/$dir/scripts"
    Rscript -e "local({ o <- setwd('$ROOT'); if (length(Sys.glob('renv/library/*/R-*/*/Seurat')) > 0) source('renv/activate.R'); setwd(o) })" \
            -e "source('$script', echo = FALSE)"
  ) > "$log" 2>&1 || { echo "FAILED: $dir/scripts/$script — see $log"; exit 1; }
done
echo "[$(date '+%H:%M:%S')] done"
