#!/bin/bash
# ---------------------------------------------------------------------
# Steps 00-05 in Environment 1 (R 4.2.0, Seurat 5.1.0, harmony 0.1.0;
# see environment/env1_steps00-05/). Run from Human-data-analysis/:
#
#   bash run/run_steps00-05.sh                     # steps 01-05, starting from
#                                                  # output/objects/00_integrated.rds
#   STEPS="00 01 02 03 04 05" bash run/run_steps00-05.sh   # from raw Cell Ranger output
#
# Steps 02-03 need ~100 GB RAM. On an HPC cluster, submit this script with
# your scheduler (e.g. `sbatch run/run_steps00-05.sh`) after loading R.
# Optional: PIPELINE_OUT (default output), RAW_10X_DIR (step 00 only).
# ---------------------------------------------------------------------
set -euo pipefail
STEPS="${STEPS:-01 02 03 04 05}"
OUT="${PIPELINE_OUT:-output}"
mkdir -p "$OUT/logs"
Rscript -e 'v <- as.character(packageVersion("harmony")); if (v != "0.1.0") message("WARNING: harmony ", v, " found; the published clusters were produced with harmony 0.1.0 and differ with harmony >= 1.0")'

declare -A SCRIPT=(
  [00]=00_integration.R [01]=01_qc_integrated_object.R [02]=02_global_clustering.R
  [03]=03_immune_subset.R [04]=04_tcell_reclustering.R [05]=05_tcell_validation.R
)
for s in $STEPS; do
  echo "=== step $s: ${SCRIPT[$s]} | $(date) ==="
  Rscript --no-save --no-restore "scripts/${SCRIPT[$s]}" 2>&1 | tee "$OUT/logs/step${s}.log"
done
echo "=== done | $(date) ==="
