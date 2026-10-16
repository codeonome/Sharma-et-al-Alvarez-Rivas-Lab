#!/bin/bash
# ---------------------------------------------------------------------
# Steps 06-09 in Environment 2 (R 4.5.1, limma, edgeR, speckle;
# see environment/env2_steps06-09/). Run from Human-data-analysis/:
#   bash run/run_steps06-09.sh
# Needs, from steps 00-05 (or the deposited files):
#   $PIPELINE_OUT/objects/04_tcell_reclustered.rds
# ---------------------------------------------------------------------
set -euo pipefail
OUT="${PIPELINE_OUT:-output}"
mkdir -p "$OUT/logs"
for s in 06_tcell_cleanup_annotation 07_differential_abundance 08_figures 09_pseudobulk_de; do
  echo "=== $s | $(date) ==="
  Rscript --no-save --no-restore "scripts/$s.R" 2>&1 | tee "$OUT/logs/$s.log"
done
echo "=== done | $(date) ==="
