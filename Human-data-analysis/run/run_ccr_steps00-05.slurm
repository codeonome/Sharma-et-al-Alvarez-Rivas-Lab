#!/bin/bash
#SBATCH --job-name=gingiva_00-05
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --account=carlaalv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/gingiva_00-05_%j.log
# ---------------------------------------------------------------------
# Steps 00-05 on UB CCR, in the environment of the published run
# (environment/ccr_steps00-05/). Submit from the pipeline/ directory:
#
#   sbatch run/run_ccr_steps00-05.slurm          # steps 01-05 (start from the
#                                                 # deposited integrated object)
#   STEPS="00 01 02 03 04 05" sbatch run/run_ccr_steps00-05.slurm
#
# Before steps 01-05, put the deposited integrated object at
#   $PIPELINE_OUT/objects/00_integrated.rds
# Optional: PIPELINE_OUT (default output), RAW_10X_DIR (step 00 only).
# ---------------------------------------------------------------------
set -euo pipefail
STEPS="${STEPS:-01 02 03 04 05}"
module load foss/2021b r/4.2.0     # ccrsoft/2023.01: R 4.2.0, Seurat 5.1.0, harmony 0.1.0
mkdir -p logs "${PIPELINE_OUT:-output}/logs"
module list 2>&1 | tee "${PIPELINE_OUT:-output}/logs/ccr_module_list.txt"

declare -A SCRIPT=(
  [00]=00_integration.R [01]=01_qc_integrated_object.R [02]=02_global_clustering.R
  [03]=03_immune_subset.R [04]=04_tcell_reclustering.R [05]=05_tcell_validation.R
)
for s in $STEPS; do
  echo "=== step $s: ${SCRIPT[$s]} | $(date) ==="
  Rscript --no-save --no-restore "scripts/${SCRIPT[$s]}" 2>&1 | tee "${PIPELINE_OUT:-output}/logs/step${s}.log"
done
echo "=== done | $(date) ==="
