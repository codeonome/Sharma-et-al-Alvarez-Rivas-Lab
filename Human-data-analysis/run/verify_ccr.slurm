#!/bin/bash
#SBATCH --job-name=gingiva_verify
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --account=carlaalv
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/gingiva_verify_%j.log
# ---------------------------------------------------------------------
# Verify steps 02-05 against the published results (objects need ~40 GB RAM).
# Submit from pipeline/:   PUBLISHED_ROOT=/path/to/repo sbatch run/verify_ccr.slurm
# ---------------------------------------------------------------------
set -euo pipefail
module load foss/2021b r/4.2.0
mkdir -p logs
Rscript verify/verify_outputs.R "${PUBLISHED_ROOT:-..}" 02 03 04 05
