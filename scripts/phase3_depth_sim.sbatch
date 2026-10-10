#!/usr/bin/env bash
#SBATCH --job-name=p3-depth
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH --output=%x_%A_%a.log
# Phase 3 grid: one array task per line of TASKS (kind depth cap bg).
# Usage: REPO_DIR=<repo> CFDNA_DATA=<data dir> sbatch --array=0-42 phase3_depth_sim.sbatch
set -euo pipefail
: "${REPO_DIR:?}" "${CFDNA_DATA:?}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASKS=()
for d in 10 30 100 300 1000; do TASKS+=("cell $d none person" "plasma $d none person" "plasma $d none mean"); done
for d in 3000 10000; do TASKS+=("plasma $d none person" "plasma $d none mean"); done
for c in 30ng_r0.5 10ng_r0.5 10ng_r0.1; do
  for d in 300 1000 3000 10000; do TASKS+=("plasma $d $c person" "plasma $d $c mean"); done
done
read -r kind depth cap bg <<<"${TASKS[$SLURM_ARRAY_TASK_ID]}"
uv run python -u scripts/phase3_depth_sim.py --kind "$kind" --depth "$depth" --cap "$cap" \
  --bg "$bg" --reps "${REPS:-100}" --out results/phase3/raw
