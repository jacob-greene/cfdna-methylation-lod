#!/usr/bin/env bash
#SBATCH --job-name=p3b
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%A_%a.log
# Phase 3b. STEP=prepare runs scripts/phase3b_prepare.py once. Otherwise one array task per
# line of TASKS (track features hi1); extra classifier flags come from ARGS.
# Usage: REPO_DIR=<repo> CFDNA_DATA=<data dir> STEP=prepare sbatch phase3b.sbatch
#        REPO_DIR=<repo> CFDNA_DATA=<data dir> sbatch --array=0-15 phase3b.sbatch
set -euo pipefail
: "${REPO_DIR:?}" "${CFDNA_DATA:?}"
cd "$REPO_DIR"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
if [[ "${STEP:-classify}" == prepare ]]; then
  exec uv run --offline python -u scripts/phase3b_prepare.py
fi
TASKS=()
for hi1 in without with; do
  for feat in atlas tumour; do
    for track in array 30 100 300; do TASKS+=("$track $feat $hi1"); done
  done
done
read -r track feat hi1 <<<"${TASKS[$SLURM_ARRAY_TASK_ID]}"
# shellcheck disable=SC2086
uv run --offline python -u scripts/phase3b_classify.py --track "$track" --features "$feat" \
  --hi1 "$hi1" --out "${OUT:-results/phase3b/raw}" ${ARGS:-}
