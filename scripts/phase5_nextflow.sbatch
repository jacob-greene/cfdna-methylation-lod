#!/usr/bin/env bash
#SBATCH --job-name=p5-nxf
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
# Phase 5: the Nextflow head process. It submits every pipeline task as its own Slurm job.
# Usage: sbatch -p campus-new --export=ALL,REPO_DIR=<repo>,CFDNA_DATA=<data dir> \
#          -o <log> scripts/phase5_nextflow.sbatch [extra nextflow args, e.g. --only_runs SRR1]
# NXF_DIR sets the Nextflow launch directory (default <data dir>/phase5/nxf; work files go in
# NXF_DIR/work).
set -euo pipefail
: "${REPO_DIR:?}" "${CFDNA_DATA:?}"
source /app/lmod/lmod/init/bash
module load Nextflow/24.10.2
export NXF_SLURM_QUEUE=${NXF_SLURM_QUEUE:-campus-new}
NXF_DIR=${NXF_DIR:-$CFDNA_DATA/phase5/nxf}
mkdir -p "$NXF_DIR"
cd "$NXF_DIR"
export NXF_HOME=$PWD/.nxf_home NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$REPO_DIR/pipeline/main.nf" -profile slurm,conda -resume "$@"
echo NXF_DONE
