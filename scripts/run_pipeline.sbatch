#!/usr/bin/env bash
#SBATCH --job-name=cfdna-lod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
# The Nextflow head process for the whole pipeline (reads, analysis, figures). It submits every
# pipeline task as its own Slurm job.
# Usage: sbatch -p campus-new --export=ALL,REPO_DIR=<repo>,CFDNA_DATA=<data dir> \
#          -o <log> scripts/run_pipeline.sbatch [extra nextflow args]
# NXF_DIR sets the Nextflow launch directory (default <data dir>/nxf; work files go in
# NXF_DIR/work).
set -euo pipefail
: "${REPO_DIR:?}" "${CFDNA_DATA:?}"
source /app/lmod/lmod/init/bash
module load Nextflow/24.10.2
export NXF_SLURM_QUEUE=${NXF_SLURM_QUEUE:-campus-new}
NXF_DIR=${NXF_DIR:-$CFDNA_DATA/nxf}
mkdir -p "$NXF_DIR"
cd "$NXF_DIR"
export NXF_HOME=$PWD/.nxf_home NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$REPO_DIR/pipeline/main.nf" -profile slurm,conda -resume "$@"
echo NXF_DONE
