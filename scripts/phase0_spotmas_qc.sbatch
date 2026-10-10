#!/usr/bin/env bash
#SBATCH --job-name=spotmas-qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.log
# Phase 0 QC of one SPOT-MAS run: bisulfite conversion (non-CpG methylation), M-bias
# (Adaptase tails) and coverage. Uses Bismark against a prebuilt hg38 bisulfite index.
# Usage: RUNS="SRR1 SRR2" BISMARK_INDEX=<dir> CFDNA_DATA=<dir> sbatch --array=0-1 phase0_spotmas_qc.sbatch
# CLIP_R2_5 / CLIP_R1_3 / CLIP_R2_3 (default 0) clip fixed lengths from the R2 5' end, the R1 3'
# end and the R2 3' end. The xGen Adaptase tail sits at the R2 5' end. CLIP_R2_3 is passed
# explicitly because fastp copies --trim_tail1 to read 2 when --trim_tail2 is unset.
# With 0, the M-bias table shows the raw end effects.
set -euo pipefail
read -ra runs <<<"$RUNS"
run=${runs[$SLURM_ARRAY_TASK_ID]}
: "${BISMARK_INDEX:?}" "${CFDNA_DATA:?}"
in=${IN_DIR:-$CFDNA_DATA/sra/fastq}
out=${OUT_DIR:-$CFDNA_DATA/qc/spotmas}/$run; mkdir -p "$out"; cd "$out"
source /etc/profile.d/lmod.sh 2>/dev/null || true

module purge; module load fastp/0.23.4-GCC-13.2.0
fastp -i "$in/${run}_1.fastq.gz" -I "$in/${run}_2.fastq.gz" \
  -o "${run}_1.trim.fq.gz" -O "${run}_2.trim.fq.gz" \
  --detect_adapter_for_pe --disable_quality_filtering --length_required 30 \
  --trim_front2 "${CLIP_R2_5:-0}" --trim_tail1 "${CLIP_R1_3:-0}" \
  --trim_tail2 "${CLIP_R2_3:-0}" \
  -w 8 -j "${run}.fastp.json" -h "${run}.fastp.html"

module purge; module load Bismark/0.24.1-GCC-12.2.0
bismark --genome "$BISMARK_INDEX" --parallel 4 -p 2 --temp_dir "$out" \
  -1 "${run}_1.trim.fq.gz" -2 "${run}_2.trim.fq.gz" -o "$out"
bam=$(ls ./*_bismark_bt2_pe.bam)
deduplicate_bismark -p --bam "$bam" --output_dir "$out"
dedup=${bam%.bam}.deduplicated.bam
bismark_methylation_extractor -p --no_overlap --comprehensive --gzip --bedGraph \
  --parallel 4 --buffer_size 16G -o "$out" "$dedup"
samtools stats -@ 8 "$dedup" > "${run}.dedup.samtools_stats.txt"
rm -f "${run}"_?.trim.fq.gz
echo done
