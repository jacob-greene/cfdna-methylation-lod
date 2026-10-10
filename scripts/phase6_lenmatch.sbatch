#!/usr/bin/env bash
#SBATCH --job-name=p6-lenmatch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
# Phase 6 length-matched spike (scripts/phase6_lenmatch.py): HCT116 BAM -> read names and
# insert sizes in the marker and off-marker blocks -> selection -> filtered BAM -> bam2pat
# (duplicates removed, as in pipeline/main.nf BAM2PAT) -> records cache.
# Usage: sbatch -p campus-new --export=ALL,REPO_DIR=<repo>,CFDNA_DATA=<data dir> \
#          -o <log> scripts/phase6_lenmatch.sbatch
set -euo pipefail
: "${REPO_DIR:?}" "${CFDNA_DATA:?}"
RUN=SRR28401050
P5=$CFDNA_DATA/phase5
OUT=$P5/lenmatch
WGBS=$CFDNA_DATA/tools/wgbs_tools-0.3.0
export PATH=$CFDNA_DATA/tools/p5env/bin:$PATH
mkdir -p "$OUT" && cd "$OUT"
cat "$P5/blocks_cpg.tsv" "$P5/offmarker_cpg.tsv" | cut -f1-3 | sort -k1,1 -k2,2n > regions.bed
samtools view -F 3844 -L regions.bed "$P5/bam/$RUN.md.bam" | cut -f1,9 | gzip > qname_tlen.tsv.gz
cd "$REPO_DIR"
uv run --offline python scripts/phase6_lenmatch.py select "$OUT/qname_tlen.tsv.gz" \
    "$OUT/keep.txt" "$OUT/tlen_qc.tsv"
cd "$OUT"
samtools view -@ 4 -b -N keep.txt -o "${RUN}_len.bam" "$P5/bam/$RUN.md.bam"
samtools index "${RUN}_len.bam"
rm -rf wt dd && mkdir -p wt/references dd && cp -r "$WGBS/src" wt/
ln -s "$(readlink -f "$WGBS/references/hg38")" wt/references/hg38
ln -sfn hg38 wt/references/default
"$WGBS/.venv/bin/python" wt/src/python/wgbs_tools.py bam2pat --genome hg38 -L regions.bed \
    -@ 4 -F 1796 -o dd "${RUN}_len.bam"
mv "dd/${RUN}_len.pat.gz" "$P5/pat/${RUN}_len.dedup.pat.gz"
samtools view -f 67 -F 1796 -q 10 -L "$P5/offmarker_cpg.tsv" "${RUN}_len.bam" \
    | awk '{t=$9<0?-$9:$9; print t}' | gzip > offmarker_tlen.txt.gz
cd "$REPO_DIR"
uv run --offline python scripts/phase6_lenmatch.py records \
    "$P5/pat/${RUN}_len.dedup.pat.gz" "$P5/records/${RUN}_len.dedup.npz"
echo LENMATCH_DONE
