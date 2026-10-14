#!/usr/bin/env bash
# FASTQ-deletion test on the bundled 1 Mb chr21 slice.
# Makes paired FASTQ from the bundled BAM, runs the pipeline from FASTQ, and checks that
# (1) the raw and trimmed FASTQ are deleted and the published BAM is intact,
# (1b) by default local --input_fastq files are kept, and (2) --keep_fastq keeps all FASTQ.
# Usage (from the repo root): pipeline/tests/fastq_cleanup_test.sh [profile, default docker]
# SAMTOOLS may name a samtools command (e.g. a `docker run ...` prefix); default `samtools`.
set -euo pipefail
profile=${1:-docker}
samtools=${SAMTOOLS:-samtools}
out=$PWD/results_test_fastq
data=pipeline/tests/data

rm -rf "$out" && mkdir -p "$out/fastq" "$out/ref"
cp "$data/ref/genome.fa" "$data/ref/genome.fa.fai" "$out/ref/"
$samtools collate -u -O "$data/HCT116_test.md.bam" \
  | $samtools fastq -n -1 "$out/fastq/HCT116_test_1.fastq.gz" -2 "$out/fastq/HCT116_test_2.fastq.gz" \
      -0 /dev/null -s /dev/null -
cp "$out"/fastq/*.fastq.gz "$out/ref/"   # untouched copies for the --keep_fastq run

fail() { echo "FAIL: $*" >&2; exit 1; }
check_bam() {
  $samtools quickcheck "$1" || fail "quickcheck $1"
  local m; m=$($samtools idxstats "$1" | awk '{s+=$3} END{print s+0}')
  [ "$m" -gt 0 ] || fail "no mapped reads in $1"
  echo "ok: $1 passes quickcheck, $m mapped reads"
}

# (1) --delete_input_fastq: raw and trimmed FASTQ deleted after BAM_CHECK (the default for
#     downloaded FASTQ). Checks use `find -type f`: Nextflow leaves dangling staging symlinks.
nextflow run pipeline/main.nf -profile "test_fastq,$profile" -work-dir "$out/work" --delete_input_fastq
check_bam "$out/bam/HCT116_test.md.bam"
[ -z "$(find "$out/fastq" -type f -name '*.fastq.gz')" ] || fail "raw FASTQ left in $out/fastq"
[ -z "$(find "$out/work" -type f -name '*.trim.fq.gz')" ] || fail "trimmed FASTQ left in $out/work"
[ "$(awk -F'\t' '$5=="deleted"' "$out/fastq_cleanup.tsv" | wc -l)" -eq 4 ] || fail "fastq_cleanup.tsv"
test -s "$out/qc_table.tsv" || fail "qc_table.tsv"
echo "ok: raw and trimmed FASTQ deleted (4 files)"

# (1b) Default with --input_fastq: the local input FASTQ are kept, only trimmed FASTQ deleted.
mkdir -p "$out/local/fastq" && cp "$out"/ref/*.fastq.gz "$out/local/fastq/"
nextflow run pipeline/main.nf -profile "test_fastq,$profile" -work-dir "$out/local/work" \
  --outdir "$out/local" --input_fastq "$out/local/fastq/*_{1,2}.fastq.gz"
check_bam "$out/local/bam/HCT116_test.md.bam"
[ "$(find "$out/local/fastq" -type f -name '*.fastq.gz' | wc -l)" -eq 2 ] || fail "local input FASTQ deleted without --delete_input_fastq"
[ -z "$(find "$out/local/work" -type f -name '*.trim.fq.gz')" ] || fail "trimmed FASTQ left in $out/local/work"
echo "ok: local input FASTQ kept, trimmed FASTQ deleted"

# (2) --keep_fastq: all FASTQ kept.
mkdir -p "$out/keep/fastq" && mv "$out"/ref/*.fastq.gz "$out/keep/fastq/"
nextflow run pipeline/main.nf -profile "test_fastq,$profile" -work-dir "$out/keep/work" \
  --keep_fastq --outdir "$out/keep" --input_fastq "$out/keep/fastq/*_{1,2}.fastq.gz"
check_bam "$out/keep/bam/HCT116_test.md.bam"
[ "$(find "$out/keep/fastq" -type f -name '*.fastq.gz' | wc -l)" -eq 2 ] || fail "--keep_fastq lost raw FASTQ"
[ "$(find "$out/keep/work" -type f -name '*.trim.fq.gz' | wc -l)" -eq 2 ] || fail "--keep_fastq lost trimmed FASTQ"
[ ! -e "$out/keep/fastq_cleanup.tsv" ] || fail "--keep_fastq wrote fastq_cleanup.tsv"
echo "ok: --keep_fastq keeps raw and trimmed FASTQ"
echo FASTQ_CLEANUP_TEST_PASS
