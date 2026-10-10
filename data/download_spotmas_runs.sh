#!/usr/bin/env bash
# Download SPOT-MAS FASTQ pairs from ENA with MD5 checks. Usage: download_spotmas_runs.sh SRR... 
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/sra/fastq; mkdir -p "$d"
report=$CFDNA_DATA/sra/PRJNA929650_ena_filereport.tsv
for run in "$@"; do
  line=$(awk -F'\t' -v r="$run" '$1==r' "$report")
  [[ -n $line ]] || { echo "unknown run $run" >&2; exit 1; }
  IFS=';' read -ra urls <<<"$(cut -f7 <<<"$line")"
  IFS=';' read -ra md5s <<<"$(cut -f9 <<<"$line")"
  for i in "${!urls[@]}"; do
    dest=$d/$(basename "${urls[$i]}")
    fetch "https://${urls[$i]}" "$dest"
    echo "${md5s[$i]}  $dest" | md5sum -c --quiet - || { echo "MD5 mismatch: $dest" >&2; exit 1; }
  done
done
