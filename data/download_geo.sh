#!/usr/bin/env bash
# GSE122126 (Moss et al. 2018): series matrices (processed betas) and the
# Illumina manifests GEO distributes in the supplementary folder.
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/geo
base=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122126
for f in GSE122126-GPL13534_series_matrix.txt.gz GSE122126-GPL21145_series_matrix.txt.gz; do
  fetch "$base/matrix/$f" "$d/$f"
done
fetch "$base/suppl/filelist.txt" "$d/filelist.txt"
# Optional raw IDATs (1.7 GB), only for the process_array.R robustness check.
if [[ ${WITH_IDAT:-0} == 1 ]]; then fetch "$base/suppl/GSE122126_RAW.tar" "$d/GSE122126_RAW.tar"; fi
