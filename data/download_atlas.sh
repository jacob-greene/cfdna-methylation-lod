#!/usr/bin/env bash
# Moss et al. 2018 reference atlas (github.com/nloyfer/meth_atlas), pinned commit.
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/meth_atlas
commit=${ATLAS_COMMIT:?set ATLAS_COMMIT}
raw=https://raw.githubusercontent.com/nloyfer/meth_atlas/$commit
for f in reference_atlas.csv full_atlas.csv.gz deconvolve.py README.md LICENSE.md examples.csv pre_process/process_array.R pre_process/ref_sample.RData; do
  fetch "$raw/$f" "$d/$f"
done
