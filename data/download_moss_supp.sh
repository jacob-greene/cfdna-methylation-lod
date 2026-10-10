#!/usr/bin/env bash
# Moss et al. 2018 (Nat Commun 9:5068) Supplementary Data 1 (sample and mixing details).
source "$(dirname "$0")/common.sh"
fetch "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-07466-6/MediaObjects/41467_2018_7466_MOESM3_ESM.xlsx" \
  "$CFDNA_DATA/moss_supp/41467_2018_7466_MOESM3_ESM.xlsx"
