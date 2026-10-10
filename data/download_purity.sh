#!/usr/bin/env bash
# Aran, Sirota & Butte 2015 (Nat Commun 6:8971), Supplementary Data 1:
# tumour purity estimates (ESTIMATE, ABSOLUTE, LUMP, IHC, CPE consensus) for TCGA samples.
source "$(dirname "$0")/common.sh"
fetch "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx" \
  "$CFDNA_DATA/purity/aran2015_supplementary_data_1.xlsx"
