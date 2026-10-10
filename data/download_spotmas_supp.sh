#!/usr/bin/env bash
# Nguyen et al. 2023 (Front Oncol, doi:10.3389/fonc.2023.1127086, PMC10200909) supplementary
# files via Europe PMC. Table_1.docx = Table S1, per-sample label, stage and subtype.
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/sra/spotmas_supp
fetch "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10200909/supplementaryFiles?includeInlineImage=false" "$d/supp.zip"
(cd "$d" && unzip -o -q supp.zip)
echo "7d5f5ad3bc0a7ec22a2b44517d0cd82d  $d/Table_1.docx" | md5sum -c -
