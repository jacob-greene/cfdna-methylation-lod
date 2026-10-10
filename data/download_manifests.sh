#!/usr/bin/env bash
# Probe annotation on hg38 (Zhou et al. 2017, InfiniumAnnotation), pinned commit.
# Illumina's own manifests (hg19) ship inside the GEO supplementary folder.
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/manifests
commit=6d176d6e7e56465403d4478a09fe11b983eeb721
raw=https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationData/$commit/Anno
fetch "$raw/EPIC/EPIC.hg38.manifest.tsv.gz" "$d/EPIC.hg38.manifest.tsv.gz"
fetch "$raw/HM450/HM450.hg38.manifest.tsv.gz" "$d/HM450.hg38.manifest.tsv.gz"
# Illumina EPIC v1.0 manifest (hg19 coordinates) as distributed with GPL21145.
fetch https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL21nnn/GPL21145/suppl/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz "$d/GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
