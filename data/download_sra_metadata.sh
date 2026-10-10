#!/usr/bin/env bash
# PRJNA929650 (SPOT-MAS breast, Nguyen et al. 2023): run-level metadata via
# Entrez E-utilities (runinfo CSV) and BioSample attributes (labels).
source "$(dirname "$0")/common.sh"
d=$CFDNA_DATA/sra; mkdir -p "$d"
eu=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
q="term=PRJNA929650%5BBioProject%5D&usehistory=y"
xml=$(curl -fsS "$eu/esearch.fcgi?db=sra&$q&retmax=0")
we=$(sed -n 's:.*<WebEnv>\(.*\)</WebEnv>.*:\1:p' <<<"$xml"); qk=$(sed -n 's:.*<QueryKey>\(.*\)</QueryKey>.*:\1:p' <<<"$xml")
curl -fsS "$eu/efetch.fcgi?db=sra&WebEnv=$we&query_key=$qk&rettype=runinfo&retmode=text" > "$d/PRJNA929650_runinfo.csv"
# BioSample attributes, fetched by accession in batches of 200.
tail -n +2 "$d/PRJNA929650_runinfo.csv" | cut -d, -f26 | sort -u > "$d/biosamples.txt"
: > "$d/PRJNA929650_biosample.xml"
split -l 200 "$d/biosamples.txt" "$d/bs_chunk_"
for c in "$d"/bs_chunk_*; do
  curl -fsS "$eu/efetch.fcgi?db=biosample&retmode=xml&id=$(paste -sd, "$c")" >> "$d/PRJNA929650_biosample.xml"
  rm "$c"; sleep 1
done
# ENA file report: FASTQ URLs, sizes and MD5s.
curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA929650&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5,instrument_model,library_layout,library_strategy&format=tsv" > "$d/PRJNA929650_ena_filereport.tsv"
wc -l "$d"/PRJNA929650_*
