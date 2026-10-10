#!/usr/bin/env bash
# TCGA open-access 450K "Methylation Beta Value" files (GDC). PROJECT defaults to TCGA-COAD.
source "$(dirname "$0")/common.sh"
project=${PROJECT:-TCGA-COAD}
d=$CFDNA_DATA/tcga/$project; mkdir -p "$d/files"
filters=$(cat <<JSON
{"op":"and","content":[
 {"op":"in","content":{"field":"cases.project.project_id","value":["$project"]}},
 {"op":"in","content":{"field":"data_type","value":["Methylation Beta Value"]}},
 {"op":"in","content":{"field":"platform","value":["Illumina Human Methylation 450"]}},
 {"op":"in","content":{"field":"access","value":["open"]}}]}
JSON
)
fields=file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.tissue_type
curl -fsS -X POST https://api.gdc.cancer.gov/files -H 'Content-Type: application/json' \
  -d "{\"filters\":$filters,\"fields\":\"$fields\",\"format\":\"TSV\",\"size\":\"5000\"}" > "$d/manifest_raw.tsv"
python3 - "$d" <<'PY'
import csv, sys
d = sys.argv[1]
rows = list(csv.DictReader(open(f"{d}/manifest_raw.tsv"), delimiter="\t"))
with open(f"{d}/manifest.tsv", "w") as out:
    out.write("id\tfilename\tmd5\tsize\tcase\tsample\tsample_type\n")
    for r in rows:
        out.write("\t".join([r["file_id"], r["file_name"], r["md5sum"], r["file_size"],
                             r["cases.0.submitter_id"], r["cases.0.samples.0.submitter_id"],
                             r["cases.0.samples.0.sample_type"]]) + "\n")
print(len(rows), "files")
PY
if [[ ${MANIFEST_ONLY:-0} == 1 ]]; then exit 0; fi
# Four parallel fetches; each file is MD5-checked against the GDC manifest.
export -f fetch; export d
tail -n +2 "$d/manifest.tsv" | cut -f1-3 | xargs -P 4 -L 1 bash -c '
  id=$0 fn=$1 md5=$2 dest=$d/files/$1
  if [[ ! -s $dest ]]; then fetch "https://api.gdc.cancer.gov/data/$id" "$dest"; fi
  echo "$md5  $dest" | md5sum -c --quiet - || { echo "MD5 mismatch: $fn (deleted; re-run to re-fetch)" >&2; rm -f "$dest"; exit 255; }'
