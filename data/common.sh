# Shared settings for download scripts. Source, do not execute.
# CFDNA_DATA: directory that holds all downloaded data (outside the repo).
set -euo pipefail
: "${CFDNA_DATA:=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/cfdna-methylation-lod-data}"
export CFDNA_DATA
mkdir -p "$CFDNA_DATA"
fetch() {  # fetch <url> <dest>; skip if dest exists and is non-empty
  local url=$1 dest=$2
  mkdir -p "$(dirname "$dest")"
  if [[ -s $dest ]]; then echo "exists: $dest"; return; fi
  curl -fsSL --retry 5 --retry-delay 5 -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}
