"""Phase 0 item 2b: atlas CpGs covered by >= 1 read in a SPOT-MAS sample.

Joins a Bismark `.bismark.cov.gz` (1-based C positions, strands not merged) to the hg38
coordinates of the deduplicated atlas CpGs. The hg38 manifest gives CpG_beg as the 0-based start
of the CpG dinucleotide, so the + strand C is at CpG_beg + 1 and the - strand C at CpG_beg + 2
(both 1-based). A CpG counts as covered when either strand has >= 1 call.

Positive control: the same join on all manifest CpGs gives the genome-wide covered fraction,
which must be close to 1 - exp(-depth) for Poisson coverage.

Usage: uv run python scripts/phase0_atlas_coverage.py <run>.bismark.cov.gz [...]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import deconv


def main() -> None:
    data = Path(os.environ.get("CFDNA_DATA",
                               Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"))
    atl = deconv.dedup_atlas(atlas_mod.read_atlas(data))
    man = atlas_mod.read_manifest(data, "450K").dropna(subset=["CpG_chrm", "CpG_beg"])
    man = man[man.index.str.startswith("cg")]
    man["CpG_beg"] = man["CpG_beg"].astype(int)
    man["is_atlas"] = man.index.isin(atl.index)
    print(f"atlas CpGs with hg38 coordinates: {man.is_atlas.sum()} of {len(atl)}")
    rows = []
    for path in sys.argv[1:]:
        cov = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 4, 5],
                          names=["chr", "pos", "meth", "unmeth"])
        cov["n"] = cov.meth + cov.unmeth
        keys = cov.set_index(["chr", "pos"])["n"]
        plus = pd.MultiIndex.from_arrays([man.CpG_chrm, man.CpG_beg + 1])
        minus = pd.MultiIndex.from_arrays([man.CpG_chrm, man.CpG_beg + 2])
        n = keys.reindex(plus).fillna(0).to_numpy() + keys.reindex(minus).fillna(0).to_numpy()
        man["calls"] = n
        for label, sub in [("atlas", man[man.is_atlas]), ("all 450K cg probes", man)]:
            rows.append({"file": Path(path).name.split("_")[0], "set": label, "cpgs": len(sub),
                         "covered_ge1": int((sub.calls >= 1).sum()),
                         "frac_covered": float((sub.calls >= 1).mean()),
                         "mean_calls": float(sub.calls.mean())})
    out = pd.DataFrame(rows)
    print(out.to_string(index=False, float_format="%.4f"))


if __name__ == "__main__":
    main()
