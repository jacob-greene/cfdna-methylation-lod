"""Phase 5, step 0: pre-register the marker CpGs and their 500 bp blocks.

Runs before any SPOT-MAS sample is scored. Inputs are array data only (Phase 3b inputs in
$CFDNA_DATA/phase3b/ and the hg38 HM450 manifest); no read data is read.

Marker sets:
  - "atlas": the 5,258 Moss atlas CpGs of Phase 3b feature set (a).
  - "tumour": the top 500 CpGs hypermethylated in TCGA-COAD tumours vs TCGA normal colon and
    unmethylated in Moss blood cell types (max blood beta <= 0.2), selected exactly as in
    Phase 3b (src/classify.select_tumour_markers), here on all 276 kept tumours. The spike
    source (HCT116) and SPOT-MAS play no part in the selection.
Blocks: a 500 bp window centred on each marker CpG ([pos - 250, pos + 250) on hg38); windows
of the same set that overlap are merged into one block. At ~0.5x coverage most single CpGs
have no read, so reads are counted per block (all CpGs inside the block).

Writes:
  pipeline/markers/phase5_markers.tsv   probe, set, chrom, pos (0-based C of the + strand),
                                        block, tumour_mean, normal_mean, blood_max
  pipeline/markers/phase5_blocks.bed    chrom, start, end, block, set, n_probes (sorted)

Usage: uv run python scripts/phase5_markers.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import classify as cl

N_TUMOUR = 500
HALF = 250
OUT = Path("pipeline/markers")


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def blocks(pos: pd.DataFrame, half: int = HALF) -> pd.Series:
    """Block id per probe: centred windows, merged where windows of one chromosome overlap."""
    pos = pos.sort_values(["chrom", "pos"])
    ids, cur, last_end, last_chr = [], -1, -1, None
    for c, p in zip(pos.chrom, pos.pos):
        if c != last_chr or p - half >= last_end:
            cur += 1
        ids.append(cur)
        last_chr, last_end = c, p + half
    return pd.Series(ids, index=pos.index)


def main() -> None:
    data = data_dir()
    p3 = data / "phase3b"
    T = pd.read_parquet(p3 / "tumours.parquet")
    N = pd.read_parquet(p3 / "normals.parquet")
    B = pd.read_parquet(p3 / "blood.parquet")
    atlas_cpgs = (p3 / "atlas_cpgs.txt").read_text().split()
    idx = cl.select_tumour_markers(T.to_numpy(float), N.to_numpy(float), B.to_numpy(float),
                                   N_TUMOUR)
    tumour_cpgs = list(T.index[idx])

    man = atlas_mod.read_manifest(data, "450K")
    rows = []
    for name, cpgs in (("atlas", atlas_cpgs), ("tumour", tumour_cpgs)):
        m = man.reindex(cpgs)
        miss = m.CpG_chrm.isna().sum()
        print(f"{name}: {len(cpgs)} CpGs, {miss} without hg38 coordinates (dropped)")
        m = m.dropna(subset=["CpG_chrm"])
        rows.append(pd.DataFrame({"probe": m.index, "set": name, "chrom": m.CpG_chrm.values,
                                  "pos": m.CpG_beg.astype(int).values}))
    mk = pd.concat(rows, ignore_index=True)
    mk["tumour_mean"] = T.reindex(mk.probe).mean(axis=1).values
    mk["normal_mean"] = N.reindex(mk.probe).mean(axis=1).values
    mk["blood_max"] = B.reindex(mk.probe).max(axis=1).values
    mk["block"] = ""
    for name, g in mk.groupby("set"):
        b = blocks(g)
        mk.loc[b.index, "block"] = [f"{name}_{i:05d}" for i in b]
    mk = mk.sort_values(["set", "block", "pos"]).reset_index(drop=True)

    bed = mk.groupby("block").agg(chrom=("chrom", "first"), start=("pos", "min"),
                                  end=("pos", "max"), set=("set", "first"),
                                  n_probes=("probe", "size")).reset_index()
    bed["start"] = (bed.start - HALF).clip(lower=0)
    bed["end"] = bed.end + HALF
    bed = bed[["chrom", "start", "end", "block", "set", "n_probes"]]
    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)
    OUT.mkdir(parents=True, exist_ok=True)
    mk.to_csv(OUT / "phase5_markers.tsv", sep="\t", index=False, float_format="%.4f")
    bed.to_csv(OUT / "phase5_blocks.bed", sep="\t", index=False, header=False)
    for name, g in bed.groupby("set"):
        w = (g.end - g.start)
        print(f"{name}: {len(g)} blocks, {int(g.n_probes.sum())} probes, "
              f"width median {int(w.median())} bp, total {w.sum() / 1e6:.2f} Mb")
    t = mk[mk.set == "tumour"]
    print(f"tumour markers: tumour mean median {t.tumour_mean.median():.3f}, normal mean median "
          f"{t.normal_mean.median():.3f}, blood max median {t.blood_max.median():.3f}")
    print("tumour/atlas probe overlap:",
          len(set(t.probe) & set(mk[mk.set == "atlas"].probe)))
    print("chromosomes:", sorted(mk.chrom.unique(), key=lambda c: (len(c), c))[-3:],
          "sex chromosomes:", int(mk.chrom.isin(["chrX", "chrY"]).sum()))
    assert np.isfinite(mk[["tumour_mean", "blood_max"]].to_numpy()).all()


if __name__ == "__main__":
    main()
