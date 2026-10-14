"""Phase 5: off-marker blocks for the batch control (pipeline/markers/phase5_offmarker.bed).

Draws 2,000 random 500 bp autosomal windows that lie at least 10 kb from every marker block,
contain no N and hold at least 5 CpGs. The batch control reads technical features and
methylation of fragments in these windows only, so it never sees a marker CpG.
Usage: CFDNA_DATA=<dir> uv run python scripts/phase5_offmarker.py
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261010
N_WINDOWS = 2000
WIDTH = 500
MIN_GAP = 10_000
MIN_CPG = 5

root = Path(__file__).resolve().parent.parent
data = Path(os.environ.get("CFDNA_DATA", root.parent / "cfdna-methylation-lod-data"))
fasta = data / "ref/hg38/genome.fa"
fai = pd.read_csv(f"{fasta}.fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "length", "offset", "bases", "bytes"]).set_index("chrom")
autosomes = [f"chr{i}" for i in range(1, 23)]
fai = fai.loc[autosomes]
blocks = pd.read_csv(root / "pipeline/markers/phase5_blocks.bed", sep="\t", header=None,
                     usecols=[0, 1, 2], names=["chrom", "start", "end"])


def fetch(fh, chrom, start, end):
    r = fai.loc[chrom]
    seq = []
    for pos in range(start, end):
        line, col = divmod(pos, r["bases"])
        fh.seek(r["offset"] + line * r["bytes"] + col)
        seq.append(fh.read(1))
    return b"".join(seq).decode().upper()


def near_block(chrom, start, end):
    b = blocks[blocks["chrom"] == chrom]
    return bool(((b["start"] < end + MIN_GAP) & (b["end"] > start - MIN_GAP)).any())


rng = np.random.default_rng(SEED)
weights = fai["length"] / fai["length"].sum()
rows = []
with open(fasta, "rb") as fh:
    while len(rows) < N_WINDOWS:
        chrom = rng.choice(autosomes, p=weights.to_numpy())
        start = int(rng.integers(0, fai.loc[chrom, "length"] - WIDTH))
        end = start + WIDTH
        if near_block(chrom, start, end):
            continue
        if any(r[0] == chrom and r[1] < end and r[2] > start for r in rows):
            continue
        seq = fetch(fh, chrom, start, end + 1)
        if "N" in seq or seq.count("CG") < MIN_CPG:
            continue
        rows.append((chrom, start, end, f"off_{len(rows):05d}", "offmarker", seq.count("CG")))

out = pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "set", "n_cpg"])
out["chrom"] = pd.Categorical(out["chrom"], autosomes, ordered=True)
out = out.sort_values(["chrom", "start"])
out.to_csv(root / "pipeline/markers/phase5_offmarker.bed", sep="\t", header=False, index=False)
print(f"{len(out)} windows, {out['end'].sub(out['start']).sum() / 1e6:.2f} Mb, "
      f"median CpGs {out['n_cpg'].median():.0f}")
