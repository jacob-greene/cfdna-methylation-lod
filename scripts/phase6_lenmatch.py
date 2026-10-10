"""Phase 6 length-matched spike (pipeline/PREREGISTERED_PHASE6.md, amendment 1).

The HCT116 spike library has longer inserts than SPOT-MAS cfDNA. This rung keeps a random
subset of HCT116 fragments whose insert-size distribution matches the pooled 14 threshold
donors. Reads cannot be trimmed in silico without losing the unsequenced fragment middle,
so fragments are selected, not trimmed.

Rule:
  - Bins: 25 equal-count bins of the threshold donors' insert sizes (frag_features, read 1,
    off-marker, MAPQ >= 10, duplicates excluded), outer edges 0 and 1000 bp.
  - p_s = 1/25 per bin; p_h = share of HCT116 fragments per bin (fragments in marker or
    off-marker blocks, one per read name, 0 < |TLEN| < 1000; others are dropped).
  - Keep probability per bin w = (p_s / p_h) / max(p_s / p_h). Draws in sorted read-name
    order with rng = default_rng([SEED, 70]).

Subcommands:
  select  <qname_tlen.tsv.gz> <keep.txt> <qc.tsv>   choose read names
  records <pat.gz> <out.npz>                          pat file -> records cache (Phase 5 rule)
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location("phase5_score", ROOT / "scripts/phase5_score.py")
p5s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5s)

N_BINS, MAX_TLEN = 25, 1000


def threshold_tlen(samples: Path) -> np.ndarray:
    smp = pd.read_csv(samples, sep="\t")
    thr = smp.run[smp.role == "threshold"].tolist()
    assert len(thr) == 14
    d = p5s.data_dir() / "phase5" / "frag_features"
    t = np.concatenate([pd.read_csv(d / f"{r}.frag_features.tsv.gz", sep="\t",
                                    usecols=["tlen"]).tlen.to_numpy() for r in thr])
    return t[(t > 0) & (t < MAX_TLEN)]


def keep_prob(ref: np.ndarray, tlen: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bin edges and keep probability per bin (module docstring rule)."""
    edges = np.quantile(ref, np.linspace(0, 1, N_BINS + 1))
    edges[0], edges[-1] = 0, MAX_TLEN
    ph = np.histogram(tlen, edges)[0] / len(tlen)
    ratio = np.where(ph > 0, (1 / N_BINS) / np.maximum(ph, 1e-300), 0.0)
    return edges, ratio / ratio.max()


def select(qfile: Path, keep: Path, qc: Path) -> None:
    q = pd.read_csv(qfile, sep="\t", header=None, names=["qname", "tlen"])
    q["tlen"] = q.tlen.abs()
    q = q.drop_duplicates("qname").sort_values("qname", kind="stable")
    q = q[(q.tlen > 0) & (q.tlen < MAX_TLEN)]
    ref = threshold_tlen(ROOT / "pipeline/samples/phase5_samples.tsv")
    edges, w = keep_prob(ref, q.tlen.to_numpy())
    b = np.clip(np.searchsorted(edges, q.tlen.to_numpy(), side="right") - 1, 0, N_BINS - 1)
    kept = np.random.default_rng([p5s.SEED, 70]).random(len(q)) < w[b]
    q.qname[kept].to_csv(keep, index=False, header=False)
    qs = [0.05, 0.25, 0.5, 0.75, 0.95]
    rows = [("threshold_donors", len(ref), *np.quantile(ref, qs)),
            ("hct116_all", len(q), *np.quantile(q.tlen, qs)),
            ("hct116_kept", int(kept.sum()), *np.quantile(q.tlen[kept], qs))]
    pd.DataFrame(rows, columns=["source", "n", *[f"q{int(x * 100)}" for x in qs]]).to_csv(
        qc, sep="\t", index=False)
    print(pd.read_csv(qc, sep="\t").to_string(index=False))


def records(pat: Path, out: Path) -> None:
    data = p5s.data_dir() / "phase5"
    blocks = p5s.load_blocks(data / "blocks_cpg.tsv")
    off = p5s.load_blocks(data / "offmarker_cpg.tsv")
    bsets = {s: blocks[blocks.set == s].reset_index(drop=True) for s in p5s.SETS}
    bsets["off"] = off
    ranges = {s: (b.startCpG.to_numpy(), b.endCpG.to_numpy()) for s, b in bsets.items()}
    print(p5s.build_records((str(pat), str(out), ranges)))


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    {"select": select, "records": records}[cmd](*map(Path, a))
