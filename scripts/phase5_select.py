"""Phase 5, step 0: pre-register the SPOT-MAS runs and the healthy-donor split.

Runs before any SPOT-MAS sample is scored. Writes pipeline/samples/phase5_samples.tsv (one
row per run: accession, label, stage, role, pool memberships, ENA bytes) and prints the
download size.

Selection (seed SEED, numpy default_rng):
  - healthy: 30 healthy controls (all SPOT-MAS donors are female). SRR23280038 is kept
    because it is already downloaded (Phase 0 QC only, never scored); 29 are drawn at random
    from the other 277.
  - cancer: 10 early-stage breast cancers (stage I, IA, IB, II, IIA or IIB). SRR23280187 is
    kept (Phase 0 QC only); 9 are drawn at random from the other early-stage runs.
Split of the 30 healthy donors (one random permutation):
  - role "threshold" (14 donors): set every threshold and the background profile p0;
  - role "eval" (16 donors): backgrounds of the spike-in evaluation and the healthy
    comparators of the external cancer check. Never used to set a threshold.
Evaluation pools: the eval donors in permuted order are cut into disjoint pools of k donors,
k = 1, 2, 4, 8 (16, 8, 4 and 2 pools). Column pool_k<k> holds the pool index.

Usage: uv run python scripts/phase5_select.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import spotmas

SEED = 20261010
N_HEALTHY, N_CANCER, N_THRESHOLD = 30, 10, 14
KEEP_HEALTHY, KEEP_CANCER = "SRR23280038", "SRR23280187"
EARLY = {"I", "IA", "IB", "II", "IIA", "IIB"}
POOL_K = (1, 2, 4, 8)
OUT = Path("pipeline/samples/phase5_samples.tsv")


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def select(runs: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    runs = runs.sort_values("Run").reset_index(drop=True)
    healthy = runs[runs.Label == "Healthy control"]
    cancer = runs[(runs.Label == "Breast cancer") & runs["Tumor Stage"].isin(EARLY)]
    h_other = healthy[healthy.Run != KEEP_HEALTHY].Run.to_numpy()
    c_other = cancer[cancer.Run != KEEP_CANCER].Run.to_numpy()
    h = [KEEP_HEALTHY] + sorted(rng.choice(h_other, N_HEALTHY - 1, replace=False))
    c = [KEEP_CANCER] + sorted(rng.choice(c_other, N_CANCER - 1, replace=False))
    order = rng.permutation(np.asarray(h, dtype=object))
    thr, ev = list(order[:N_THRESHOLD]), list(order[N_THRESHOLD:])
    rows = [{"run": r, "role": "threshold", "eval_rank": -1} for r in thr]
    rows += [{"run": r, "role": "eval", "eval_rank": i} for i, r in enumerate(ev)]
    rows += [{"run": r, "role": "cancer", "eval_rank": -1} for r in c]
    out = pd.DataFrame(rows)
    for k in POOL_K:
        out[f"pool_k{k}"] = np.where(out.eval_rank >= 0, out.eval_rank // k, -1)
    meta = runs.set_index("Run")
    out["label"] = out.run.map(meta.Label)
    out["stage"] = out.run.map(meta["Tumor Stage"])
    out["sample_name"] = out.run.map(meta.sample_name)
    return out


def main() -> None:
    data = data_dir()
    runs = spotmas.read_runs(data)
    sel = select(runs, np.random.default_rng(SEED))
    ena = pd.read_csv(data / "sra" / "PRJNA929650_ena_filereport.tsv", sep="\t")
    ena = ena.set_index("run_accession")
    sel["fastq_bytes"] = sel.run.map(
        lambda r: sum(int(b) for b in str(ena.loc[r, "fastq_bytes"]).split(";")))
    sel["read_pairs"] = sel.run.map(ena.read_count)
    assert sel.run.is_unique and len(sel) == N_HEALTHY + N_CANCER
    assert not set(sel[sel.role == "threshold"].run) & set(sel[sel.role == "eval"].run)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    sel.drop(columns="eval_rank").to_csv(OUT, sep="\t", index=False)
    print(sel.groupby("role").agg(n=("run", "size"), gb=("fastq_bytes", lambda b: b.sum() / 1e9),
                                  pairs_median=("read_pairs", "median")))
    print(f"total {sel.fastq_bytes.sum() / 1e9:.2f} GB; already downloaded: "
          f"{KEEP_HEALTHY}, {KEEP_CANCER}")


if __name__ == "__main__":
    main()
