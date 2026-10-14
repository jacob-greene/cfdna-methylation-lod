"""Phase 5: thresholds, detection rates, LoD, held-out FPR, AUC and the external check.

Reads results/phase5/scores.tsv.gz (scripts/phase5_score.py). Rules (pipeline/PREREGISTERED.md):
  - untrained scores (LLR_binom, LLR_bb, UXM): threshold = one-sided 95 % normal prediction
    bound over the threshold pools of the same k (src/phase5.prediction_bound_pools);
  - transfer models (T_*): the Phase 3b valid-rule threshold stored with the model;
  - LLR at a spiked level f uses the alternative f_alt = f; f = 0 rows give the FPR for each
    f_alt;
  - LoD: lowest f detected in >= 95 % of eval tests, interpolated on log f
    (src/simulate.lod_interp).

Writes results/phase5/{thresholds,detection,lod,fpr,external}.tsv.

Usage: uv run python scripts/phase5_summary.py
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import phase5 as p5
from src import simulate as sim

F_GRID = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
KEY = ["rung", "score", "set", "f_alt", "k"]


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def main() -> None:
    out = Path("results/phase5")
    d = pd.read_csv(out / "scores.tsv.gz", sep="\t")
    d["f_alt"] = d.f_alt.fillna(-1.0)  # -1 = score without an alternative level

    thr = pd.DataFrame(
        [(*key, p5.prediction_bound_pools(g.value.to_numpy(), int(key[-1])), len(g))
         # Transfer models keep their stored Phase 3b threshold (below), not a pool bound.
         for key, g in d[(d.kind == "threshold") & ~d.score.str.startswith("T_")].groupby(KEY)],
        columns=[*KEY, "threshold", "n_pools"])
    tdir = data_dir() / "phase5" / "transfer"
    t_rows = []
    for pk in sorted(tdir.glob("*.pkl")):
        with open(pk, "rb") as fh:
            m = pickle.load(fh)
        for k in sorted(d.k.unique()):
            t_rows.append(("real", f"T_{m['model']}", m["features"], -1.0, k, m["threshold"],
                           0))
    thr = pd.concat([thr, pd.DataFrame(t_rows, columns=[*KEY, "threshold", "n_pools"])])
    thr.to_csv(out / "thresholds.tsv", sep="\t", index=False, float_format="%.6g")

    assert not thr.duplicated(KEY).any(), "a score has more than one threshold"
    ev = d[d.kind != "threshold"].merge(thr, on=KEY, how="left")
    assert ev.threshold.notna().all(), "a score has no threshold"
    ev["detected"] = ev.value > ev.threshold
    e = ev[ev.kind == "eval"]

    # FPR: f = 0 eval tests, per f_alt.
    fpr = (e[e.f == 0].groupby(KEY).agg(fpr=("detected", "mean"), n=("detected", "size"))
           .reset_index())
    fpr.to_csv(out / "fpr.tsv", sep="\t", index=False, float_format="%.4g")

    # Detection at spiked levels: LLR rows only at f_alt = f.
    sp = e[(e.f > 0) & ((e.f_alt == -1) | np.isclose(e.f_alt, e.f))]
    det = (sp.groupby(["rung", "score", "set", "k", "f"])
           .agg(rate=("detected", "mean"), n=("detected", "size"),
                median_score=("value", "median"))
           .reset_index())
    # AUC: spiked tests vs f = 0 tests of the same k (LLR scored at the same f_alt).
    aucs = []
    for (rung, score, s, k, f), g in sp.groupby(["rung", "score", "set", "k", "f"]):
        fa = f if score.startswith("LLR") else -1.0
        neg = e[(e.f == 0) & (e.rung == rung) & (e.score == score) & (e.set == s)
                & (e.k == k) & np.isclose(e.f_alt, fa)].value
        aucs.append(p5.auc(g.value, neg))
    det["auc"] = aucs
    det.to_csv(out / "detection.tsv", sep="\t", index=False, float_format="%.4g")

    lod_rows = []
    for (rung, score, s, k), g in det.groupby(["rung", "score", "set", "k"]):
        g = g.sort_values("f")
        f0 = fpr[(fpr.rung == rung) & (fpr.score == score) & (fpr.set == s) & (fpr.k == k)]
        lod_rows.append((rung, score, s, k, sim.lod_interp(g.f, g.rate),
                         float(f0.fpr.max()), int(f0.n.max()), int(g.n.min())))
    lod = pd.DataFrame(lod_rows, columns=["rung", "score", "set", "k", "lod", "fpr_max",
                                          "n_null", "n_per_f"])
    lod.to_csv(out / "lod.tsv", sep="\t", index=False, float_format="%.4g")

    # External check: 10 cancers vs the 16 eval donors, k = 1.
    c = ev[ev.kind == "cancer"]
    neg = e[(e.f == 0) & (e.k == 1) & (e.rung == "real")]
    ext = []
    for (score, s, fa), g in c.groupby(["score", "set", "f_alt"]):
        n = neg[(neg.score == score) & (neg.set == s) & np.isclose(neg.f_alt, fa)]
        ext.append((score, s, fa, len(g), int(g.detected.sum()), len(n), int(n.detected.sum()),
                    p5.auc(g.value, n.value)))
    pd.DataFrame(ext, columns=["score", "set", "f_alt", "n_cancer", "cancer_detected",
                               "n_healthy", "healthy_detected", "auc"]).to_csv(
        out / "external.tsv", sep="\t", index=False, float_format="%.4g")
    print(lod.to_string(index=False))


if __name__ == "__main__":
    main()
