"""Phase 5 batch control: can off-marker features tell spike fragments from background?

Pre-registered (pipeline/PREREGISTERED.md, "Batch control"): cross-validated logistic
regression, spike-source fragments vs threshold-donor fragments, on
  (a) fragment-level technical features: insert size, mapping quality, mismatches (NM) and
      read length (read 1 of proper pairs in the off-marker blocks, duplicates removed,
      MAPQ >= 10; pipeline FRAG_FEATURES);
  (b) off-marker methylation: per fragment, its methylated share and its number of CpG calls
      in the off-marker blocks (duplicates removed).
An AUC well above 0.5 is reported as a batch confound. Matched depth: each class contributes
the same number of fragments, and every threshold donor contributes the same number.

Interpretation note for (b): the spike source is a cancer cell line, so off-marker
methylation can differ for biological reasons as well as for batch reasons.

Writes results/phase5/batch_control.tsv.

Usage: uv run python scripts/phase5_batch.py [--per-donor 2000]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import phase5 as p5

SEED = 20261010
TECH = ["tlen", "mapq", "nm", "read_len"]


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def off_meth(cache: Path) -> pd.DataFrame:
    """Per fragment: methylated share and number of calls in the off-marker blocks."""
    z = np.load(cache)
    df = pd.DataFrame({"row": z["off_row"], "m": z["off_m"], "u": z["off_u"],
                       "w": z["off_w"]})
    g = df.groupby("row").agg(m=("m", "sum"), u=("u", "sum"), w=("w", "first"))
    g = g.loc[g.index.repeat(g.w)]
    return pd.DataFrame({"meth_share": g.m / (g.m + g.u), "n_calls": g.m + g.u})


def cv_auc(X: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """5-fold CV AUC of standardized logistic regression; also the fold-AUC SD."""
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    s = cross_val_predict(model, X, y, cv=cv, method="decision_function")
    folds = [p5.auc(s[te][y[te] == 1], s[te][y[te] == 0]) for _, te in cv.split(X, y)]
    return p5.auc(s[y == 1], s[y == 0]), float(np.std(folds, ddof=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-donor", type=int, default=2000)
    ap.add_argument("--samples", type=Path, default=Path("pipeline/samples/phase5_samples.tsv"))
    ap.add_argument("--spike", default="SRR28401050")
    args = ap.parse_args()
    data = data_dir() / "phase5"
    rng = np.random.default_rng(SEED)
    smp = pd.read_csv(args.samples, sep="\t")
    donors = smp.run[smp.role == "threshold"].tolist()
    n_bg = args.per_donor * len(donors)

    rows = []
    for name, load, cols in (
        ("technical", lambda r: pd.read_csv(data / "frag_features" / f"{r}.frag_features.tsv.gz",
                                            sep="\t", na_values="NA").dropna(), TECH),
        ("off_marker_methylation",
         lambda r: off_meth(data / "records" / f"{r}.dedup.npz"), ["meth_share", "n_calls"]),
    ):
        bg = []
        for d in donors:
            x = load(d)
            if len(x) < args.per_donor:
                raise SystemExit(f"{d}: {len(x)} fragments < {args.per_donor}")
            bg.append(x.iloc[rng.choice(len(x), args.per_donor, replace=False)])
        sp = load(args.spike)
        sp = sp.iloc[rng.choice(len(sp), n_bg, replace=False)]
        X = np.concatenate([pd.concat(bg)[cols].to_numpy(float), sp[cols].to_numpy(float)])
        y = np.r_[np.zeros(n_bg, int), np.ones(n_bg, int)]
        a, sd = cv_auc(X, y)
        means = {f"{c}_bg": float(pd.concat(bg)[c].mean()) for c in cols}
        means |= {f"{c}_spike": float(sp[c].mean()) for c in cols}
        rows.append({"features": name, "n_per_class": n_bg, "cv_auc": a, "fold_auc_sd": sd,
                     **means})
        print(name, round(a, 4), means, flush=True)
    pd.DataFrame(rows).to_csv(Path("results/phase5/batch_control.tsv"), sep="\t", index=False,
                              float_format="%.4g")


if __name__ == "__main__":
    main()
