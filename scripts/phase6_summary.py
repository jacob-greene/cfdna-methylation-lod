"""Phase 6: thresholds, detection, LoD, FPR with intervals, matched-FPR comparison to Phase 5.

Reads results/phase6/{threshold_scores.tsv, eval_scores.tsv.gz, effective_share.tsv}
(scripts/phase6_score.py) and results/phase5/scores.tsv.gz. Rules
(pipeline/PREREGISTERED_PHASE6.md):
  - "bound" threshold: one-sided 95 % normal prediction bound over the threshold pools of the
    same k (src/phase5.prediction_bound_pools), the Phase 5 rule;
  - "model" threshold (FRAG only): Z > 1.645 (O4);
  - LoD: lowest f detected in >= 95 % of eval tests, interpolated on log f;
  - LoD on effective share: the LoD mapped through the mean effective share per f
    (interpolated on log-log axes);
  - matched FPR (descriptive): for Phase 5 and Phase 6 scores alike, the threshold is the
    highest eval f = 0 score of the same k (eval FPR = 0 by construction).

Writes results/phase6/{thresholds,fpr,detection,lod,matched,external}.tsv.

Usage: uv run python scripts/phase6_summary.py [results_dir]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import phase5 as p5
from src import simulate as sim

KEY = ["rung", "score", "set", "k"]
Z_MODEL = stats.norm.ppf(0.95)


def clopper_pearson(x: int, n: int, a: float = 0.05) -> tuple[float, float]:
    lo = stats.beta.ppf(a / 2, x, n - x + 1) if x else 0.0
    hi = stats.beta.ppf(1 - a / 2, x + 1, n - x) if x < n else 1.0
    return float(lo), float(hi)


def lod_effective(lod: float, eff: pd.Series) -> float:
    """Map a nominal-f LoD to the effective share (eff: mean share indexed by f)."""
    if not np.isfinite(lod):
        return float("nan")
    x, y = np.log10(eff.index.to_numpy(float)), np.log10(eff.to_numpy(float))
    return float(10 ** np.interp(np.log10(lod), x, y))


def detect_table(e: pd.DataFrame, keys: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """FPR and detection/AUC per keys + f from eval rows with a 'detected' column."""
    fp = []
    for key, g in e[e.f == 0].groupby(keys):
        x, n = int(g.detected.sum()), len(g)
        fp.append((*key, x / n, x, n, *clopper_pearson(x, n)))
    fpr = pd.DataFrame(fp, columns=[*keys, "fpr", "fp", "n", "fpr_lo95", "fpr_hi95"])
    det = []
    for key, g in e[e.f > 0].groupby([*keys, "f"]):
        neg = e[(e.f == 0)].merge(pd.DataFrame([key[:-1]], columns=keys), on=keys)
        det.append((*key, g.detected.mean(), len(g), g.value.median(),
                    p5.auc(g.value, neg.value)))
    return fpr, pd.DataFrame(det, columns=[*keys, "f", "rate", "n", "median_score", "auc"])


def lods(det: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    return pd.DataFrame([(*key, sim.lod_interp(g.f, g.rate)) for key, g in
                         det.groupby(keys)], columns=[*keys, "lod"])


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/phase6"
    th = pd.read_csv(out / "threshold_scores.tsv", sep="\t")
    ev = pd.read_csv(out / "eval_scores.tsv.gz", sep="\t")
    eff = pd.read_csv(out / "effective_share.tsv", sep="\t")

    thr = [(*key, "bound", p5.prediction_bound_pools(g.value.dropna().to_numpy(), key[-1]),
            g.value.notna().sum()) for key, g in th.groupby(KEY)]
    thr += [(r, "FRAG", s, k, "model", Z_MODEL, 0) for r, s, k in
            th[th.score == "FRAG"][["rung", "set", "k"]].drop_duplicates().itertuples(
                index=False)]
    thr = pd.DataFrame(thr, columns=[*KEY, "rule", "threshold", "n_pools"])
    thr.to_csv(out / "thresholds.tsv", sep="\t", index=False, float_format="%.6g")

    e = ev[ev.kind == "eval"].merge(thr, on=KEY)
    e["detected"] = e.value > e.threshold
    k6 = [*KEY, "rule"]
    fpr, det = detect_table(e, k6)
    lod = lods(det, k6)
    share = {"WZ": "share_calls", "FRAG": "share_frag_n3_selected"}
    mean_eff = eff.groupby(["set", "k", "f"])[list(share.values())].mean()
    lod["lod_effective"] = [
        lod_effective(r.lod, mean_eff.loc[(r.set, r.k)][share[r.score]])
        for r in lod.itertuples()]
    lod = lod.merge(fpr, on=k6)
    fpr.to_csv(out / "fpr.tsv", sep="\t", index=False, float_format="%.4g")
    det.to_csv(out / "detection.tsv", sep="\t", index=False, float_format="%.4g")
    lod.to_csv(out / "lod.tsv", sep="\t", index=False, float_format="%.4g")

    # Matched FPR: Phase 5 (real rung) and Phase 6 scores, threshold = max eval null.
    p5d = pd.read_csv(ROOT / "results/phase5/scores.tsv.gz", sep="\t")
    p5d = p5d[(p5d.rung == "real") & (p5d.kind == "eval")]
    p5d = p5d[p5d.f_alt.isna() | np.isclose(p5d.f_alt, p5d.f) | (p5d.f == 0)].copy()
    p5d["f_alt"] = p5d.f_alt.fillna(-1.0)
    p5d["phase"] = 5
    p6d = ev[ev.kind == "eval"].assign(f_alt=-1.0, phase=6)
    allx = pd.concat([p5d, p6d], ignore_index=True)
    mrows = []
    for (ph, rung, score, s, k), g in allx.groupby(["phase", *KEY]):
        for f in sorted(g.f[g.f > 0].unique()):
            fa = f if score.startswith("LLR") else -1.0
            neg = g[(g.f == 0) & np.isclose(g.f_alt, fa)].value
            pos = g[(g.f == f) & np.isclose(g.f_alt, fa)].value
            mrows.append((ph, rung, score, s, k, f, float((pos > neg.max()).mean()),
                          len(neg), p5.auc(pos, neg)))
    m = pd.DataFrame(mrows, columns=["phase", *KEY, "f", "rate_at_max_null", "n_null",
                                     "auc"])
    mlod = (m.groupby(["phase", *KEY])
            .apply(lambda g: sim.lod_interp(g.f, g.rate_at_max_null), include_groups=False)
            .rename("lod_at_max_null").reset_index())
    m = m.merge(mlod, on=["phase", *KEY])
    m.to_csv(out / "matched.tsv", sep="\t", index=False, float_format="%.4g")

    # External check: 10 cancers vs the 16 eval donors, k = 1.
    c = ev[ev.kind == "cancer"].merge(thr, on=KEY)
    c["detected"] = c.value > c.threshold
    neg = e[(e.f == 0) & (e.k == 1)]
    ext = []
    for key, g in c.groupby(k6):
        n = neg.merge(pd.DataFrame([key], columns=k6), on=k6)
        ext.append((*key, len(g), int(g.detected.sum()), len(n), int(n.detected.sum()),
                    p5.auc(g.value, n.value)))
    pd.DataFrame(ext, columns=[*k6, "n_cancer", "cancer_detected", "n_healthy",
                               "healthy_detected", "auc"]).to_csv(
        out / "external.tsv", sep="\t", index=False, float_format="%.4g")
    print(lod[lod.rung == "real"].to_string(index=False))


if __name__ == "__main__":
    main()
