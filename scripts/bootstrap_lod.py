"""Donor bootstrap: 95 % intervals for detection rates, FPR and LoD of the real-read tests.

The 5 seeds of a test redraw only the spike fragments; the background reads are the same. So
the independent unit is the eval pool (one donor at k = 1, two disjoint donors at k = 2), not
the test. This script resamples eval pools with replacement and keeps all seeds of a pool
together. Thresholds stay fixed: they come from the 14 threshold donors, which are not
resampled. Each replicate recomputes the detection rate per f, the FPR (f = 0) and the LoD
(src/simulate.lod_interp, as in the summaries). Intervals are percentile intervals. A
replicate whose rate stays below 95 % at the top of the grid gives a censored LoD; when more
than 2.5 % of replicates are censored, the upper bound is reported as censored (NaN).

The interpolated LoD assumes the rate is linear in log f between grid points. The grid has no
point between 3 % and 10 %, so the interpolated interval shows donor sampling only, not the
grid gap. The grid interval does not interpolate: in each replicate, the LoD lies between the
last grid f below 95 % and the next grid f. The grid interval runs from the 2.5th percentile
of the lower ends to the 97.5th percentile of the upper ends. It is the interval to quote.
The FPR interval is Clopper-Pearson over pools (one spike-free test per pool), because a
bootstrap of 0 of 16 gives a zero-width interval.

Arms (tumour marker set, k = 1 and 2, "bound" threshold rule unless noted):
  - hct116: FRAG and WZ (Phase 6, rung "real"); UXM and LLR_binom (Phase 5, rungs "real"
    and "binom_real", LLR scored at f_alt = f; the f = 0 FPR is the largest over f_alt);
  - patient: FRAG, every GSE186007 run used as the spike.

Every LoD is also given on the effective share, the pre-registered primary scale (amendment
A1 of pipeline/PREREGISTERED_PHASE6.md): the mean share of informative fragments in the
FRAG-selected blocks that come from the spike, per k and f. One mapping (the FRAG share) is
used for every score, so all curves of one arm share one x-axis. For a patient run it is the
run's effective share times its tumour fraction (TF): tumour DNA at the markers.

Reads results/phase5/{scores.tsv.gz, thresholds.tsv}, results/phase6/{eval_scores.tsv.gz,
thresholds.tsv, effective_share.tsv}, results/phase6_gse/{eval_scores.tsv.gz,
effective_share.tsv, tf.tsv}. Writes results/bootstrap/{detection.tsv, lod.tsv,
spike_fragments.tsv}. In detection.tsv, x_spike is the effective share of the spike and
x_marker the tumour DNA at the markers (x_spike times TF; equal for HCT116).

Usage: uv run python scripts/bootstrap_lod.py [n_boot]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import simulate as sim

RES = ROOT / "results"
OUT = RES / "bootstrap"
SEED = 20261014
N_BOOT = 2000
KS = (1, 2)
SET = "tumour"
SHARE = "share_frag_n3_selected"


def clopper_pearson(x: int, n: int, a: float = 0.05) -> tuple[float, float]:
    lo = stats.beta.ppf(a / 2, x, n - x + 1) if x else 0.0
    hi = stats.beta.ppf(1 - a / 2, x + 1, n - x) if x < n else 1.0
    return float(lo), float(hi)


def lod_effective(lod: float, eff: pd.Series) -> float:
    """Map a nominal-f LoD to the effective share on log-log axes (as scripts/phase6_summary)."""
    if not np.isfinite(lod):
        return float("nan")
    x, y = np.log10(eff.index.to_numpy(float)), np.log10(eff.to_numpy(float))
    return float(10 ** np.interp(np.log10(lod), x, y))


def curve_tests() -> pd.DataFrame:
    """One row per eval test: arm, spike, score, rule, k, pool, f, seed, detected."""
    p6 = pd.read_csv(RES / "phase6/eval_scores.tsv.gz", sep="\t")
    p6 = p6[(p6.kind == "eval") & (p6.rung == "real") & (p6.set == SET) & p6.k.isin(KS)]
    p6 = p6.assign(pool=p6.pool.astype(int), spike="HCT116", arm="hct116", f_alt=-1.0)
    t6 = pd.read_csv(RES / "phase6/thresholds.tsv", sep="\t")
    t6 = t6[t6.rung == "real"].drop(columns=["n_pools"])
    a = p6.merge(t6, on=["rung", "score", "set", "k"])

    p5 = pd.read_csv(RES / "phase5/scores.tsv.gz", sep="\t")
    p5 = p5[(p5.kind == "eval") & p5.rung.isin(["real", "binom_real"]) & (p5.set == SET)
            & p5.score.isin(["UXM", "LLR_binom"]) & p5.k.isin(KS)]
    p5 = p5.assign(pool=p5.pool.astype(int), f_alt=p5.f_alt.fillna(-1.0), spike="HCT116",
                   arm="hct116", rule="bound")
    t5 = pd.read_csv(RES / "phase5/thresholds.tsv", sep="\t").drop(columns=["n_pools"])
    b = p5.merge(t5, on=["rung", "score", "set", "f_alt", "k"])
    # Spiked tests: LLR only at f_alt = f. Spike-free tests: every f_alt (FPR per f_alt).
    b = b[(b.f == 0) | (b.f_alt == -1) | np.isclose(b.f_alt, b.f)]

    g = pd.read_csv(RES / "phase6_gse/eval_scores.tsv.gz", sep="\t")
    g = g[(g.set == SET) & (g.score == "FRAG") & g.k.isin(KS)]
    null = g[g.spike == "none"]
    g = pd.concat([g[g.spike != "none"]] + [null.assign(spike=s) for s in
                                             g.spike[g.spike != "none"].unique()])
    g = g.assign(rung="real", arm="patient", f_alt=-1.0, pool=g.pool.astype(int))
    c = g.merge(t6[t6.score == "FRAG"], on=["rung", "score", "set", "k"])

    d = pd.concat([a, b, c], ignore_index=True)
    # A (curve, f) cell with an unscored test is dropped, as in the summaries.
    cur = ["arm", "spike", "rung", "score", "rule", "k"]
    bad = d.groupby([*cur, "f"]).value.transform(lambda v: v.isna().any())
    d = d[~bad].copy()
    d["detected"] = d.value > d.threshold
    return d[[*cur, "f_alt", "pool", "f", "seed", "detected"]]


def boot_curve(g: pd.DataFrame, rng: np.random.Generator, n_boot: int):
    """Point and bootstrap rates per f (f = 0: max FPR over f_alt), and the LoD per replicate."""
    pools = np.sort(g.pool.unique())
    cols = [(f, fa) for f, fa in sorted(set(zip(g.f, g.f_alt, strict=True)))]
    hits = np.zeros((len(pools), len(cols)))
    tests = np.zeros_like(hits)
    pi = {p: i for i, p in enumerate(pools)}
    ci = {c: j for j, c in enumerate(cols)}
    for (p, f, fa), h in g.groupby(["pool", "f", "f_alt"]).detected:
        hits[pi[p], ci[(f, fa)]] = h.sum()
        tests[pi[p], ci[(f, fa)]] = len(h)
    idx = np.vstack([np.arange(len(pools)), rng.integers(0, len(pools), (n_boot, len(pools)))])
    rate = hits[idx].sum(1) / tests[idx].sum(1)  # (1 + n_boot, n_cols)
    fs = sorted({f for f, _ in cols})
    by_f = np.stack([rate[:, [j for j, (f, _) in enumerate(cols) if f == ff]].max(1)
                     for ff in fs], axis=1)
    pos = [i for i, f in enumerate(fs) if f > 0]
    fpos = np.array(fs)[pos]
    lods = np.array([sim.lod_interp(fpos, r[pos]) for r in by_f])
    return fs, by_f, lods, np.array([bracket(fpos, r[pos]) for r in by_f]), len(pools)


def bracket(levels: np.ndarray, rates: np.ndarray, target: float = 0.95):
    """Grid points around the LoD: (last f below target, next f); 0 below the grid, inf above."""
    if rates[-1] < target:
        return levels[-1], np.inf
    below = np.where(rates < target)[0]
    if len(below) == 0:
        return 0.0, levels[0]
    return levels[below[-1]], levels[below[-1] + 1]


def interval(x: np.ndarray) -> tuple[float, float, float]:
    """2.5 and 97.5 percentiles, censored (NaN) values ranked above every finite value."""
    v = np.where(np.isfinite(x), x, np.inf)
    with np.errstate(invalid="ignore"):
        lo, hi = np.percentile(v, [2.5, 97.5], method="inverted_cdf")
    share = float(np.mean(~np.isfinite(x)))
    return (float(lo) if np.isfinite(lo) else np.nan, float(hi) if np.isfinite(hi) else np.nan,
            share)


def main() -> None:
    n_boot = int(sys.argv[1]) if len(sys.argv) > 1 else N_BOOT
    OUT.mkdir(exist_ok=True)
    tests = curve_tests()
    e6 = pd.read_csv(RES / "phase6/effective_share.tsv", sep="\t")
    e6 = e6[(e6.rung == "real") & (e6.set == SET) & (e6.f > 0)]
    eff_h = e6.groupby(["k", "f"])[SHARE].mean()
    eg = pd.read_csv(RES / "phase6_gse/effective_share.tsv", sep="\t")
    eg = eg[(eg.set == SET) & (eg.f > 0)]
    eff_g = eg.groupby(["spike", "k", "f"])[SHARE].mean()
    tf = pd.read_csv(RES / "phase6_gse/tf.tsv", sep="\t").set_index("run").tf

    det, lod = [], []
    cur = ["arm", "spike", "rung", "score", "rule", "k"]
    for i, (key, g) in enumerate(tests.groupby(cur, sort=True)):
        arm, spike, _, _, _, k = key
        rng = np.random.default_rng([SEED, i])
        fs, by_f, lods, br, n_pools = boot_curve(g, rng, n_boot)
        if arm == "hct116":
            eff, scale = eff_h.loc[k], 1.0
        else:
            eff, scale = eff_g.loc[(spike, k)], float(tf[spike])
        for j, f in enumerate(fs):
            if f > 0:
                lo, hi = np.percentile(by_f[1:, j], [2.5, 97.5])
            else:  # Clopper-Pearson: one spike-free test per pool (per f_alt)
                lo, hi = clopper_pearson(round(by_f[0, j] * n_pools), n_pools)
            x = float(eff.loc[f]) if f > 0 else 0.0
            det.append((*key, f, x, x * scale, by_f[0, j], lo, hi, n_pools))
        lo, hi, share = interval(lods[1:])
        g_lo = float(np.percentile(br[1:, 0], 2.5, method="inverted_cdf"))
        g_hi = float(np.percentile(br[1:, 1], 97.5, method="inverted_cdf"))
        g_hi = g_hi if np.isfinite(g_hi) else np.nan
        e = [0.0 if v == 0 else lod_effective(v, eff) * scale
             for v in (lods[0], lo, hi, g_lo, g_hi)]
        x0 = round(by_f[0, 0] * n_pools)
        lod.append((*key, n_pools, lods[0], lo, hi, g_lo, g_hi, share, scale, *e,
                    x0 / n_pools, x0, *clopper_pearson(x0, n_pools)))
    det = pd.DataFrame(det, columns=[*cur, "f", "x_spike", "x_marker", "rate", "rate_lo95",
                                     "rate_hi95", "n_pools"])
    lod = pd.DataFrame(lod, columns=[
        *cur, "n_pools", "lod", "lod_lo95", "lod_hi95", "grid_lo95", "grid_hi95",
        "share_censored", "tf", "lod_marker", "lod_marker_lo95", "lod_marker_hi95",
        "grid_marker_lo95", "grid_marker_hi95", "fpr", "fp", "fpr_lo95", "fpr_hi95"])
    det.to_csv(OUT / "detection.tsv", sep="\t", index=False, float_format="%.4g")
    lod.to_csv(OUT / "lod.tsv", sep="\t", index=False, float_format="%.4g")

    # Tumour-like spike fragments per test (HCT116, FRAG-selected blocks): mean and range.
    sf = (e6.assign(x_marker=e6[SHARE]).groupby(["k", "f"])
          .agg(x_marker=(SHARE, "mean"), mean=("spike_tumour_like_selected", "mean"),
               lo=("spike_tumour_like_selected", "min"),
               hi=("spike_tumour_like_selected", "max")).reset_index())
    sf.to_csv(OUT / "spike_fragments.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 200)
    show = lod[(lod.rule == "bound") & ((lod.arm == "hct116") | (lod.tf > 0))]
    print(show.drop(columns=["arm", "rule"]).to_string(index=False))


if __name__ == "__main__":
    main()
