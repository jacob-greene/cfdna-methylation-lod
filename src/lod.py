"""Limit of blank (LoB) and limit of detection (LoD) helpers, following CLSI EP17-A2.

All fractions are on the 0-1 scale. The measurand is the NNLS colon fraction of an array.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize


def nnls_fractions(y: np.ndarray, a: np.ndarray) -> np.ndarray:
    """NNLS of one beta vector against an aligned atlas matrix (CpG x cell type).

    Rows where y or any atlas value is NaN are dropped; coefficients are rescaled to sum 1,
    exactly as `src.deconv.nnls_one` and the authors' deconvolve.py do.
    """
    keep = ~np.isnan(y) & ~np.isnan(a).any(axis=1)
    coef, _ = optimize.nnls(a[keep], y[keep])
    total = coef.sum()
    return coef / total if total > 0 else coef


def mix(background: np.ndarray, spike: np.ndarray, f: float) -> np.ndarray:
    """Beta-level in-silico mix: f * spike + (1 - f) * background (NaN if either is NaN)."""
    return f * spike + (1.0 - f) * background


def lob_nonparametric(values, q: float = 0.95) -> float:
    """CLSI EP17-A2 nonparametric LoB: the value at rank 0.5 + q * n of the sorted blanks.

    Ranks are 1-based; a fractional rank interpolates linearly; a rank above n returns the
    maximum. EP17 asks for at least 60 blank results; below that the estimate is near the
    maximum and its uncertainty is large.
    """
    v = np.sort(np.asarray(values, dtype=float))
    n = len(v)
    rank = 0.5 + q * n
    if rank >= n:
        return float(v[-1])
    lo = int(np.floor(rank))
    return float(v[lo - 1] + (rank - lo) * (v[lo] - v[lo - 1]))


def lob_by_spike(d, lob_backgrounds):
    """LoB per spike replicate from blank rows only.

    d: one row per measurement with columns spike, background, f and colon. Blank rows are
    those at f == 0 whose background is in lob_backgrounds; spiked rows and other
    backgrounds never enter the LoB.
    """
    blanks = d[(d.f == 0) & d.background.isin(lob_backgrounds)]
    return blanks.groupby("spike").colon.apply(lob_nonparametric)


def bootstrap_ci(values, stat, n_boot: int, rng: np.random.Generator,
                 level: float = 0.95) -> tuple[float, float]:
    """Percentile bootstrap CI of stat(values), resampling values with replacement."""
    v = np.asarray(values, dtype=float)
    boots = np.array([stat(rng.choice(v, size=len(v), replace=True)) for _ in range(n_boot)])
    a = (1 - level) / 2
    return float(np.quantile(boots, a)), float(np.quantile(boots, 1 - a))


def lod_parametric(lob: float, groups: list[np.ndarray], z: float = 1.645) -> dict:
    """CLSI EP17-A2 parametric LoD = LoB + c_p * SD_S.

    SD_S is the SD pooled over K low-level groups (N results in total), and
    c_p = z / (1 - 1 / (4 * (N - K))) corrects the SD estimate for small N.
    """
    groups = [np.asarray(g, dtype=float) for g in groups]
    n_total = sum(len(g) for g in groups)
    k = len(groups)
    ss = sum(((g - g.mean()) ** 2).sum() for g in groups)
    sd = float(np.sqrt(ss / (n_total - k)))
    cp = z / (1 - 1 / (4 * (n_total - k)))
    return {"lod": lob + cp * sd, "sd_pooled": sd, "c_p": cp, "n": n_total, "k": k}


def empirical_lod(levels, rates, target: float = 0.95) -> float:
    """Lowest level whose detection rate is >= target and stays >= target at every higher level.

    Returns NaN when no level qualifies.
    """
    levels = np.asarray(levels, dtype=float)
    rates = np.asarray(rates, dtype=float)
    order = np.argsort(levels)
    levels, rates = levels[order], rates[order]
    ok = rates >= target
    for i in range(len(levels)):
        if ok[i:].all():
            return float(levels[i])
    return float("nan")
