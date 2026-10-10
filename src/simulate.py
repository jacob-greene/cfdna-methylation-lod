"""Sequencing-depth simulation of tumour detection at marker CpGs (Phase 3).

Model: at each CpG, the methylated read count of a sample at tumour fraction f is
Binomial(n, f * beta_tumour + (1 - f) * beta_background), where n is the depth D, or the
number of unique input molecules when that is smaller (the input cap).

Two untrained detection scores:
  - NNLS colon fraction of the simulated betas (k / n) against the reference atlas;
  - summed binomial log-likelihood ratio (LLR) of "tumour at fraction f" against
    "background only", with the background profile taken from training blanks.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

EPS = 1e-3  # probabilities are clipped to [EPS, 1 - EPS] inside the LLR


def demix(beta_plasma: np.ndarray, beta_background: np.ndarray, c: float) -> np.ndarray:
    """Pure-tumour profile from a cancer plasma whose tumour fraction is c.

    Solves beta_plasma = c * beta_tumour + (1 - c) * beta_background, clipped to [0, 1].
    """
    return np.clip((beta_plasma - (1.0 - c) * beta_background) / c, 0.0, 1.0)


def rank_markers(beta_tumour: np.ndarray, train: np.ndarray, sd_floor: float = 0.02
                 ) -> np.ndarray:
    """CpG indices ordered by tumour-vs-blank contrast, strongest first.

    Score = |beta_tumour - mean(train)| / sqrt(sd(train)^2 + sd_floor^2), where train is a
    blank x CpG matrix of training blanks only. The SD term down-ranks CpGs that vary
    between people.
    """
    mu = train.mean(axis=0)
    sd = train.std(axis=0, ddof=1) if len(train) > 1 else np.zeros_like(mu)
    score = np.abs(beta_tumour - mu) / np.sqrt(sd ** 2 + sd_floor ** 2)
    return np.argsort(-score, kind="stable")


def molecules(depth: int, cap_mean: float | None, n_cpg: int, size: int,
              rng: np.random.Generator) -> np.ndarray:
    """Informative reads per CpG: D, or min(D, Poisson(cap_mean)) unique molecules.

    Reads beyond the unique-molecule count are PCR duplicates and carry no new information.
    Returns an int array of shape (size, n_cpg); values are at least 1.
    """
    if cap_mean is None:
        return np.full((size, n_cpg), depth, dtype=np.int64)
    m = rng.poisson(cap_mean, size=(size, n_cpg))
    return np.maximum(np.minimum(m, depth), 1)


def simulate_counts(beta_background: np.ndarray, beta_tumour: np.ndarray | None, f: float,
                    n: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Methylated read counts, shape like n: Binomial(n, f * tumour + (1 - f) * background)."""
    p = beta_background if beta_tumour is None or f == 0 else (
        f * beta_tumour + (1.0 - f) * beta_background)
    return rng.binomial(n, np.clip(p, 0.0, 1.0))


def llr(k: np.ndarray, n: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> np.ndarray:
    """Summed binomial LLR over CpGs (last axis) of p1 against p0."""
    p0 = np.clip(p0, EPS, 1 - EPS)
    p1 = np.clip(p1, EPS, 1 - EPS)
    w1 = np.log(p1 / p0)
    w0 = np.log((1 - p1) / (1 - p0))
    return (k * w1 + (n - k) * w0).sum(axis=-1)


def nnls_colon(betas: np.ndarray, atlas: np.ndarray, i_colon: int) -> np.ndarray:
    """NNLS colon fraction (coefficients rescaled to sum 1) of each row of betas."""
    out = np.empty(len(betas))
    for i, y in enumerate(betas):
        coef, _ = optimize.nnls(atlas, y)
        s = coef.sum()
        out[i] = coef[i_colon] / s if s > 0 else 0.0
    return out


def prediction_bound(per_person: list[np.ndarray], q: float = 0.95) -> float:
    """One-sided normal prediction bound for one new result from a new person.

    per_person: one array of null scores per training person. Between-person spread comes
    from the per-person means; within-person spread is the pooled SD around them. A new
    person's single result exceeds the bound with probability 1 - q under normality.
    """
    from scipy import stats

    means = np.array([np.mean(s) for s in per_person])
    n = len(means)
    within = np.mean([np.var(s, ddof=1) for s in per_person if len(s) > 1] or [0.0])
    sd = np.sqrt(np.var(means, ddof=1) * (1 + 1 / n) + within)
    return float(means.mean() + stats.t.ppf(q, n - 1) * sd)


def lod_from_rates(levels, rates, target: float = 0.95) -> float:
    """Lowest level with rate >= target at it and at every higher level (NaN if none)."""
    levels = np.asarray(levels, dtype=float)
    rates = np.asarray(rates, dtype=float)
    order = np.argsort(levels)
    levels, rates = levels[order], rates[order]
    for i in range(len(levels)):
        if (rates[i:] >= target).all():
            return float(levels[i])
    return float("nan")


def lod_interp(levels, rates, target: float = 0.95) -> float:
    """f where the detection rate first crosses target, interpolated on log10 f.

    Uses the last crossing from below, so a non-monotone rate curve gives the conservative
    value. NaN if the top level stays below target; the lowest level if all pass.
    """
    levels = np.asarray(levels, dtype=float)
    rates = np.asarray(rates, dtype=float)
    order = np.argsort(levels)
    x, r = np.log10(levels[order]), rates[order]
    if r[-1] < target:
        return float("nan")
    below = np.where(r < target)[0]
    if len(below) == 0:
        return float(levels[order][0])
    i = below[-1]
    t = (target - r[i]) / (r[i + 1] - r[i])
    return float(10 ** (x[i] + t * (x[i + 1] - x[i])))
