"""Phase 6 helpers: a null model fitted on threshold donors, and two optimised scores.

All fitted values come from per-donor block summaries of the threshold donors
(pipeline/PREREGISTERED_PHASE6.md). A summary of one sample and one block set holds, per block:
  M, N  methylated calls, all calls;
  G     sum over fragments of n (1 + (n - 1) rho): the variance of M under within-fragment
        co-methylation is p (1 - p) G (rho = Phase 5 intra-fragment correlation);
  cM    fragments with >= 3 calls, all methylated;
  cU    fragments with >= 3 calls, all unmethylated;
  n3    fragments with >= 3 calls.
Summaries add over donors and over the spike fragments of a mix.

Scores:
  WZ   weighted z-score: sum_b D_b (x_b - p0_b) / v_b / sqrt(sum_b D_b^2 / v_b),
       D_b = T_b - p0_b, v_b = q (1 - q) (G_b / N_b^2 + phi_b / k), q = clip(p0_b), phi_b the
       between-donor excess variance per block, shrunk to the set mean.
  FRAG tumour-like fragment count in selected blocks: Z = (C - E) / sqrt(E + s2 E^2 / k),
       C = sum c_b, E = sum n3_b r_b; c_b = cM_b where D_b > 0 else cU_b; r_b = background
       rate of tumour-like fragments; s2 = between-donor overdispersion (nested leave-one-out).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import phase5 as p5

FIELDS = ("M", "N", "G", "cM", "cU", "n3")
Q_CLIP = 0.01  # p0 is clipped to [Q_CLIP, 1 - Q_CLIP] in variances
PHI_NU = 14.0  # shrinkage weight (in donors) of the set-mean phi
MIN_DONORS_PHI = 3  # blocks with fewer covered donors use the set mean
SEL_DELTA = 0.3  # FRAG blocks: abs(T - p0) >= SEL_DELTA
SEL_RATE = 0.02  # FRAG blocks: background tumour-like rate <= SEL_RATE
RATE_PSEUDO = 10.0  # pseudo-fragments at the set-mean rate added to each block's rate
MIN_CALLS = p5.UXM_MIN_CALLS


def summarise(rec: p5.Records, nb: int, rho: float) -> np.ndarray:
    """Block summary (len(FIELDS) x nb) of one sample's records for one set."""
    n = rec.m + rec.u
    q3 = n >= MIN_CALLS
    w = rec.w.astype(float)
    cols = (w * rec.m, w * n, w * n * (1 + (n - 1) * rho),
            w * (q3 & (rec.u == 0)), w * (q3 & (rec.m == 0)), w * q3)
    return np.stack([np.bincount(rec.block, weights=c, minlength=nb) for c in cols])


def _get(S: np.ndarray, name: str) -> np.ndarray:
    return S[FIELDS.index(name)]


def tumour_like(S: np.ndarray, direction: np.ndarray) -> np.ndarray:
    return np.where(direction > 0, _get(S, "cM"), _get(S, "cU"))


@dataclass
class Model:
    p0: np.ndarray
    delta: np.ndarray
    phi: np.ndarray  # shrunk between-donor excess variance per block
    phi_mean: float
    sel: np.ndarray  # FRAG blocks
    rate: np.ndarray  # FRAG background rate per block
    s2: float  # FRAG overdispersion


def _p0(Ss: list[np.ndarray]) -> np.ndarray:
    M = sum(_get(S, "M") for S in Ss)
    N = sum(_get(S, "N") for S in Ss)
    return np.where(N > 0, M / np.maximum(N, 1), np.nan)


def _rate_sel(Ss, p0, T):
    delta = T - p0
    d = np.sign(delta)
    c = sum(tumour_like(S, d) for S in Ss)
    n3 = sum(_get(S, "n3") for S in Ss)
    cand = np.isfinite(p0) & (np.abs(delta) >= SEL_DELTA)
    rbar = c[cand].sum() / max(n3[cand].sum(), 1)
    rate = (c + RATE_PSEUDO * rbar) / (n3 + RATE_PSEUDO)
    return rate, cand & (rate <= SEL_RATE), d


def _phi(Ss, p0):
    """Per-block between-donor excess variance (scaled by q (1 - q)), leave-one-out p0."""
    M_all = sum(_get(S, "M") for S in Ss)
    N_all = sum(_get(S, "N") for S in Ss)
    ex = []
    for S in Ss:
        M, N, G = _get(S, "M"), _get(S, "N"), _get(S, "G")
        Nr = N_all - N
        q = np.clip((M_all - M) / np.maximum(Nr, 1), Q_CLIP, 1 - Q_CLIP)
        ok = (N > 0) & (Nr > 0)
        x = M / np.maximum(N, 1)
        ex.append(np.where(ok, (x - q) ** 2 / (q * (1 - q)) - G / np.maximum(N, 1) ** 2,
                           np.nan))
    ex = np.array(ex)
    n_cov = np.isfinite(ex).sum(0)
    raw = np.where(n_cov >= MIN_DONORS_PHI, np.nansum(ex, 0) / np.maximum(n_cov, 1), np.nan)
    ok = np.isfinite(raw) & np.isfinite(p0)
    phi_mean = float(max(raw[ok].mean(), 0.0))
    phi = np.where(ok, (n_cov * np.nan_to_num(raw) + PHI_NU * phi_mean) / (n_cov + PHI_NU),
                   phi_mean)
    return np.maximum(phi, 0.0), phi_mean


def frag_counts(S, model: Model) -> tuple[float, float]:
    """(C, E): tumour-like fragments in the selected blocks and their background expectation."""
    d = np.sign(model.delta)
    c = tumour_like(S, d)
    return float(c[model.sel].sum()), float((_get(S, "n3") * model.rate)[model.sel].sum())


def fit(Ss: list[np.ndarray], T: np.ndarray) -> Model:
    """Null model from the block summaries of single donors (threshold donors only)."""
    p0 = _p0(Ss)
    rate, sel, _ = _rate_sel(Ss, p0, T)
    phi, phi_mean = _phi(Ss, p0)
    # FRAG overdispersion: each donor scored with rates and selection fitted without it.
    terms = []
    for i, S in enumerate(Ss):
        rest = Ss[:i] + Ss[i + 1:]
        p0r = _p0(rest)
        rr, sr, _ = _rate_sel(rest, p0r, T)
        C, E = frag_counts(S, Model(p0r, T - p0r, phi, phi_mean, sr, rr, 0.0))
        if E > 0:
            terms.append(((C - E) ** 2 - E) / E ** 2)
    s2 = float(max(np.mean(terms), 0.0)) if terms else 0.0
    return Model(p0, T - p0, phi, phi_mean, sel, rate, s2)


def wz(S: np.ndarray, model: Model, k: int, M: np.ndarray | None = None) -> float:
    """Weighted z-score. M overrides the methylated calls (binomial rung, same weights)."""
    M = _get(S, "M") if M is None else M
    N, G = _get(S, "N"), _get(S, "G")
    ok = (N > 0) & np.isfinite(model.p0)
    q = np.clip(model.p0[ok], Q_CLIP, 1 - Q_CLIP)
    v = q * (1 - q) * (G[ok] / N[ok] ** 2 + model.phi[ok] / k)
    D = model.delta[ok]
    x = M[ok] / N[ok]
    return float((D * (x - model.p0[ok]) / v).sum() / np.sqrt((D ** 2 / v).sum()))


def frag(S: np.ndarray, model: Model, k: int) -> float:
    C, E = frag_counts(S, model)
    return float((C - E) / np.sqrt(E + model.s2 * E ** 2 / k)) if E > 0 else np.nan
