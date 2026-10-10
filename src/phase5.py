"""Phase 5 helpers: fragment records from wgbstools pat files, read-level mixing and scores.

A fragment record is one (pat row, block) pair: the methylated (m) and unmethylated (u) CpG
calls that one fragment has inside one block, with the row's multiplicity (w, the pat count
column). A fragment that spans two blocks of one set gives two records with the same row id.

Mixing (pipeline/PREREGISTERED.md): f is the share of fragments that come from the spike
source. Background fragments of a pool of k donors are kept whole; each spike fragment is kept
with probability q = f / (1 - f) * N_bg / N_spike, where N counts the fragments in the
off-marker blocks (2,000 random autosomal windows, a genome-wide depth proxy).

Scores use methylation calls at marker blocks only:
  - binomial LLR: sum over blocks of M_b log(p1/p0) + U_b log((1-p1)/(1-p0)),
    p1 = p0 + f_alt (T_b - p0), T_b = mean TCGA-COAD tumour beta of the block's probes;
  - beta-binomial LLR: each record's calls in a block are one beta-binomial draw with mean
    p0 or p1 and intra-fragment correlation rho (estimated on threshold donors);
  - UXM-style: tumour blocks, share of records with >= 3 calls that are all methylated.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

EPS = 1e-3  # probabilities are clipped to [EPS, 1 - EPS], as in src/simulate.llr
UXM_MIN_CALLS = 3


@dataclass
class Records:
    """Fragment records of one source for one block set (parallel arrays)."""

    row: np.ndarray  # pat row id (fragment identity for sampling)
    block: np.ndarray  # block index within the set
    m: np.ndarray
    u: np.ndarray
    w: np.ndarray  # multiplicity (fragments)

    @staticmethod
    def concat(parts: list[Records]) -> Records:
        return Records(*(np.concatenate([getattr(p, k) for p in parts])
                         for k in ("row", "block", "m", "u", "w")))

    def take(self, mask: np.ndarray, w: np.ndarray | None = None) -> Records:
        return Records(self.row[mask], self.block[mask], self.m[mask], self.u[mask],
                       self.w[mask] if w is None else w)


def read_pat(path: Path) -> pd.DataFrame:
    """wgbstools pat file: chrom, start CpG index (1-based), pattern, count."""
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 2, 3],
                         names=["chrom", "start", "pat", "count"],
                         dtype={"chrom": str, "start": np.int64, "pat": str, "count": np.int64})
    return df


def pat_records(pat: pd.DataFrame, starts: np.ndarray, ends: np.ndarray) -> Records:
    """Records of the pat rows that overlap the blocks [starts, ends) (CpG indices).

    Blocks must be sorted and disjoint. Blocks are indexed by their position in the input.
    """
    order = np.argsort(starts, kind="stable")
    s_sorted, e_sorted = starts[order], ends[order]
    assert (s_sorted[1:] >= e_sorted[:-1]).all(), "blocks overlap"
    st = pat.start.to_numpy()
    ln = pat.pat.str.len().to_numpy()
    first = np.searchsorted(e_sorted, st, side="right")  # first block whose end > start
    last = np.searchsorted(s_sorted, st + ln, side="left") - 1  # last block with start < end
    rows, blocks, ms, us, ws = [], [], [], [], []
    pats, counts = pat.pat.to_numpy(), pat["count"].to_numpy()
    for i in np.nonzero(last >= first)[0]:
        p, s = pats[i], st[i]
        for j in range(first[i], last[i] + 1):
            seg = p[max(s_sorted[j], s) - s: min(e_sorted[j], s + ln[i]) - s]
            m, u = seg.count("C"), seg.count("T")
            if m + u:
                rows.append(i)
                blocks.append(order[j])
                ms.append(m)
                us.append(u)
                ws.append(counts[i])
    a = lambda x: np.asarray(x, dtype=np.int64)
    return Records(a(rows), a(blocks), a(ms), a(us), a(ws))


def n_fragments(rec: Records) -> int:
    """Fragments with at least one call in the set (a row counts once even if it spans two)."""
    _, first = np.unique(rec.row, return_index=True)
    return int(rec.w[first].sum())


def n_fragments_pool(parts: list[Records]) -> int:
    """Fragments of a donor pool: the sum over donors. Pat row ids repeat across donors, so
    n_fragments(Records.concat(parts)) would merge different donors' fragments."""
    return sum(n_fragments(p) for p in parts)


def thin_rows(counts: np.ndarray, q: float, rng: np.random.Generator) -> np.ndarray:
    """Fragments kept per pat row when each fragment is kept with probability q."""
    if not 0 <= q <= 1:
        raise ValueError(f"spike sampling probability {q:.4f} is outside [0, 1]")
    return rng.binomial(counts, q)


def keep_rows(rec: Records, kept: np.ndarray) -> Records:
    """Records of the kept fragments; kept = fragments kept per pat row (thin_rows)."""
    k = kept[rec.row]
    return rec.take(k > 0, k[k > 0])


def block_sums(rec: Records, n_blocks: int) -> tuple[np.ndarray, np.ndarray]:
    """Methylated and unmethylated calls per block."""
    M = np.bincount(rec.block, weights=rec.w * rec.m, minlength=n_blocks)
    U = np.bincount(rec.block, weights=rec.w * rec.u, minlength=n_blocks)
    return M, U


def p_alt(p0: np.ndarray, tumour: np.ndarray, f: float) -> np.ndarray:
    return p0 + f * (tumour - p0)


def llr_binomial(M: np.ndarray, U: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> float:
    p0, p1 = np.clip(p0, EPS, 1 - EPS), np.clip(p1, EPS, 1 - EPS)
    return float((M * np.log(p1 / p0) + U * np.log((1 - p1) / (1 - p0))).sum())


def _bb_logpmf_core(m, n, p, rho):
    """Beta-binomial log pmf without the binomial coefficient (it cancels in an LLR)."""
    p = np.clip(p, EPS, 1 - EPS)
    s = (1 - rho) / rho
    a, b = p * s, (1 - p) * s
    return special.betaln(m + a, n - m + b) - special.betaln(a, b)


def llr_betabinomial(rec: Records, p0: np.ndarray, p1: np.ndarray, rho: float) -> float:
    n = rec.m + rec.u
    pb0, pb1 = p0[rec.block], p1[rec.block]
    d = _bb_logpmf_core(rec.m, n, pb1, rho) - _bb_logpmf_core(rec.m, n, pb0, rho)
    return float((rec.w * d).sum())


def fit_rho(rec: Records, p0: np.ndarray) -> float:
    """MLE of one intra-fragment correlation rho, block means fixed at p0 (records n >= 2)."""
    n = rec.m + rec.u
    keep = n >= 2
    m, n, w, p = rec.m[keep], n[keep], rec.w[keep], p0[rec.block[keep]]

    def nll(lr):
        rho = 1 / (1 + np.exp(-lr))
        return -(w * _bb_logpmf_core(m, n, p, rho)).sum()

    res = optimize.minimize_scalar(nll, bounds=(-10, 5), method="bounded")
    return float(1 / (1 + np.exp(-res.x)))


def uxm_share(rec: Records, blocks_mask: np.ndarray | None = None) -> float:
    """Share of records with >= 3 calls whose calls are all methylated (weighted)."""
    n = rec.m + rec.u
    keep = n >= UXM_MIN_CALLS
    if blocks_mask is not None:
        keep &= blocks_mask[rec.block]
    tot = rec.w[keep].sum()
    return float(rec.w[keep & (rec.u == 0)].sum() / tot) if tot else np.nan


def prediction_bound_pools(null_scores: np.ndarray, k: int, n_donors: int = 14,
                           q: float = 0.95) -> float:
    """Pre-registered threshold: one-sided normal prediction bound over threshold pools of k.

    mean + t(q, df) * sd * sqrt(1 + 1 / n_eff), df = max(floor(n_donors / k) - 1, 1),
    n_eff = df + 1 (the number of disjoint pools the donors can form).
    """
    s = np.asarray(null_scores, dtype=float)
    df = max(n_donors // k - 1, 1)
    return float(s.mean() + stats.t.ppf(q, df) * s.std(ddof=1) * np.sqrt(1 + 1 / (df + 1)))


def auc(pos, neg) -> float:
    """Mann-Whitney AUC; ties count one half."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    r = stats.rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))
