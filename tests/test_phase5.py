import numpy as np
import pandas as pd
import pytest
from scipy import stats

from src import phase5 as p5
from src import simulate as sim


def _pat(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "pat", "count"])


def test_pat_records_splits_calls_by_block_and_skips_missing():
    starts, ends = np.array([20, 10]), np.array([25, 15])  # unsorted input order
    pat = _pat([
        ("chr1", 13, "CCT.TT", 1),   # CpGs 13-18: block [10,15) gets 13,14 = "CC"
        ("chr1", 18, "TTCCCT", 2),   # CpGs 18-23: block [20,25) gets 20-23 = "CCCT"
        ("chr1", 14, "TCCCCCCT", 1),  # spans both: 14 = "T"; 20,21 = "CT"
        ("chr1", 30, "CCCC", 1),     # no block
        ("chr1", 9, "C.", 1),        # CpG 10 has no call
    ])
    r = p5.pat_records(pat, starts, ends)
    got = sorted(zip(r.row, r.block, r.m, r.u, r.w))
    assert got == [(0, 1, 2, 0, 1), (1, 0, 3, 1, 2), (2, 0, 1, 1, 1), (2, 1, 0, 1, 1)]
    assert p5.n_fragments(r) == 4  # row 2 counts once; row 1 has count 2


def test_pat_records_rejects_overlapping_blocks():
    with pytest.raises(AssertionError):
        p5.pat_records(_pat([("chr1", 1, "C", 1)]), np.array([0, 5]), np.array([10, 12]))


def test_thinning_keeps_rows_together_and_scales_with_q():
    rec = p5.Records(row=np.array([0, 0, 1]), block=np.array([0, 1, 1]),
                     m=np.array([1, 1, 2]), u=np.array([0, 1, 0]), w=np.array([5, 5, 3]))
    counts = np.array([5, 3])
    rng = np.random.default_rng(1)
    assert len(p5.keep_rows(rec, p5.thin_rows(counts, 0.0, rng)).row) == 0
    full = p5.keep_rows(rec, p5.thin_rows(counts, 1.0, rng))
    assert (full.w == rec.w).all()
    kept = np.array([2, 0])
    part = p5.keep_rows(rec, kept)
    assert list(part.row) == [0, 0] and list(part.w) == [2, 2]
    big = np.full(10000, 10)
    assert abs(p5.thin_rows(big, 0.3, rng).sum() / 1e5 - 0.3) < 0.01
    with pytest.raises(ValueError):
        p5.thin_rows(counts, 1.2, rng)


def test_binomial_llr_matches_phase3_llr():
    rng = np.random.default_rng(2)
    p0 = rng.uniform(0.05, 0.95, 50)
    p1 = p5.p_alt(p0, rng.uniform(0, 1, 50), 0.1)
    n = rng.integers(0, 6, 50)
    k = rng.binomial(n, p0)
    assert np.isclose(p5.llr_binomial(k, n - k, p0, p1), sim.llr(k, n, p0, p1))
    assert p5.llr_binomial(k, n - k, p0, p0) == 0


def test_betabinomial_llr_tends_to_binomial_as_rho_vanishes():
    rng = np.random.default_rng(3)
    nb = 30
    p0 = rng.uniform(0.1, 0.9, nb)
    p1 = p5.p_alt(p0, np.ones(nb), 0.2)
    rec = p5.Records(row=np.arange(200), block=rng.integers(0, nb, 200),
                     m=rng.integers(0, 3, 200), u=rng.integers(0, 3, 200), w=np.ones(200, int))
    M, U = p5.block_sums(rec, nb)
    assert np.isclose(p5.llr_betabinomial(rec, p0, p1, 1e-8), p5.llr_binomial(M, U, p0, p1),
                      rtol=1e-5)
    # n = 1 records do not depend on rho
    one = rec.take(rec.m + rec.u == 1)
    assert np.isclose(p5.llr_betabinomial(one, p0, p1, 0.5), p5.llr_betabinomial(one, p0, p1,
                                                                                  0.01))


def test_fit_rho_recovers_simulated_correlation():
    rng = np.random.default_rng(4)
    nb, n_rec, rho = 100, 20000, 0.3
    p0 = rng.uniform(0.1, 0.9, nb)
    blk = rng.integers(0, nb, n_rec)
    s = (1 - rho) / rho
    pf = rng.beta(p0[blk] * s, (1 - p0[blk]) * s)
    n = rng.integers(1, 7, n_rec)
    m = rng.binomial(n, pf)
    rec = p5.Records(np.arange(n_rec), blk, m, n - m, np.ones(n_rec, int))
    assert abs(p5.fit_rho(rec, p0) - rho) < 0.03


def test_uxm_share_counts_only_fully_methylated_records_with_three_calls():
    rec = p5.Records(row=np.arange(4), block=np.zeros(4, int), m=np.array([3, 2, 4, 1]),
                     u=np.array([0, 1, 0, 0]), w=np.array([1, 1, 2, 5]))
    assert p5.uxm_share(rec) == pytest.approx(3 / 4)  # record 3 has < 3 calls
    assert np.isnan(p5.uxm_share(rec.take(np.array([False, False, False, True]))))


def test_prediction_bound_pools_uses_preregistered_df():
    s = np.array([1.0, 2.0, 3.0, 4.0])
    sd = s.std(ddof=1)
    assert p5.prediction_bound_pools(s, k=1) == pytest.approx(
        2.5 + stats.t.ppf(0.95, 13) * sd * np.sqrt(1 + 1 / 14))
    assert p5.prediction_bound_pools(s, k=8) == pytest.approx(
        2.5 + stats.t.ppf(0.95, 1) * sd * np.sqrt(1 + 1 / 2))


def test_auc_matches_definition():
    assert p5.auc([3, 4], [1, 2]) == 1.0
    assert p5.auc([1, 2], [1, 2]) == 0.5
