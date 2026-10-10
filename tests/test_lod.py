import numpy as np
import pandas as pd
import pytest

from src import lod


def test_lob_rank_interpolates():
    v = np.arange(1, 101, dtype=float)  # rank 0.5 + 95 = 95.5 -> halfway between 95 and 96
    assert lod.lob_nonparametric(v) == pytest.approx(95.5)


def test_lob_small_n_is_max():
    assert lod.lob_nonparametric([0.3, 0.1, 0.2]) == 0.3


def test_lob_uses_only_given_blanks():
    blanks = np.array([0.01, 0.02, 0.015])
    spiked = np.array([0.5, 0.6])
    assert lod.lob_nonparametric(blanks) == lod.lob_nonparametric(blanks.copy())
    assert lod.lob_nonparametric(blanks) < spiked.min()


def test_lod_parametric_matches_hand_calc():
    g = [np.array([1.0, 3.0]), np.array([2.0, 4.0])]
    out = lod.lod_parametric(0.5, g)
    sd = np.sqrt((2 + 2) / (4 - 2))
    cp = 1.645 / (1 - 1 / (4 * 2))
    assert out["sd_pooled"] == pytest.approx(sd)
    assert out["lod"] == pytest.approx(0.5 + cp * sd)


def test_empirical_lod_requires_all_higher_levels():
    assert lod.empirical_lod([0, 1, 2, 3], [0, 1.0, 0.5, 1.0]) == 3
    assert np.isnan(lod.empirical_lod([0, 1], [0, 0.5]))


def test_nnls_recovers_known_mix():
    rng = np.random.default_rng(0)
    a = rng.uniform(0, 1, (300, 4))
    w = np.array([0.6, 0.3, 0.1, 0.0])
    y = a @ w
    y[5] = np.nan
    assert lod.nnls_fractions(y, a) == pytest.approx(w, abs=1e-8)


def test_mix_is_linear():
    bg, sp = np.array([0.2, np.nan]), np.array([0.8, 0.5])
    out = lod.mix(bg, sp, 0.25)
    assert out[0] == pytest.approx(0.35) and np.isnan(out[1])


def test_lob_by_spike_uses_only_blank_rows_of_the_lob_set():
    rng = np.random.default_rng(4)
    rows = []
    for spike in ("s1", "s2"):
        for bg in ("b1", "b2", "b3", "e1", "e2"):
            for f in (0.0, 0.01, 0.03):
                # Spiked rows and non-LoB backgrounds sit far above the LoB-set blanks.
                high = f > 0 or bg.startswith("e")
                rows.append({"spike": spike, "background": bg, "f": f,
                             "colon": rng.uniform(0.5, 0.9) if high
                             else rng.uniform(0, 0.02)})
    d = pd.DataFrame(rows)
    got = lod.lob_by_spike(d, ["b1", "b2", "b3"])
    for spike in ("s1", "s2"):
        blanks = d[(d.spike == spike) & (d.f == 0) & d.background.isin(["b1", "b2", "b3"])]
        assert got[spike] == lod.lob_nonparametric(blanks.colon)
        assert got[spike] < 0.02
