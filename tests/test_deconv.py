import numpy as np
import pandas as pd

from src import deconv

SEED = 20261008


def _toy_atlas(rng, n=400, k=6):
    return pd.DataFrame(rng.uniform(0, 1, (n, k)), index=[f"cg{i:08d}" for i in range(n)],
                        columns=[f"t{j}" for j in range(k)])


def test_nnls_nonnegative_and_sums_to_one():
    rng = np.random.default_rng(SEED)
    atl = _toy_atlas(rng)
    y = pd.Series(rng.uniform(0, 1, len(atl)), index=atl.index)
    f, n, _ = deconv.nnls_one(y, atl)
    assert (f >= 0).all()
    assert abs(f.sum() - 1) < 1e-9
    assert n == len(atl)


def test_recovers_known_proportions():
    rng = np.random.default_rng(SEED)
    atl = _toy_atlas(rng)
    props = np.array([0.5, 0.3, 0.1, 0.07, 0.03, 0.0])
    y = deconv.synthetic_mix(atl, props, noise_sd=0.01, rng=rng)
    f, _, _ = deconv.nnls_one(y, atl)
    assert np.max(np.abs(f.to_numpy() - props)) < 0.02


def test_missing_values_dropped_per_sample():
    rng = np.random.default_rng(SEED)
    atl = _toy_atlas(rng)
    y = pd.Series(atl.to_numpy() @ np.full(6, 1 / 6), index=atl.index)
    y.iloc[:50] = np.nan
    _, n, _ = deconv.nnls_one(y, atl)
    assert n == len(atl) - 50


def test_dedup_keeps_unique_ids():
    atl = pd.DataFrame({"a": [0.1, 0.2, 0.1]}, index=["cg1", "cg2", "cg1"])
    assert list(deconv.dedup_atlas(atl).index) == ["cg1", "cg2"]
