import numpy as np

from src import simulate as sim


def test_counts_reproducible_under_seed():
    bg = np.linspace(0.05, 0.95, 50)
    t = 1 - bg
    n = np.full((4, 50), 30)
    a = sim.simulate_counts(bg, t, 0.1, n, np.random.default_rng(7))
    b = sim.simulate_counts(bg, t, 0.1, n, np.random.default_rng(7))
    assert np.array_equal(a, b)


def test_demix_inverts_mixing():
    rng = np.random.default_rng(0)
    tumour, bg = rng.uniform(0, 1, 100), rng.uniform(0, 1, 100)
    plasma = 0.6 * tumour + 0.4 * bg
    assert np.allclose(sim.demix(plasma, bg, 0.6), tumour)


def test_llr_zero_under_identical_hypotheses_and_positive_for_tumour():
    p0 = np.full(200, 0.1)
    p1 = 0.2 * np.ones(200) + 0.8 * p0  # tumour fully methylated, f = 0.2
    n = np.full((500, 200), 100)
    rng = np.random.default_rng(1)
    k0 = sim.simulate_counts(p0, None, 0.0, n, rng)
    k1 = sim.simulate_counts(p0, np.ones(200), 0.2, n, rng)
    assert np.allclose(sim.llr(k0, n, p0, p0), 0)
    assert sim.llr(k1, n, p0, p1).mean() > 0 > sim.llr(k0, n, p0, p1).mean()


def test_input_cap_limits_molecules():
    m = sim.molecules(10_000, 300.0, 1000, 5, np.random.default_rng(2))
    assert m.max() < 1000 and abs(m.mean() - 300) < 5
    assert (sim.molecules(100, None, 10, 2, np.random.default_rng(2)) == 100).all()


def test_rank_markers_prefers_stable_contrast():
    train = np.array([[0.1, 0.1, 0.5], [0.1, 0.5, 0.5], [0.1, 0.9, 0.5]])
    tumour = np.array([0.9, 0.9, 0.6])
    assert sim.rank_markers(tumour, train)[0] == 0


def test_lod_helpers():
    lv, r = [1e-3, 1e-2, 1e-1], [0.1, 0.9, 1.0]
    assert sim.lod_from_rates(lv, r) == 1e-1
    x = sim.lod_interp(lv, r)
    assert 1e-2 < x < 1e-1
    assert np.isnan(sim.lod_interp(lv, [0.1, 0.2, 0.5]))


def test_prediction_bound_covers_new_person_at_nominal_rate():
    rng = np.random.default_rng(3)
    hits = []
    for _ in range(2000):
        people = rng.normal(0, 1, 11)
        per = [m + rng.normal(0, 0.1, 20) for m in people]
        new = rng.normal(0, 1) + rng.normal(0, 0.1)
        hits.append(new > sim.prediction_bound(per))
    assert 0.03 < np.mean(hits) < 0.07
