import numpy as np
import pytest

from src import phase5 as p5
from src import phase6 as p6


def _rec(block, m, u, w=None):
    a = lambda x: np.asarray(x, dtype=np.int64)
    w = np.ones(len(block)) if w is None else w
    return p5.Records(a(range(len(block))), a(block), a(m), a(u), a(w))


def test_summarise_counts_calls_design_and_tumour_like_fragments():
    rec = _rec([0, 0, 1, 1], [3, 1, 0, 2], [0, 1, 4, 0], [1, 2, 1, 1])
    S = p6.summarise(rec, 2, rho=0.5)
    get = lambda f: S[p6.FIELDS.index(f)]
    assert get("M").tolist() == [3 + 2 * 1, 2]
    assert get("N").tolist() == [3 + 2 * 2, 4 + 2]
    # G = w n (1 + (n - 1) rho): 3 * 2 + 2 * 2 * 1.5 = 12; 4 * 2.5 + 2 * 1.5 = 13
    assert get("G").tolist() == [12, 13]
    assert get("cM").tolist() == [1, 0]  # only the 3-call all-methylated record
    assert get("cU").tolist() == [0, 1]
    assert get("n3").tolist() == [1, 1]


def test_tumour_like_follows_direction():
    S = np.zeros((len(p6.FIELDS), 2))
    S[p6.FIELDS.index("cM")] = [5, 7]
    S[p6.FIELDS.index("cU")] = [1, 2]
    assert p6.tumour_like(S, np.array([1.0, -1.0])).tolist() == [5, 2]


def _donor(rng, p, depth, nb):
    """Summary of a donor whose fragments have 4 calls, all ~ Bernoulli(p_b) independently."""
    blocks = np.repeat(np.arange(nb), depth)
    m = rng.binomial(4, p[blocks])
    return p6.summarise(_rec(blocks, m, 4 - m), nb, rho=0.0)


def test_wz_is_near_standard_normal_under_a_binomial_null():
    rng = np.random.default_rng(0)
    nb = 400
    p = rng.uniform(0.05, 0.95, nb)
    T = np.clip(p + rng.choice([-0.4, 0.4], nb), 0, 1)
    donors = [_donor(rng, p, 30, nb) for _ in range(14)]
    model = p6.fit(donors, T)
    assert model.phi_mean < 0.02  # no between-donor variance was simulated
    z = np.array([p6.wz(_donor(rng, p, 30, nb), model, 1) for _ in range(200)])
    assert abs(z.mean()) < 0.3 and 0.7 < z.std() < 1.3


def test_wz_rises_with_a_tumour_shift():
    rng = np.random.default_rng(1)
    nb = 400
    p = rng.uniform(0.05, 0.3, nb)
    T = np.full(nb, 0.9)
    model = p6.fit([_donor(rng, p, 30, nb) for _ in range(14)], T)
    z0 = p6.wz(_donor(rng, p, 30, nb), model, 1)
    z1 = p6.wz(_donor(rng, p + 0.05 * (T - p), 30, nb), model, 1)
    assert z1 > z0 + 5


def test_frag_selects_quiet_blocks_and_counts_excess():
    rng = np.random.default_rng(2)
    nb = 50
    p = np.full(nb, 0.02)
    p[:5] = 0.6  # noisy blocks: healthy fragments often all methylated
    T = np.full(nb, 0.9)
    donors = [_donor(rng, p, 40, nb) for _ in range(14)]
    model = p6.fit(donors, T)
    assert not model.sel[:5].any() and model.sel[5:].all()
    S = _donor(rng, p, 40, nb)
    C0, E0 = p6.frag_counts(S, model)
    S[p6.FIELDS.index("cM"), 10] += 6
    S[p6.FIELDS.index("n3"), 10] += 6
    C1, E1 = p6.frag_counts(S, model)
    assert C1 == C0 + 6 and E1 > E0
    assert p6.frag(S, model, 1) > 2


def test_frag_overdispersion_is_zero_without_donor_variation_and_positive_with_it():
    rng = np.random.default_rng(3)
    nb = 200
    p = np.full(nb, 0.3)  # all-methylated rate 0.3^4 = 0.008: selected
    T = np.full(nb, 0.9)
    flat = p6.fit([_donor(rng, p, 50, nb) for _ in range(14)], T)
    varied = p6.fit([_donor(rng, p * rng.choice([0.7, 1.3]), 50, nb) for _ in range(14)], T)
    assert flat.s2 < 0.05 < varied.s2


def test_eval_stage_guard_refuses_uncommitted_changes(monkeypatch):
    import importlib.util
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("phase6_score",
                                                  root / "scripts/phase6_score.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    class R:
        def __init__(self, code, out=""):
            self.returncode, self.stdout = code, out

    calls = {"diff": 1}
    def fake_run(cmd, **kw):
        if "ls-files" in cmd:
            return R(0)
        if "diff" in cmd:
            return R(calls["diff"])
        return R(0, "abc\n")
    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    with pytest.raises(SystemExit):
        mod.preregistered_sha()
    calls["diff"] = 0
    assert mod.preregistered_sha() == "abc"
