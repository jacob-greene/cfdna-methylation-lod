import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from src import classify as cl


def test_spike_scales_by_purity_and_keeps_background_at_nan():
    bg = np.array([0.1, 0.2, 0.3])
    tumour = np.array([0.9, np.nan, 0.3])
    out = cl.spike(bg, tumour, purity=0.5, f=0.1)  # share 0.2 of the tumour sample
    assert np.allclose(out, [0.1 + 0.2 * 0.8, 0.2, 0.3])
    assert np.allclose(cl.spike(bg, tumour, 0.5, 0.0), bg)
    with pytest.raises(ValueError):
        cl.spike(bg, tumour, purity=0.6, f=0.7)


def test_paired_folds_never_share_a_source():
    groups = [["a1", "a2"], ["b"], ["c"], ["d1", "d2"]]
    cases = [f"case{i}" for i in range(23)]
    folds = cl.paired_folds(groups, cases, np.random.default_rng(0))
    assert len(folds) == 4
    seen_cases, seen_bg = [], []
    for fd in folds:
        assert not set(fd["test_bg"]) & set(fd["train_bg"])
        assert not set(fd["test_cases"]) & set(fd["train_cases"])
        # a pool pair is never split between train and test
        for g in groups:
            assert set(g) <= set(fd["test_bg"]) or not set(g) & set(fd["test_bg"])
        seen_cases += fd["test_cases"]
        seen_bg += fd["test_bg"]
    assert sorted(seen_cases) == sorted(cases)
    assert sorted(seen_bg) == sorted(b for g in groups for b in g)


def test_tumour_markers_use_only_the_given_tumours_and_skip_blood_methylated():
    rng = np.random.default_rng(1)
    n = 200
    normals = rng.uniform(0, 0.1, (n, 5))
    blood = rng.uniform(0, 0.1, (n, 3))
    train = rng.uniform(0, 0.1, (n, 8))
    train[:10] = 0.9  # hyper in training tumours
    blood[5:10] = 0.8  # but 5 of them are methylated in blood
    train[100:110] = 0.0  # hyper only in the (unseen) test tumours: must not be picked
    idx = cl.select_tumour_markers(train, normals, blood, k=5)
    assert set(idx) == set(range(5))
    assert not set(cl.select_tumour_markers(train, normals, blood, k=50)) & set(range(100, 110))
    test = np.zeros((n, 4))
    test[100:110] = 0.95
    # positive control: the same rows are picked once the test tumours are (wrongly) included
    leaked = cl.select_tumour_markers(np.c_[train, test], normals, blood, k=15)
    assert set(range(100, 110)) <= set(leaked)


def test_auc_matches_sklearn_with_ties():
    rng = np.random.default_rng(2)
    pos = rng.integers(0, 5, 40).astype(float)
    neg = rng.integers(0, 4, 30).astype(float)
    y = np.r_[np.ones(40), np.zeros(30)]
    assert np.isclose(cl.auc(pos, neg), roc_auc_score(y, np.r_[pos, neg]))


def test_top_k_jaccard():
    assert cl.top_k_jaccard([[1, 2], [1, 2]]) == 1.0
    assert np.isclose(cl.top_k_jaccard([[1, 2], [2, 3]]), 1 / 3)


def test_weighted_auc_equals_auc_at_unit_weights_and_duplication():
    rng = np.random.default_rng(3)
    pos = rng.integers(0, 6, 25).astype(float)
    neg = rng.integers(0, 5, 20).astype(float)
    assert np.isclose(cl.weighted_auc(pos, np.ones(25), neg, np.ones(20)), cl.auc(pos, neg))
    wp = rng.integers(0, 3, 25)
    wn = rng.integers(0, 3, 20)
    assert np.isclose(cl.weighted_auc(pos, wp, neg, wn),
                      cl.auc(np.repeat(pos, wp), np.repeat(neg, wn)))
