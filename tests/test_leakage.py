"""Leakage tests for the Phase 3b fold design (scripts/phase3b_classify.py).

The script builds every fold with the helpers tested here, so a leak in them is a leak in
the reported results.
"""

import numpy as np
import pandas as pd

from src import classify as cl

GROUPS = [["p1a", "p1b"], ["p2a", "p2b"], ["h2"], ["h3"], ["h4"]]
CASES = [f"case{i}" for i in range(31)]
F_TRAIN = [1e-3, 1e-2, 1e-1]


class Recording(dict):
    """A dict that records every key read through [] (the way the script reads tumours)."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.read = set()

    def __getitem__(self, key):
        self.read.add(key)
        return super().__getitem__(key)


def _training_meta(fold):
    """The training pool of one outer fold, as the script lays it out: nulls, then spikes."""
    rows = [(b, "", 0.0) for b in fold["train_bg"] for _ in range(2)]
    rows += cl.training_triples(fold, F_TRAIN)
    return pd.DataFrame(rows, columns=["person", "case", "f"])


def test_no_source_in_both_train_and_test_folds():
    folds = cl.paired_folds(GROUPS, CASES, np.random.default_rng(0))
    for fd in folds:
        meta = _training_meta(fd)
        # Outer level: the training pool holds no test person and no test tumour case.
        assert not set(meta.person) & set(fd["test_bg"])
        assert not set(meta.case) & set(fd["test_cases"])
        # Inner level: same rule inside the training part.
        groups = [g for g in GROUPS if not set(g) & set(fd["test_bg"])]
        inner = cl.paired_folds(groups, fd["train_cases"], np.random.default_rng(1))
        masks = cl.inner_masks(meta, inner)
        covered = np.zeros(len(meta), dtype=bool)
        for jf, (trm, te) in zip(inner, masks):
            assert not (trm & te).any()
            tr_rows, te_rows = meta[trm], meta[te]
            assert not set(tr_rows.person) & set(te_rows.person)
            assert not (set(tr_rows.case) - {""}) & (set(te_rows.case) - {""})
            assert set(te_rows.person) <= set(jf["test_bg"])
            assert set(te_rows.case) - {""} <= set(jf["test_cases"])
            covered |= te
        # Every training null is scored out of fold exactly once (it sets the threshold).
        assert covered[meta.f == 0].all()


def test_fold_feature_selection_never_sees_test_samples():
    rng = np.random.default_rng(2)
    n = 300
    normals = rng.uniform(0, 0.1, (n, 4))
    blood = rng.uniform(0, 0.1, (n, 3))
    tumours = {c: rng.uniform(0, 0.1, n) for c in CASES}
    for v in tumours.values():
        v[:40] = 0.0  # unmethylated in every tumour, unless a test tumour sets them below
    folds = cl.paired_folds(GROUPS, CASES, np.random.default_rng(0))
    for fd in folds:
        for c in fd["test_cases"]:  # strong markers that only the test tumours carry
            tumours[c] = tumours[c].copy()
            tumours[c][:40] = 1.0
        rec = Recording(tumours)
        feat = cl.fold_features("tumour", fd, rec, normals, blood, None, k=20)
        assert rec.read == set(fd["train_cases"])
        assert not set(feat) & set(range(40))
        # Changing the test tumours must not change the selected features.
        poisoned = Recording({c: (np.full(n, np.nan) if c in fd["test_cases"] else v)
                              for c, v in tumours.items()})
        same = cl.fold_features("tumour", fd, poisoned, normals, blood, None, k=20)
        assert np.array_equal(feat, same)
        for c in fd["test_cases"]:
            tumours[c][:40] = 0.0
    atlas_idx = np.arange(5)
    rec = Recording(tumours)
    assert cl.fold_features("atlas", folds[0], rec, normals, blood, atlas_idx, 20) is atlas_idx
    assert not rec.read
