"""The vectorized group-bootstrap AUCs of scripts/phase3b_summary.py match a direct loop."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import phase3b_summary as ps

from src import classify as cl


def _toy(seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for k, g in enumerate(["A", "B", "C", "D"]):
        shift = rng.normal()  # each fold model has its own scale offset
        for i in range(k + 1):  # unequal null counts per group
            rows.append((k, g, "neg", "", shift + rng.normal()))
        for c in rng.choice(list("pqrstuvw"), 5, replace=False):
            rows.append((k, g, "pos", c, shift + 1 + rng.normal()))
    d = pd.DataFrame(rows, columns=["fold", "group", "kind", "source", "score"])
    return d[d.kind == "pos"], d[d.kind == "neg"]


def test_boot_auc_matches_loop():
    pos, neg = _toy()
    groups, cases = sorted(neg.group.unique()), sorted(pos.source.unique())
    fold_group = neg.groupby("fold").group.first()
    wg_m, wc_m = ps.boot_weights(np.random.default_rng(1), groups, cases, 20)
    point, pooled, point_w, within = ps.boot_auc(pos, neg, wg_m, wc_m, groups, cases,
                                                 fold_group)
    assert np.isclose(point, cl.auc(pos.score, neg.score))
    assert np.isclose(point_w, np.mean([cl.auc(pos[pos.fold == k].score,
                                               neg[neg.fold == k].score) for k in range(4)]))
    for b in range(20):
        wg = pd.Series(wg_m[b], index=groups)
        wc = pd.Series(wc_m[b], index=cases)
        ref = cl.weighted_auc(pos.score, pos.group.map(wg) * pos.source.map(wc), neg.score,
                              neg.group.map(wg))
        assert np.isclose(pooled[b], ref, equal_nan=True)
        num = den = 0.0
        for k, g in fold_group.items():
            pk, nk = pos[pos.fold == k], neg[neg.fold == k]
            a = cl.weighted_auc(pk.score, pk.source.map(wc), nk.score, np.ones(len(nk)))
            if wg[g] > 0 and not np.isnan(a):
                num, den = num + wg[g] * a, den + wg[g]
        assert np.isclose(within[b], num / den if den else np.nan, equal_nan=True)
