"""Phase 3b helpers: tumour spikes, grouped folds, fold-internal markers, models and metrics.

Every split groups by source: a healthy background person (a pool's rep1/rep2 pair counts as
one group) or a TCGA tumour case. Nothing here sees a test source while it is fitting.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np


def spike(bg: np.ndarray, tumour: np.ndarray, purity: float, f: float) -> np.ndarray:
    """Beta-level tumour spike at tumour-cell DNA fraction f.

    The tumour sample is purity * tumour cells + (1 - purity) * normal cells, so adding a
    share f / purity of it adds a tumour-cell fraction f. A CpG where the tumour is NaN keeps
    the background value (no tumour signal there).
    """
    w = f / purity
    if not 0 <= w <= 1:
        raise ValueError(f"f / purity = {w:.3f} is outside [0, 1]")
    d = np.where(np.isnan(tumour), 0.0, tumour - bg)
    return bg + w * d


def paired_folds(bg_groups: list[list[str]], cases: list[str],
                 rng: np.random.Generator) -> list[dict]:
    """Outer folds: fold k tests background group k on tumour-case group k.

    Cases are shuffled and split into len(bg_groups) near-equal groups. Training uses the
    other background groups and the other case groups, so no person and no tumour is in
    both train and test.
    """
    k = len(bg_groups)
    perm = rng.permutation(np.asarray(cases, dtype=object))
    case_groups = [list(g) for g in np.array_split(perm, k)]
    folds = []
    for i in range(k):
        folds.append({
            "test_bg": list(bg_groups[i]),
            "train_bg": [b for j, g in enumerate(bg_groups) if j != i for b in g],
            "test_cases": case_groups[i],
            "train_cases": [c for j, g in enumerate(case_groups) if j != i for c in g],
        })
    return folds


def select_tumour_markers(tumours: np.ndarray, normals: np.ndarray, blood: np.ndarray,
                          k: int, blood_max: float = 0.2) -> np.ndarray:
    """Feature set (b): CpGs hypermethylated in tumour vs normal colon and unmethylated in blood.

    tumours: CpG x training-tumour betas (NaN allowed); normals: CpG x normal betas;
    blood: CpG x blood cell-type means. Score = mean tumour - max(mean normal, max blood),
    over CpGs whose max blood beta is <= blood_max. Returns the indices of the top k.
    """
    t = np.nanmean(tumours, axis=1)
    n = np.nanmean(normals, axis=1)
    b = np.nanmax(blood, axis=1)
    score = t - np.maximum(n, b)
    score[~(b <= blood_max) | np.isnan(score)] = -np.inf
    order = np.argsort(-score, kind="stable")
    return order[:k]


def fold_features(features: str, fold: dict, tumours: dict, normals: np.ndarray,
                  blood: np.ndarray, atlas_idx: np.ndarray, k: int) -> np.ndarray:
    """Feature columns for one outer fold.

    "atlas": the fixed atlas CpGs. "tumour": the top-k tumour markers from the fold's training
    cases only (tumours maps case -> beta vector); the test cases' tumours are never read.
    """
    if features == "atlas":
        return atlas_idx
    tr = np.stack([tumours[c] for c in fold["train_cases"]], axis=1)
    return select_tumour_markers(tr, normals, blood, k)


def training_triples(fold: dict, levels: list[float]) -> list[tuple]:
    """Candidate training positives of one fold: training person x training case x f."""
    return [(b, c, f) for b in fold["train_bg"] for c in fold["train_cases"] for f in levels]


def inner_masks(meta, inner: list[dict]) -> list[tuple[np.ndarray, np.ndarray]]:
    """Row masks (train, test) of the training pool for each inner fold.

    meta: one row per training sample with columns person, case ("" for nulls) and f
    (0 for nulls). Inner test rows: the inner test persons' nulls and their spikes of inner
    test cases. Inner train rows: every other person's nulls and spikes of inner training
    cases. A spike that pairs a train person with a test case (or the reverse) is in neither.
    """
    out = []
    for jf in inner:
        pb, pc = meta.person.isin(jf["test_bg"]), meta.case.isin(jf["test_cases"])
        te = pb & ((meta.f == 0) | pc)
        trm = ~pb & ((meta.f == 0) | meta.case.isin(jf["train_cases"]))
        out.append((trm.to_numpy(), te.to_numpy()))
    return out


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Area under the ROC curve (Mann-Whitney, ties count 1/2)."""
    pos, neg = np.asarray(pos, dtype=float), np.asarray(neg, dtype=float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv))
    sv = allv[order]
    i = 0
    while i < len(sv):  # average ranks over ties
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    rp = ranks[:len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def top_k_jaccard(sets: list) -> float:
    """Mean pairwise Jaccard index of feature sets (NaN if fewer than two)."""
    sets = [set(s) for s in sets]
    vals = [len(a & b) / len(a | b) for a, b in combinations(sets, 2) if a | b]
    return float(np.mean(vals)) if vals else float("nan")


def reliability(prob: np.ndarray, y: np.ndarray, bins: int = 10) -> tuple[np.ndarray, ...]:
    """Reliability curve: per probability bin, mean predicted, observed rate and count."""
    prob, y = np.asarray(prob, dtype=float), np.asarray(y, dtype=float)
    edges = np.linspace(0, 1, bins + 1)
    idx = np.clip(np.digitize(prob, edges[1:-1]), 0, bins - 1)
    mp = np.array([prob[idx == b].mean() if (idx == b).any() else np.nan for b in range(bins)])
    ob = np.array([y[idx == b].mean() if (idx == b).any() else np.nan for b in range(bins)])
    n = np.array([(idx == b).sum() for b in range(bins)])
    return mp, ob, n


def weighted_auc(pos: np.ndarray, wp: np.ndarray, neg: np.ndarray, wn: np.ndarray) -> float:
    """AUC with case weights (ties count 1/2), as used by the group bootstrap.

    Equals auc(pos, neg) when all weights are 1. Weight w means the item was drawn w times.
    """
    pos, neg = np.asarray(pos, dtype=float), np.asarray(neg, dtype=float)
    wp, wn = np.asarray(wp, dtype=float), np.asarray(wn, dtype=float)
    tp, tn = wp.sum(), wn.sum()
    if tp == 0 or tn == 0:
        return float("nan")
    o = np.argsort(neg, kind="mergesort")
    sn, cw = neg[o], np.concatenate([[0.0], np.cumsum(wn[o])])
    lo = np.searchsorted(sn, pos, side="left")
    hi = np.searchsorted(sn, pos, side="right")
    below = cw[lo]
    ties = cw[hi] - cw[lo]
    return float((wp * (below + 0.5 * ties)).sum() / (tp * tn))
