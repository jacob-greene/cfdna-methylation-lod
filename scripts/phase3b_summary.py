"""Phase 3b summary: sensitivity, FPR, LoD, AUC with group-bootstrap CIs, calibration, overlap.

Reads results/phase3b/raw/scores_*.tsv.gz (scripts/phase3b_classify.py) and writes to
results/phase3b/:
  sensitivity.tsv   detection rate per (track, features, HI1, model, f), with a 95 % CI
  fpr.tsv           false-positive rate of held-out blanks, overall and for the worst person
  lod.tsv           LoD per model (f where sensitivity crosses 95 %, log-interpolated), with
                    the FPR, a flag when FPR > 5 % (then the LoD is not valid under EP17), the
                    share of group-bootstrap draws with FPR > 5 %, and (array track) the exact
                    one-sided 95 % upper bound of the FPR over the null people
  auc.tsv           AUC per f, and the paired difference XGB - LR_L2, with 95 % CIs. Two
                    kinds: "pooled" ranks the scores of all outer folds together; "within"
                    is the mean of the per-fold AUCs, so no AUC mixes two models' scales
  b2b_control.tsv   share of blank-into-blank spikes called tumour (no tumour present)
  calibration.tsv   Brier score and reliability bins of the trained models
  overlap.tsv       mean pairwise Jaccard index of the top-50 features across outer folds
  external.tsv      detection rate on the held-out Moss colon plasmas per f
  normal_control.tsv  share of TCGA normal-colon spikes called tumour (batch control)
  figures/phase3b_classifiers.png

CIs come from a two-way group bootstrap: background groups (a pool pair is one group) and
tumour cases are resampled with replacement, the same draw for every model (paired).
LLR rows use the alternative f_alt = f for positives and f_alt = f for the nulls.

Usage: uv run python scripts/phase3b_summary.py [--raw results/phase3b/raw] [--boot 500]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import classify as cl
from src import simulate as sim

SEED = 20261008
KEYS = ["track", "features", "hi1", "model"]
PAIR = {"GSM3455776": "P1", "GSM3455777": "P1", "GSM3455778": "P2", "GSM3455779": "P2",
        "GSM3455780": "P3", "GSM3455781": "P3", "GSM3455782": "P4", "GSM3455783": "P4"}
LABEL = {"GSM3455807": "HI1", "GSM3455840": "HI2", "GSM3455846": "HI3",
         "GSM3455848": "HI4"}
TRACKS = ["array", "30", "100", "300"]
MODELS = ["NNLS", "LLR", "LR_L2", "LR_L1", "XGB"]
COLORS = {"NNLS": "#7f7f7f", "LLR": "#bcbd22", "LR_L2": "#1f78b4", "LR_L1": "#a6cee3",
          "XGB": "#e31a1c"}


def per_f_view(d: pd.DataFrame, f: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Positives at f and the matching nulls (LLR: the f_alt = f statistic)."""
    if d.model.iloc[0] == "LLR":
        d = d[d.f_alt == f]
    return d[(d.kind == "pos") & (d.f == f)], d[d.kind == "neg"]


def boot_weights(rng, groups, cases, n):
    """n group-bootstrap draws as count matrices: (n, groups) and (n, cases)."""
    wg = np.zeros((n, len(groups)))
    wc = np.zeros((n, len(cases)))
    for i in range(n):
        wg[i] = np.bincount(rng.choice(len(groups), len(groups)), minlength=len(groups))
        wc[i] = np.bincount(rng.choice(len(cases), len(cases)), minlength=len(cases))
    return wg, wc


def codes(values: pd.Series, levels: list) -> np.ndarray:
    return pd.Categorical(values, categories=levels).codes


def below(p, n):
    """Count of n below each p, ties counted 1/2 (the AUC kernel)."""
    sn = np.sort(n)
    return (np.searchsorted(sn, p, "left") + np.searchsorted(sn, p, "right")) / 2



def boot_auc(pos, nul, wg_m, wc_m, groups, cases, fold_group):
    """Point and bootstrap AUCs: pooled over folds, and the mean within-fold AUC.

    wg_m, wc_m: bootstrap count matrices (boot_weights). fold_group: outer fold -> its test
    group. Pooled: positives weighted by group x case draws, nulls by group draws.
    Within: per fold, positives weighted by case draws against that fold's nulls;
    folds averaged with their group's draw weight. Vectorized over draws.
    """
    ps, ns = pos.score.to_numpy(float), nul.score.to_numpy(float)
    pg, pc = codes(pos.group, groups), codes(pos.source, cases)
    ng = codes(nul.group, groups)
    cnt = np.stack([below(ps, ns[ng == g]) for g in range(len(groups))], axis=1)
    wp = wg_m[:, pg] * wc_m[:, pc]
    wn = wg_m[:, ng].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = (wp * (wg_m @ cnt.T)).sum(axis=1) / (wp.sum(axis=1) * wn)
    point_p = cnt.sum(axis=1).sum() / (len(ps) * len(ns))
    num_w, den_w, pt = np.zeros(len(wg_m)), np.zeros(len(wg_m)), []
    for k, g in fold_group.items():
        pk, nk = (pos.fold == k).to_numpy(), nul.fold == k
        if not pk.any() or not nk.any():
            continue
        h = below(ps[pk], ns[nk.to_numpy()]) / nk.sum()
        wck = wc_m[:, pc[pk]]
        tot = wck.sum(axis=1)
        gi = groups.index(g)
        with np.errstate(invalid="ignore", divide="ignore"):
            ak = np.where(tot > 0, wck @ h / tot, np.nan)
        ok = ~np.isnan(ak)
        num_w[ok] += wg_m[ok, gi] * ak[ok]
        den_w[ok] += wg_m[ok, gi]
        pt.append(h.mean())
    with np.errstate(invalid="ignore", divide="ignore"):
        within = num_w / den_w
    return point_p, pooled, float(np.mean(pt)) if pt else np.nan, within


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("results/phase3b/raw"))
    ap.add_argument("--out", type=Path, default=Path("results/phase3b"))
    ap.add_argument("--fig", type=Path, default=Path("figures/phase3b_classifiers.png"))
    ap.add_argument("--boot", type=int, default=500)
    args = ap.parse_args()
    files = sorted(p for p in args.raw.glob("scores_*.tsv.gz") if "_f" not in p.stem)
    raw = pd.concat([pd.read_csv(p, sep="\t", dtype={"track": str}) for p in files],
                    ignore_index=True)
    raw["detected"] = raw.score > raw.threshold
    raw["group"] = raw.person.map(lambda p: PAIR.get(p, LABEL.get(p, p)))
    rng = np.random.default_rng(SEED)

    sens, fpr, lods, aucs, cal, ext, normal, b2b = [], [], [], [], [], [], [], []
    for key, d in raw.groupby(KEYS, sort=False):
        neg = d[d.kind == "neg"]
        # FPR: for LLR, one threshold per f_alt; report the range.
        by = neg.groupby("f_alt", dropna=False) if key[3] == "LLR" else [(np.nan, neg)]
        fr = [(fa, len(g), g.detected.mean(), g.groupby("group").detected.mean().idxmax(),
               g.groupby("group").detected.mean().max()) for fa, g in by]
        fpr += [(*key, *r) for r in fr]
        fmin, fmax = min(r[2] for r in fr), max(r[2] for r in fr)

        # Group bootstrap draws shared by every f of this configuration.
        groups = sorted(d.group.unique())
        cases = sorted(d.loc[d.kind == "pos", "source"].unique())
        WG, WC = boot_weights(rng, groups, cases, args.boot)
        levels = sorted(d.loc[d.kind == "pos", "f"].unique())
        rates = []
        for f in levels:
            pos, nul = per_f_view(d, f)
            r = pos.detected.mean()
            rates.append(r)
            w = WG[:, codes(pos.group, groups)] * WC[:, codes(pos.source, cases)]
            tot = w.sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                bs = np.where(tot > 0, w @ pos.detected.to_numpy(float) / tot, np.nan)
            sens.append((*key, f, len(pos), r, np.nanquantile(bs, 0.025),
                         np.nanquantile(bs, 0.975)))
        lod_v = sim.lod_interp(levels, rates)
        # FPR uncertainty: the unit is the person. Group bootstrap of the (worst f_alt) FPR;
        # on the array track (one null result per person) also the exact one-sided bound.
        nb = [(g.groupby("group").detected.agg(["sum", "count"])) for _, g in by]
        bfpr = np.max([(WG[:, codes(t.index, groups)] @ t["sum"].to_numpy(float))
                       / (WG[:, codes(t.index, groups)] @ t["count"].to_numpy(float))
                       for t in nb], axis=0)
        n_people = neg.person.nunique()
        if key[0] == "array":
            x = int(neg.groupby("f_alt", dropna=False).detected.sum().max())
            upper = 1.0 if x == n_people else float(stats.beta.ppf(0.95, x + 1, n_people - x))
        else:
            upper = np.nan
        lods.append((*key, lod_v, sim.lod_from_rates(levels, rates), fmin, fmax,
                     fmax <= 0.05, d.group.nunique(), n_people, len(cases),
                     float((bfpr > 0.05).mean()), upper))

        # External Moss plasmas.
        e = d[d.kind == "ext"]
        if key[3] == "LLR":
            e = e[e.f_alt == e.f]
        for (src, f), g in e.groupby(["source", "f"]):
            ext.append((*key, src, f, len(g), g.detected.mean()))

        # Batch control: TCGA normal-colon spikes called tumour.
        nr = d[d.kind == "nrm"]
        if key[3] == "LLR":
            nr = nr[nr.f_alt == nr.f]
        for f, g in nr.groupby("f"):
            normal.append((*key, f, len(g), g.detected.mean()))

        # Blank-into-blank control: another person's blank spiked into a test blank.
        bb = d[d.kind == "b2b"]
        if key[3] == "LLR":
            bb = bb[bb.f_alt == bb.f]
        for f, g in bb.groupby("f"):
            b2b.append((*key, f, len(g), g.detected.mean(),
                        g.groupby("source").detected.mean().idxmax(),
                        g.groupby("source").detected.mean().max()))

        # Calibration of the trained models (test positives and nulls only).
        if key[3] not in ("NNLS", "LLR") and d.prob.notna().any():
            t = d[d.kind.isin(["pos", "neg"])]
            yv = (t.kind == "pos").to_numpy(dtype=float)
            p = t.prob.to_numpy()
            mp, ob, n = cl.reliability(p, yv)
            cal.append((*key, "brier", float(np.mean((p - yv) ** 2)), yv.mean(), len(yv)))
            for b in range(len(mp)):
                cal.append((*key, f"bin{b}", mp[b], ob[b], n[b]))

    sens = pd.DataFrame(sens, columns=[*KEYS, "f", "n", "sensitivity", "ci_lo", "ci_hi"])
    fpr = pd.DataFrame(fpr, columns=[*KEYS, "f_alt", "n", "fpr", "worst_group",
                                     "worst_group_fpr"])
    lods = pd.DataFrame(lods, columns=[*KEYS, "lod_interp", "lod_grid", "fpr_min", "fpr_max",
                                       "specificity_ok", "n_bg_groups", "n_null_people",
                                       "n_tumour_cases", "boot_share_fpr_gt5",
                                       "fpr_upper95_exact"])
    b2b = pd.DataFrame(b2b, columns=[*KEYS, "f", "n", "called_tumour", "worst_donor",
                                     "worst_donor_called"])
    ext = pd.DataFrame(ext, columns=[*KEYS, "plasma", "f", "n", "detection_rate"])
    cal = pd.DataFrame(cal, columns=[*KEYS, "item", "a", "b", "n"])
    normal = pd.DataFrame(normal, columns=[*KEYS, "f", "n", "called_tumour"])

    # AUC per f and paired difference XGB - LR_L2 (same bootstrap draws for both models).
    for (track, feats, hi1), d in raw.groupby(["track", "features", "hi1"], sort=False):
        groups = sorted(d.group.unique())
        cases = sorted(d.loc[d.kind == "pos", "source"].unique())
        WG, WC = boot_weights(rng, groups, cases, args.boot)
        # Outer fold k tests one background group, so a group weight is a fold weight.
        fold_group = d[d.kind == "neg"].groupby("fold").group.first()

        for f in sorted(d.loc[d.kind == "pos", "f"].unique()):
            boot, point = {}, {}
            for m, dm in d.groupby("model"):
                pos, nul = per_f_view(dm, f)
                a, bs, aw, bw = boot_auc(pos, nul, WG, WC, groups, cases, fold_group)
                boot[m] = {"pooled": np.array(bs), "within": np.array(bw)}
                point[m] = {"pooled": a, "within": aw}
                aucs.append((track, feats, hi1, m, "pooled", f, a, np.nanquantile(bs, 0.025),
                             np.nanquantile(bs, 0.975)))
                aucs.append((track, feats, hi1, m, "within", f, aw, np.nanquantile(bw, 0.025),
                             np.nanquantile(bw, 0.975)))
            if "XGB" in boot and "LR_L2" in boot:
                for kind in ("pooled", "within"):
                    dd = boot["XGB"][kind] - boot["LR_L2"][kind]
                    diff = point["XGB"][kind] - point["LR_L2"][kind]
                    aucs.append((track, feats, hi1, "XGB-LR_L2", kind, f, diff,
                                 np.nanquantile(dd, 0.025), np.nanquantile(dd, 0.975)))
    aucs = pd.DataFrame(aucs, columns=[*KEYS, "auc_kind", "f", "auc", "ci_lo", "ci_hi"])

    top = pd.concat([pd.read_csv(p, sep="\t", dtype={"track": str})
                     for p in sorted(args.raw.glob("top_*.tsv")) if "_f" not in p.stem])
    ov = [(*k, g.fold.nunique(), cl.top_k_jaccard([s.cpg for _, s in g.groupby("fold")]))
          for k, g in top.groupby(KEYS)]
    ov = pd.DataFrame(ov, columns=[*KEYS, "n_folds", "top50_jaccard"])

    args.out.mkdir(parents=True, exist_ok=True)
    for name, df in [("sensitivity", sens), ("fpr", fpr), ("lod", lods), ("auc", aucs),
                     ("calibration", cal), ("overlap", ov), ("external", ext),
                     ("normal_control", normal), ("b2b_control", b2b)]:
        df.to_csv(args.out / f"{name}.tsv", sep="\t", index=False)

    # Figure: sensitivity vs f, one line per model, HI1 excluded; rows = feature sets.
    s = sens[sens.hi1 == "without"]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
    for i, feats in enumerate(["atlas", "tumour"]):
        for j, track in enumerate(TRACKS):
            ax = axes[i, j]
            for m in MODELS:
                g = s[(s.features == feats) & (s.track == track) & (s.model == m)]
                if g.empty:
                    continue
                ok = lods[(lods.features == feats) & (lods.track == track) & (lods.model == m)
                          & (lods.hi1 == "without")].specificity_ok.all()
                ax.plot(g.f * 100, g.sensitivity, "-o" if ok else "--o", color=COLORS[m],
                        mfc=COLORS[m] if ok else "white", ms=4, label=m)
                ax.fill_between(g.f * 100, g.ci_lo, g.ci_hi, color=COLORS[m], alpha=0.15)
            ax.axhline(0.95, color="k", lw=0.6, ls=":")
            ax.set_xscale("log")
            ttl = "array betas" if track == "array" else f"reads, D = {track}"
            ax.set_title(f"{feats} features, {ttl}", fontsize=10)
            if i == 1:
                ax.set_xlabel("tumour-cell DNA fraction f (%)")
            if j == 0:
                ax.set_ylabel("sensitivity at the valid-rule threshold")
    handles = {}
    for ax in axes.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            handles.setdefault(lab, h)
    fig.legend(handles.values(), handles.keys(), loc="lower center", ncol=len(handles),
               fontsize=9)
    fig.suptitle("Phase 3b: TCGA-COAD tumour spikes into held-out healthy blanks (HI1 "
                 "excluded). Dashed, open markers: held-out FPR > 5 %.", fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    args.fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.fig, dpi=150)
    print(lods.to_string())


if __name__ == "__main__":
    main()
