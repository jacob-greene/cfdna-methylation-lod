"""Phase 3b: logistic regression vs XGBoost for tumour detection, with valid-rule thresholds.

One run = one (track, feature set, HI1 setting). It writes every held-out test score;
scripts/phase3b_summary.py turns them into sensitivity, LoD, FPR, AUC and calibration.

Task: tumour present at fraction f vs f = 0 (Phase 3b).

Tumour sources: TCGA-COAD primary tumours with CPE purity >= 0.6 (one file per case), spiked
at a share f / purity so that f is the tumour-cell DNA fraction (src/classify.spike). The
Moss colon-cancer plasmas 806, 684 and 2010 are an external test set: never trained on, and
spiked at a share f / c(plasma), where c is the plasma's NNLS colon fraction.

Backgrounds: the healthy cfDNA blanks. HI1 is excluded in the headline run ("--hi1
without") and kept in the sensitivity run ("--hi1 with"). A pool's rep1/rep2 pair is one
group.

Tracks: "array" = beta-level mixtures; "30", "100", "300" = the Phase 3 binomial read model
at that depth (methylated reads ~ Binomial(D, beta)), features = simulated betas k / D.

Feature sets: "atlas" = the Moss atlas CpGs; "tumour" = 500 CpGs hypermethylated in the
training TCGA tumours vs TCGA normals and unmethylated in Moss blood cell types, selected
inside each training fold (src/classify.select_tumour_markers).

Models: NNLS colon fraction (atlas features only) and summed binomial LLR (read tracks only)
as untrained baselines; logistic regression with L2 and L1 penalties on standardized
features; XGBoost (hist). The three trained models get the same tuning budget: 6
configurations, scored by the pooled out-of-fold AUC over the inner folds for positives at
f = 1 % and 3 %.

Splitting: outer fold k tests background group k on tumour-case group k (src/classify.
paired_folds); the inner folds repeat that design inside the training part. No person and no
tumour case is in both train and test, at either level.

Threshold (valid rule): each training person's f = 0 scores come from the inner-fold model
that did not see that person's group. The threshold is the one-sided 95 % normal prediction
bound over people (src/simulate.prediction_bound). It is applied unchanged to the test fold.

Test score of a trained model (revised after an independent re-analysis): the mean margin
of the inner-fold models of the chosen configuration. These are the models whose out-of-fold
scores set the threshold, so threshold and test score share one scale. A model refit on the
whole training fold has larger margins, which put the threshold too low. Two diagnostic
variants are kept: "<model>_refit" (the refit model with the inner threshold, the first
version's rule) and "<model>_refitcal" (the refit model with the inner null scores mapped onto
its scale by a linear fit over shared training positives, as in that re-analysis).

Controls: (1) batch: the 38 TCGA solid-tissue normals are spiked like tumours (purity 1,
f = 1, 3, 10 %) into the test blanks; a score that calls them tumour reads TCGA batch or
colon tissue, not cancer. (2) blank-into-blank: each test blank is spiked with another
person's blank (not in its group, never HI1) at f = 1, 3, 10 %. The spike holds no tumour,
so a specific model calls it at about its FPR. (3) --permute shuffles the training labels
(first version's leakage control; kept for reference only).

Usage: uv run python scripts/phase3b_classify.py --track 100 --features atlas --hi1 without
       [--folds 0,1] [--max-pos 4000] [--threads 4] [--out results/phase3b/raw]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import classify as cl
from src import deconv
from src import simulate as sim

SEED = 20261008
# sklearn >= 1.8 deprecates penalty= in favour of l1_ratio=; the fits are identical.
warnings.filterwarnings("ignore", module="sklearn")
COLON = "Colon_epithelial_cells"
POOL_PAIRS = [("GSM3455776", "GSM3455777"), ("GSM3455778", "GSM3455779"),
              ("GSM3455780", "GSM3455781"), ("GSM3455782", "GSM3455783")]
HI1, HI2, HI3, HI4 = "GSM3455807", "GSM3455840", "GSM3455846", "GSM3455848"
PLASMAS = {"pl806": "GSM3455808", "pl684": "GSM3455809", "pl2010": "GSM3455842"}
F_TRAIN = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]  # study array f grid without 0
F_EVAL = sorted(F_TRAIN + [2e-2, 5e-2])  # two added levels to place the 95 % crossing
F_TUNE = (1e-2, 3e-2)
F_NORMAL = [1e-2, 3e-2, 1e-1]  # batch control: TCGA normal-colon spikes
F_B2B = F_NORMAL  # blank-into-blank control
N_TUMOUR_MARKERS = 500
NULL_TRAIN, NULL_TEST, EXT_REPS = 40, 200, 20

# Six configurations per model. The L1 grid starts higher: at C <= 1e-3 every weight is 0.
L2_GRID = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0]
L1_GRID = [3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]
XGB_GRID = [{"max_depth": d, "learning_rate": lr} for d in (2, 4, 6) for lr in (0.05, 0.2)]


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def make_model(name: str, cfg, threads: int):
    if name == "LR_L2":
        return LogisticRegression(C=cfg, penalty="l2", solver="lbfgs", max_iter=3000,
                                  class_weight="balanced")
    if name == "LR_L1":
        return LogisticRegression(C=cfg, penalty="l1", solver="liblinear", max_iter=3000,
                                  class_weight="balanced", random_state=SEED)
    return XGBClassifier(n_estimators=300, subsample=0.8, colsample_bytree=0.5,
                         tree_method="hist", n_jobs=threads, random_state=SEED,
                         eval_metric="logloss", **cfg)


class Fitted:
    """A trained model with its scaler; margin() is the score the threshold applies to."""

    def __init__(self, name, cfg, X, y, threads):
        self.name = name
        self.scaler = StandardScaler().fit(X) if name.startswith("LR") else None
        Xs = self.scaler.transform(X) if self.scaler else X
        self.m = make_model(name, cfg, threads)
        if name == "XGB":
            self.m.set_params(scale_pos_weight=float((y == 0).sum() / max((y == 1).sum(), 1)))
        self.m.fit(Xs, y)

    def _x(self, X):
        return self.scaler.transform(X) if self.scaler else X

    def margin(self, X):
        if self.name == "XGB":
            return self.m.predict(self._x(X), output_margin=True)
        return self.m.decision_function(self._x(X))

    def prob(self, X):
        return self.m.predict_proba(self._x(X))[:, 1]

    def top(self, k=50):
        if self.name == "XGB":
            g = self.m.get_booster().get_score(importance_type="gain")
            imp = np.zeros(self.m.n_features_in_)
            for key, v in g.items():
                imp[int(key[1:])] = v
        else:
            imp = np.abs(self.m.coef_[0])
        return np.argsort(-imp, kind="stable")[:k]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", choices=["array", "30", "100", "300"], required=True)
    ap.add_argument("--features", choices=["atlas", "tumour"], required=True)
    ap.add_argument("--hi1", choices=["with", "without"], default="without")
    ap.add_argument("--folds", default=None, help="comma-separated outer folds (test runs)")
    ap.add_argument("--max-pos", type=int, default=4000)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    ap.add_argument("--models", default="LR_L2,LR_L1,XGB")
    ap.add_argument("--out", type=Path, default=Path("results/phase3b/raw"))
    ap.add_argument("--permute", action="store_true", help="shuffle training labels (control)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data = data_dir() / "phase3b"
    depth = None if args.track == "array" else int(args.track)
    trained = args.models.split(",") if args.models else []

    geo = pd.read_parquet(data / "geo.parquet")
    tum = pd.read_parquet(data / "tumours.parquet")
    nor = pd.read_parquet(data / "normals.parquet").to_numpy(dtype=float)
    blood = pd.read_parquet(data / "blood.parquet").to_numpy(dtype=float)
    tinfo = pd.read_csv(data / "tumours.tsv", sep="\t").set_index("id")
    cpgs = geo.index
    a_cpgs = (data / "atlas_cpgs.txt").read_text().split()
    a_idx = cpgs.get_indexer(a_cpgs)
    atl = deconv.dedup_atlas(atlas_mod.read_atlas(data.parent)).loc[a_cpgs]
    a_std = atl.to_numpy(dtype=float)
    ic = list(atl.columns).index(COLON)

    bg_groups = [list(p) for p in POOL_PAIRS] + [[HI2], [HI3], [HI4]]
    if args.hi1 == "with":
        bg_groups.append([HI1])
    group_of = {b: i for i, g in enumerate(bg_groups) for b in g}
    B = {g: geo[g].to_numpy(dtype=float) for g in geo.columns}
    case_ids = tinfo.case.tolist()
    T = {c: tum[fid].to_numpy(dtype=float) for fid, c in zip(tinfo.index, tinfo.case)}
    purity = dict(zip(tinfo.case, tinfo.CPE))
    c_pl = {k: float(sim.nnls_colon(B[g][a_idx][None], a_std, ic)[0])
            for k, g in PLASMAS.items()}
    task = [SEED, ["array", "30", "100", "300"].index(args.track),
            ["atlas", "tumour"].index(args.features), int(args.hi1 == "with")]
    folds = cl.paired_folds(bg_groups, case_ids, np.random.default_rng(task))
    run = [int(x) for x in args.folds.split(",")] if args.folds else range(len(folds))
    print(f"track={args.track} features={args.features} hi1={args.hi1} CpGs={len(cpgs)} "
          f"atlas={len(a_idx)} cases={len(case_ids)} groups={len(bg_groups)} "
          f"c={ {k: round(v, 4) for k, v in c_pl.items()} }", flush=True)

    def observe(p: np.ndarray, rng, reps: int = 1) -> np.ndarray:
        """Array track: the betas themselves. Read tracks: Binomial(D, p) / D."""
        p = np.clip(np.broadcast_to(p, (reps, len(p)) if p.ndim == 1 else p.shape), 0, 1)
        return p.copy() if depth is None else rng.binomial(depth, p) / depth

    rows, tops, cfg_rows = [], [], []
    for k in run:
        fd = folds[k]
        rng = np.random.default_rng(task + [k])
        feat = cl.fold_features(args.features, fd, T, nor, blood, a_idx, N_TUMOUR_MARKERS)
        nf = len(feat)
        # Work on the CpGs this fold needs: its features plus the atlas CpGs (for NNLS).
        use = np.union1d(feat, a_idx)
        fpos, apos = np.searchsorted(use, feat), np.searchsorted(use, a_idx)
        Bu = {g: v[use] for g, v in B.items()}
        Tu = {c: v[use] for c, v in T.items()}
        Nu = nor[use]

        # Training pool: nulls per person; positives over persons x cases x f (subsampled).
        Xs, meta = [], []
        for b in fd["train_bg"]:
            Xs.append(observe(Bu[b], rng, 1 if depth is None else NULL_TRAIN))
            meta += [(b, "", 0.0)] * len(Xs[-1])
        trip = cl.training_triples(fd, F_TRAIN)
        pick = rng.choice(len(trip), size=min(args.max_pos, len(trip)), replace=False)
        P = np.stack([cl.spike(Bu[trip[i][0]], Tu[trip[i][1]], purity[trip[i][1]],
                               trip[i][2]) for i in pick])
        Xs.append(observe(P, rng))
        meta += [trip[i] for i in pick]
        Xfull = np.concatenate(Xs)  # columns = use
        M = pd.DataFrame(meta, columns=["person", "case", "f"])
        y = (M.f > 0).to_numpy().astype(int)
        null = y == 0  # true nulls; thresholds always use these
        if args.permute:  # leakage control: labels carry no information
            y = rng.permutation(y)
        X = Xfull[:, fpos]

        # Inner folds: same paired design inside the training part.
        tr_groups = sorted({group_of[b] for b in fd["train_bg"]})
        inner = cl.paired_folds([bg_groups[g] for g in tr_groups], fd["train_cases"], rng)
        inner_masks = cl.inner_masks(M, inner)

        # Test set: held-out persons x held-out cases x F_EVAL; nulls; external plasmas.
        tX, tmeta = [], []
        for b in fd["test_bg"]:
            tX.append(observe(Bu[b], rng, 1 if depth is None else NULL_TEST))
            tmeta += [("neg", b, "", 0.0)] * len(tX[-1])
            for c in fd["test_cases"]:
                tX.append(observe(np.stack([cl.spike(Bu[b], Tu[c], purity[c], f)
                                            for f in F_EVAL]), rng))
                tmeta += [("pos", b, c, f) for f in F_EVAL]
            # Batch control: TCGA solid-tissue normals spiked like tumours (purity 1).
            for j in range(Nu.shape[1]):
                tX.append(observe(np.stack([cl.spike(Bu[b], Nu[:, j], 1.0, f)
                                            for f in F_NORMAL]), rng))
                tmeta += [("nrm", b, f"normal{j}", f) for f in F_NORMAL]
            for pname, g in PLASMAS.items():
                fs = [f for f in F_EVAL if f <= c_pl[pname]]
                for f in fs:
                    tX.append(observe(cl.spike(Bu[b], Bu[g], c_pl[pname], f), rng,
                                      1 if depth is None else EXT_REPS))
                    tmeta += [("ext", b, pname, f)] * len(tX[-1])
        # Blank-into-blank control. Its own generator, so the draws above (and so every
        # NNLS and LLR score) stay identical to the first version.
        rng_b = np.random.default_rng(task + [k, 2])
        for b in fd["test_bg"]:
            for g in bg_groups:
                if b in g or HI1 in g:
                    continue
                for d in g:
                    tX.append(observe(np.stack([cl.spike(Bu[b], Bu[d], 1.0, f)
                                                for f in F_B2B]), rng_b))
                    tmeta += [("b2b", b, d, f) for f in F_B2B]
        tXfull = np.concatenate(tX)
        TM = pd.DataFrame(tmeta, columns=["kind", "person", "source", "f"])
        tXf = tXfull[:, fpos]

        def emit(model, score, thr, prob=None, f_alt=np.nan, TM=TM, k=k):
            for i, r in enumerate(TM.itertuples(index=False)):
                rows.append((args.track, args.features, args.hi1, k, model, r.kind, r.person,
                             r.source, r.f, f_alt, float(score[i]), thr,
                             np.nan if prob is None else float(prob[i])))

        # Baseline NNLS colon fraction (atlas CpGs; nothing is fitted to blanks).
        if args.features == "atlas":
            per = [sim.nnls_colon(Xfull[(M.person == b).to_numpy() & null][:, apos], a_std, ic)
                   for b in fd["train_bg"]]
            thr = sim.prediction_bound(per)
            emit("NNLS", sim.nnls_colon(tXfull[:, apos], a_std, ic), thr)

        # Baseline LLR (read tracks): p0 = training-blank mean, p1 = p0 + f * mean tumour
        # shift per unit tumour-cell fraction. Inner leave-one-group-out for the threshold.
        if depth is not None:
            def llr_parts(persons, cases, feat=feat):
                p0 = np.mean([B[b][feat] for b in persons], axis=0)
                dl = np.nanmean(np.stack([np.where(np.isnan(T[c][feat]), 0.0,
                                                   (T[c][feat] - p0) / purity[c])
                                          for c in cases]), axis=0)
                return p0, dl
            p0, dl = llr_parts(fd["train_bg"], fd["train_cases"])
            n_t = np.full(tXf.shape, depth)
            k_t = np.rint(tXf * depth)
            for f in F_EVAL:
                per = []
                for g in tr_groups:
                    rest = [b for b in fd["train_bg"] if group_of[b] != g]
                    p0j, dlj = llr_parts(rest, fd["train_cases"])
                    for b in bg_groups[g]:
                        xb = X[(M.person == b).to_numpy() & null]
                        per.append(sim.llr(np.rint(xb * depth), np.full(xb.shape, depth), p0j,
                                           np.clip(p0j + f * dlj, 0, 1)))
                thr = sim.prediction_bound(per)
                s = sim.llr(k_t, n_t, p0, np.clip(p0 + f * dl, 0, 1))
                emit("LLR", s, thr, f_alt=f)

        # Trained models: tune on the inner folds, threshold from inner out-of-fold nulls.
        for name in trained:
            grid = {"XGB": XGB_GRID, "LR_L2": L2_GRID, "LR_L1": L1_GRID}[name]
            best = None
            for cfg in grid:
                oof = np.full(len(y), np.nan)
                fits = []
                for trm, te in inner_masks:
                    fits.append(Fitted(name, cfg, X[trm], y[trm], args.threads))
                    oof[te] = fits[-1].margin(X[te])
                neg = oof[null & ~np.isnan(oof)]
                pos = oof[M.f.isin(F_TUNE).to_numpy() & ~np.isnan(oof)]
                a = cl.auc(pos, neg)
                cfg_rows.append((args.track, args.features, args.hi1, k, name, str(cfg), a))
                if best is None or a > best[0]:
                    best = (a, cfg, oof, fits)
            _, cfg, oof, fits = best
            per = [oof[(M.person == b).to_numpy() & null] for b in fd["train_bg"]]
            thr = sim.prediction_bound(per)
            emit(name, np.mean([m.margin(tXf) for m in fits], axis=0), thr,
                 prob=np.mean([m.prob(tXf) for m in fits], axis=0))
            fm = Fitted(name, cfg, X, y, args.threads)
            s_ref = fm.margin(tXf)
            emit(f"{name}_refit", s_ref, thr, prob=fm.prob(tXf))
            mapped = np.full(len(y), np.nan)
            for m, (trm, te) in zip(fits, inner_masks):
                sh = trm & ~null
                slope, icpt = np.polyfit(m.margin(X[sh]), fm.margin(X[sh]), 1)
                mapped[te] = slope * oof[te] + icpt
            thr_c = sim.prediction_bound([mapped[(M.person == b).to_numpy() & null]
                                          for b in fd["train_bg"]])
            emit(f"{name}_refitcal", s_ref, thr_c)
            tops += [(args.track, args.features, args.hi1, k, name, r, cpgs[feat[i]])
                     for r, i in enumerate(fm.top())]
        print(f"  fold {k}: train n={len(y)} (neg {int((y == 0).sum())}), test n={len(TM)}, "
              f"features {nf}, {time.time() - t0:.0f}s", flush=True)

    cols = ["track", "features", "hi1", "fold", "model", "kind", "person", "source", "f",
            "f_alt", "score", "threshold", "prob"]
    tag = f"{args.track}_{args.features}_{args.hi1}"
    if args.folds:
        tag += "_f" + args.folds.replace(",", "-")
    if args.permute:
        tag += "_perm"
    pd.DataFrame(rows, columns=cols).to_csv(args.out / f"scores_{tag}.tsv.gz", sep="\t",
                                            index=False)
    pd.DataFrame(tops, columns=["track", "features", "hi1", "fold", "model", "rank", "cpg"]
                 ).to_csv(args.out / f"top_{tag}.tsv", sep="\t", index=False)
    pd.DataFrame(cfg_rows, columns=["track", "features", "hi1", "fold", "model", "config",
                                    "inner_auc"]).to_csv(args.out / f"tune_{tag}.tsv",
                                                         sep="\t", index=False)
    print(f"wrote {tag} rows={len(rows)} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
