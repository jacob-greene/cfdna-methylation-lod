"""Phase 5: train the Phase 3b models once on the whole D = 30 simulation, for the transfer test.

Pre-registered rule (pipeline/PREREGISTERED.md, "Transfer models"): LR_L2, LR_L1 and XGB,
trained on the Phase 3b D = 30 read simulation with the non-HI1 blanks and all 276 kept
tumours, with the Phase 3b tuning. Applied unchanged to real reads with the Phase 3b
valid-rule threshold.

This is Phase 3b's training step with the outer fold = everything: the inner folds
(src/classify.paired_folds over all background groups and all tumour cases) tune the 6
configurations by pooled out-of-fold AUC at f = 1 % and 3 %; the chosen configuration's
inner-fold models score new samples by their mean margin; the threshold is the one-sided 95 %
normal prediction bound over each person's out-of-fold f = 0 margins (valid rule).

Feature sets: "atlas" = Phase 3b atlas CpGs; "tumour" = the top 500 tumour markers selected on
all 276 tumours (identical to pipeline/markers/phase5_markers.tsv, checked below).

Writes $CFDNA_DATA/phase5/transfer/<features>_<model>.pkl: dict with the fitted inner-fold
models, the threshold, the feature probe names, the training-pool feature means (used for
probes without a read call) and the tuning table.

Usage: uv run python scripts/phase5_transfer_train.py --features atlas [--threads 8]
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase3b_classify import (
    F_TRAIN,
    F_TUNE,
    HI2,
    HI3,
    HI4,
    L1_GRID,
    L2_GRID,
    N_TUMOUR_MARKERS,
    NULL_TRAIN,
    POOL_PAIRS,
    SEED,
    XGB_GRID,
    Fitted,
    data_dir,
)

from src import classify as cl
from src import simulate as sim

DEPTH = 30


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", choices=["atlas", "tumour"], required=True)
    ap.add_argument("--models", default="LR_L2,LR_L1,XGB")
    ap.add_argument("--max-pos", type=int, default=4000)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--markers", type=Path, default=Path("pipeline/markers/phase5_markers.tsv"))
    args = ap.parse_args()
    t0 = time.time()
    data = data_dir() / "phase3b"
    out = data.parent / "phase5" / "transfer"
    out.mkdir(parents=True, exist_ok=True)

    geo = pd.read_parquet(data / "geo.parquet")
    tum = pd.read_parquet(data / "tumours.parquet")
    nor = pd.read_parquet(data / "normals.parquet").to_numpy(dtype=float)
    blood = pd.read_parquet(data / "blood.parquet").to_numpy(dtype=float)
    tinfo = pd.read_csv(data / "tumours.tsv", sep="\t").set_index("id")
    cpgs = geo.index
    a_idx = cpgs.get_indexer((data / "atlas_cpgs.txt").read_text().split())
    bg_groups = [list(p) for p in POOL_PAIRS] + [[HI2], [HI3], [HI4]]
    group_of = {b: i for i, g in enumerate(bg_groups) for b in g}
    B = {g: geo[g].to_numpy(dtype=float) for g in geo.columns}
    case_ids = tinfo.case.tolist()
    T = {c: tum[fid].to_numpy(dtype=float) for fid, c in zip(tinfo.index, tinfo.case)}
    purity = dict(zip(tinfo.case, tinfo.CPE))
    every = {"train_bg": [b for g in bg_groups for b in g], "train_cases": case_ids}
    feat = cl.fold_features(args.features, every, T, nor, blood, a_idx, N_TUMOUR_MARKERS)
    probes = list(cpgs[feat])
    if args.features == "tumour":
        reg = pd.read_csv(args.markers, sep="\t")
        reg = set(reg.probe[reg.set == "tumour"])
        assert reg <= set(probes), "pre-registered tumour markers differ from this selection"
        print(f"tumour features: {len(probes)}; pre-registered with hg38 position: {len(reg)}")

    task = [SEED, 5, ["atlas", "tumour"].index(args.features)]
    rng = np.random.default_rng(task)
    Xs, meta = [], []
    for b in every["train_bg"]:
        p = np.broadcast_to(np.clip(B[b][feat], 0, 1), (NULL_TRAIN, len(feat)))
        Xs.append(rng.binomial(DEPTH, p) / DEPTH)
        meta += [(b, "", 0.0)] * NULL_TRAIN
    trip = cl.training_triples(every, F_TRAIN)
    pick = rng.choice(len(trip), size=min(args.max_pos, len(trip)), replace=False)
    P = np.stack([cl.spike(B[trip[i][0]][feat], T[trip[i][1]][feat], purity[trip[i][1]],
                           trip[i][2]) for i in pick])
    Xs.append(rng.binomial(DEPTH, np.clip(P, 0, 1)) / DEPTH)
    meta += [trip[i] for i in pick]
    X = np.concatenate(Xs)
    M = pd.DataFrame(meta, columns=["person", "case", "f"])
    y = (M.f > 0).to_numpy().astype(int)
    null = y == 0
    inner = cl.paired_folds(bg_groups, case_ids, rng)
    masks = cl.inner_masks(M, inner)
    print(f"features={args.features} n={len(y)} (neg {int(null.sum())}) inner folds "
          f"{len(masks)} groups {len(group_of)}", flush=True)

    for name in args.models.split(","):
        grid = {"XGB": XGB_GRID, "LR_L2": L2_GRID, "LR_L1": L1_GRID}[name]
        best, tune = None, []
        for cfg in grid:
            oof = np.full(len(y), np.nan)
            fits = []
            for trm, te in masks:
                fits.append(Fitted(name, cfg, X[trm], y[trm], args.threads))
                oof[te] = fits[-1].margin(X[te])
            neg = oof[null & ~np.isnan(oof)]
            pos = oof[M.f.isin(F_TUNE).to_numpy() & ~np.isnan(oof)]
            a = cl.auc(pos, neg)
            tune.append((str(cfg), a))
            if best is None or a > best[0]:
                best = (a, cfg, oof, fits)
        a, cfg, oof, fits = best
        per = [oof[(M.person == b).to_numpy() & null] for b in every["train_bg"]]
        thr = sim.prediction_bound(per)
        obj = {"model": name, "features": args.features, "config": cfg, "inner_auc": a,
               "fits": fits, "threshold": thr, "probes": probes,
               "train_mean": X.mean(axis=0), "tune": tune, "depth": DEPTH,
               "n_train": len(y), "n_null": int(null.sum())}
        with open(out / f"{args.features}_{name}.pkl", "wb") as fh:
            pickle.dump(obj, fh)
        print(f"  {name}: cfg={cfg} inner AUC={a:.4f} threshold={thr:.4f} "
              f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
