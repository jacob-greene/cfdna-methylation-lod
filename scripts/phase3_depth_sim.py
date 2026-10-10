"""Phase 3: sequencing-depth simulation of colon-tumour detection at the atlas CpGs.

One run = one (spike kind, depth D, input cap, background mode). It writes per-fold
detection counts for both untrained scores; scripts/phase3_summary.py aggregates them.

Model (Phase 3): methylated reads per CpG ~ Binomial(n, f * beta_t + (1 - f) * beta_b)
  - n = D, or min(D, unique molecules) under an input cap (see --cap);
  - beta_b: one real healthy blank array per simulated sample (between-person variation kept);
  - beta_t, two spike kinds:
      cell:   a colon-epithelium replicate (GSM3455852/63/64). The detector's tumour model and
              the NNLS atlas colon column are the mean of the other two replicates.
      plasma: the pl806 colon-cancer plasma, de-mixed to a pure-tumour profile with
              c(pl806) and the training-blank mean. Truth = one pl806 array (250 ng or 50 ng);
              the detector's tumour model = the other array.
    f is the colon / tumour DNA fraction, as in Phase 2.

Scores: NNLS colon fraction of the simulated betas; summed binomial LLR of the alternative
"tumour at fraction f" (tumour model, f fixed at the evaluated level) against "training-blank
mean". Each f has its own LLR and its own threshold.

Threshold rules: "pooled" = EP17 0.95 rank of the training blanks' null simulations, scored
by the fold's own model (in-sample). "valid" = each training blank's null simulations scored by
a model rebuilt without it (inner leave-one-out), then a normal prediction bound over people:
mean + t(0.95, n - 1) * sqrt(SD_between^2 * (1 + 1/n) + SD_within^2).

Folds (leave-one-blank-out): each blank in turn is the evaluation background. Markers,
background profile p0, de-mixing and the 95 %-specificity thresholds use only the other blanks
(a pool's rep1/rep2 twin is also left out). Thresholds = EP17 0.95 rank of f = 0 simulations
of the training blanks. Run with all 12 blanks and with HI1 removed.

Markers: the atlas CpGs present on EPIC with no missing value in the blanks or spikes, ranked
by |beta_t model - training mean| / sqrt(training SD^2 + 0.02^2); top 100, top 500, all.

Input cap: unique molecules per CpG ~ Poisson(ng * 303 * recovery); 303 haploid genomes per
ng = 3.3 pg per haploid genome (Moss et al. 2018, Methods). Reads beyond that are duplicates.

Background mode "mean" (control): every simulated sample uses the training-blank mean as
beta_b, so between-person variation is removed and only binomial noise remains.

Usage: uv run python scripts/phase3_depth_sim.py --kind plasma --depth 100 [--cap none]
       [--bg person] [--reps 100] [--out results/phase3/raw]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import deconv, lod
from src import simulate as sim

SEED = 20261008
COLON = "Colon_epithelial_cells"
POOLS = ["GSM3455776", "GSM3455777", "GSM3455778", "GSM3455779",
         "GSM3455780", "GSM3455781", "GSM3455782", "GSM3455783"]
POOL_PAIRS = [("GSM3455776", "GSM3455777"), ("GSM3455778", "GSM3455779"),
              ("GSM3455780", "GSM3455781"), ("GSM3455782", "GSM3455783")]
HI1 = "GSM3455807"
INDIV = [HI1, "GSM3455840", "GSM3455846", "GSM3455848"]
BLANKS = POOLS + INDIV
COLON_REPS = ["GSM3455852", "GSM3455863", "GSM3455864"]
PL806 = ["GSM3455808", "GSM3455810"]
F_GRID = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]  # study f grid
# Added between grid points so the 95 % crossing is interpolated over a short span.
F_FINE = [5e-4, 2e-3, 5e-3, 1.5e-2, 2e-2]
F_ALL = sorted(F_GRID + F_FINE)
MARKER_SETS = {"100": 100, "500": 500, "all": None}
CAPS = {  # name -> (ng input, recovery after bisulfite)
    "none": None,
    "30ng_r0.5": (30, 0.5),
    "10ng_r0.5": (10, 0.5),
    "10ng_r0.1": (10, 0.1),
}
GE_PER_NG = 303  # haploid genome equivalents per ng (Moss et al. 2018, Methods)


def sel(order: np.ndarray, m: int | None) -> np.ndarray:
    """Top-m marker indices of a ranking (all when m is None)."""
    return order if m is None else order[:m]


def default_data() -> Path:
    env = os.environ.get("CFDNA_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def load(data: Path):
    atl = deconv.dedup_atlas(atlas_mod.read_atlas(data))
    epic = pd.read_parquet(data / "geo" / "betas_epic.parquet").reindex(atl.index)
    need = BLANKS + COLON_REPS + PL806
    ok = epic[need].notna().all(axis=1) & atl.notna().all(axis=1)
    atl, epic = atl[ok], epic.loc[ok, need].clip(0, 1)
    return atl, epic


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["cell", "plasma"], required=True)
    ap.add_argument("--depth", type=int, required=True)
    ap.add_argument("--cap", choices=list(CAPS), default="none")
    ap.add_argument("--bg", choices=["person", "mean"], default="person")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--data", type=Path, default=default_data())
    ap.add_argument("--out", type=Path, default=Path("results/phase3/raw"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    atl, epic = load(args.data)
    cols = list(atl.columns)
    ic = cols.index(COLON)
    a_std = atl.to_numpy(dtype=float)
    B = {g: epic[g].to_numpy(dtype=float) for g in epic.columns}
    n_cpg = len(atl)
    # c(pl806): mean NNLS colon fraction of the two pl806 arrays, standard atlas (Phase 2).
    c806 = float(np.mean(sim.nnls_colon(np.stack([B[g] for g in PL806]), a_std, ic)))
    cap = CAPS[args.cap]
    cap_mean = None if cap is None else cap[0] * GE_PER_NG * cap[1]
    twin = {a: b for a, b in POOL_PAIRS} | {b: a for a, b in POOL_PAIRS}
    truths = COLON_REPS if args.kind == "cell" else PL806
    task_key = [SEED, {"cell": 1, "plasma": 2}[args.kind], args.depth,
                list(CAPS).index(args.cap), {"person": 0, "mean": 1}[args.bg]]
    print(f"kind={args.kind} D={args.depth} cap={args.cap} ({cap_mean}) bg={args.bg} "
          f"CpGs={n_cpg} c806={c806:.4f}", flush=True)

    rows = []
    for hi1 in ("with", "without"):
        bset = BLANKS if hi1 == "with" else [g for g in BLANKS if g != HI1]
        for it, truth in enumerate(truths):
            for ie, ev in enumerate(bset):
                rng = np.random.default_rng(task_key + [int(hi1 == "with"), it, ie])
                train = [g for g in bset if g not in (ev, twin.get(ev))]
                T = np.stack([B[g] for g in train])
                mu = T.mean(axis=0)
                if args.kind == "cell":
                    others = [r for r in COLON_REPS if r != truth]
                    model = np.mean([B[r] for r in others], axis=0)
                    beta_t = B[truth]
                    a_fold = a_std.copy()
                    a_fold[:, ic] = model
                else:
                    other = next(g for g in PL806 if g != truth)
                    model = sim.demix(B[other], mu, c806)
                    beta_t = sim.demix(B[truth], mu, c806)
                    a_fold = a_std
                order = sim.rank_markers(model, T)
                bg_eval = mu if args.bg == "mean" else B[ev]
                bg_train = [mu] * len(train) if args.bg == "mean" else [B[g] for g in train]

                # Simulate once on all CpGs; marker sets are column subsets.
                n_null = np.concatenate([sim.molecules(args.depth, cap_mean, n_cpg,
                                                       args.reps, rng) for _ in bg_train])
                k_null = np.concatenate([
                    sim.simulate_counts(b, None, 0.0, n_null[i * args.reps:(i + 1) * args.reps],
                                        rng) for i, b in enumerate(bg_train)])
                ev_n, ev_k = {}, {}
                for f in [0.0, *F_ALL]:
                    n = sim.molecules(args.depth, cap_mean, n_cpg, args.reps, rng)
                    ev_n[f] = n
                    ev_k[f] = sim.simulate_counts(bg_eval, beta_t, f, n, rng)

                # Inner leave-one-out: each training blank j is scored by a model built
                # without j (and its twin), so its null scores are out-of-sample.
                inner = []
                for j, g in enumerate(train):
                    tr_j = [h for h in train if h not in (g, twin.get(g))]
                    mu_j = np.stack([B[h] for h in tr_j]).mean(axis=0)
                    model_j = model if args.kind == "cell" else sim.demix(
                        B[next(h for h in PL806 if h != truth)], mu_j, c806)
                    order_j = sim.rank_markers(model_j, np.stack([B[h] for h in tr_j]))
                    sl = slice(j * args.reps, (j + 1) * args.reps)
                    inner.append((sl, mu_j, model_j, order_j))

                for mname, m in MARKER_SETS.items():
                    idx = sel(order, m)
                    a_m = a_fold[idx]
                    # Thresholds per rule. pooled: EP17 0.95 rank of in-sample null scores.
                    # valid: normal prediction bound over out-of-sample null scores.
                    # NNLS colon fraction: one threshold for all f.
                    per = []
                    for sl, _, _, o in inner:
                        ii = sel(o, m)
                        per.append(sim.nnls_colon(k_null[sl][:, ii] / n_null[sl][:, ii],
                                                  a_fold[ii], ic))
                    thr = {"pooled": lod.lob_nonparametric(
                        sim.nnls_colon(k_null[:, idx] / n_null[:, idx], a_m, ic)),
                           "valid": sim.prediction_bound(per)}
                    for f in [0.0, *F_ALL]:
                        s = sim.nnls_colon(ev_k[f][:, idx] / ev_n[f][:, idx], a_m, ic)
                        for rule, th in thr.items():
                            rows.append((hi1, truth, ev, mname, "NNLS", f, len(s),
                                         int((s > th).sum()), th, float(np.median(s)), rule))
                    # LLR: alternative fixed at each evaluated f; threshold per f.
                    p0 = mu[idx]
                    for f in F_ALL:
                        p1 = f * model[idx] + (1 - f) * p0
                        per = []
                        for sl, mu_j, model_j, o in inner:
                            ii = sel(o, m)
                            per.append(sim.llr(k_null[sl][:, ii], n_null[sl][:, ii], mu_j[ii],
                                               f * model_j[ii] + (1 - f) * mu_j[ii]))
                        thr = {"pooled": lod.lob_nonparametric(
                            sim.llr(k_null[:, idx], n_null[:, idx], p0, p1)),
                               "valid": sim.prediction_bound(per)}
                        for fe in (0.0, f):
                            s = sim.llr(ev_k[fe][:, idx], ev_n[fe][:, idx], p0, p1)
                            for rule, th in thr.items():
                                rows.append((hi1, truth, ev, mname,
                                             "LLR" if fe else "LLR_null", f, len(s),
                                             int((s > th).sum()), th, float(np.median(s)),
                                             rule))
            print(f"  hi1={hi1} truth={truth} done {time.time() - t0:.0f}s", flush=True)

    out = pd.DataFrame(rows, columns=["hi1", "truth", "eval_blank", "markers", "score", "f",
                                      "n", "n_detected", "threshold", "median_score",
                                      "rule"])
    for k, v in [("kind", args.kind), ("depth", args.depth), ("cap", args.cap),
                 ("bg", args.bg), ("reps", args.reps), ("n_cpg_all", n_cpg),
                 ("cap_mean", cap_mean if cap_mean else np.nan)]:
        out[k] = v
    name = f"{args.kind}_D{args.depth}_{args.cap}_{args.bg}.tsv.gz"
    out.to_csv(args.out / name, sep="\t", index=False)
    print(f"wrote {name} rows={len(out)} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
