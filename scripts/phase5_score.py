"""Phase 5: read-level mixing, scores, thresholds and detection (pipeline/PREREGISTERED.md).

Inputs (pipeline outputs under $CFDNA_DATA/phase5): pat/<run>.{dedup,withdup}.pat.gz,
blocks_cpg.tsv, offmarker_cpg.tsv; transfer models from scripts/phase5_transfer_train.py.

Steps:
  1. Fragment records per run and block set (atlas, tumour, off-marker), cached as .npz.
  2. Background profile p0 per block = pooled methylation of the 14 threshold donors;
     beta-binomial rho per set fitted on the threshold donors.
  3. Threshold pools: all k-subsets of the threshold donors for k = 1, else 50 random
     k-subsets (seed 20261010). Each pool is scored at f = 0 with p0 rebuilt WITHOUT the
     pool's donors (out-of-sample null scores). Threshold = src/phase5.prediction_bound_pools.
  4. Eval tests: each eval pool x f grid x 5 seeds, spike fragments sampled from the spike
     source (src/phase5 module docstring). f = 0 rows give the held-out FPR.
  5. Binomial gap ladder at matched k (binomial LLR only):
       rung 1 "binom_mean": Binomial counts, every block with the mean real calls per block;
       rung 2 "binom_real": Binomial counts at each test's real calls per block;
       rung 3 "real":       real fragments, duplicates removed (the main result);
       rung 4 "real_dup":   real fragments, duplicates kept.
     Rungs 1-2 draw methylated calls ~ Binomial(n_b, (1 - f) p0_b + f h_b), h_b = spike
     source block methylation; their thresholds come from the same draws at f = 0.
  6. External check: the 10 cancers vs the 16 eval donors at k = 1, k = 1 thresholds.

Writes results/phase5/scores.tsv.gz (one row per test x score) and
results/phase5/run_info.tsv (p0 coverage, rho, spike q range).

Usage: uv run python scripts/phase5_score.py [--seeds 5] [--workers 8]
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import phase5 as p5

SEED = 20261010
KS = [1, 2, 4, 8]
F_GRID = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
N_THRESHOLD_POOLS = 50
SETS = ["atlas", "tumour"]
TRANSFER = ["LR_L2", "LR_L1", "XGB"]


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def load_blocks(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None,
                       names=["chrom", "start", "end", "startCpG", "endCpG", "block", "set", "n"])


def build_records(job) -> str:
    """Parse one pat file into records per set; cache as .npz. Returns the cache path."""
    pat_path, cache, ranges = job
    if not Path(cache).exists():
        pat = p5.read_pat(Path(pat_path))
        out = {"counts": pat["count"].to_numpy()}
        for s, (st, en) in ranges.items():
            r = p5.pat_records(pat, st, en)
            for k in ("row", "block", "m", "u", "w"):
                out[f"{s}_{k}"] = getattr(r, k)
        np.savez_compressed(cache + ".tmp.npz", **out)
        os.replace(cache + ".tmp.npz", cache)
    return cache


def load_cache(path: str, sets) -> tuple[dict, np.ndarray]:
    z = np.load(path)
    recs = {s: p5.Records(*(z[f"{s}_{k}"] for k in ("row", "block", "m", "u", "w")))
            for s in sets}
    return recs, z["counts"]


class Scorer:
    """Every pre-registered score of one sample, given its records per set."""

    def __init__(self, nb, tumour_mean, rho, transfer, probe_block):
        self.nb, self.T, self.rho = nb, tumour_mean, rho
        self.transfer, self.probe_block = transfer, probe_block

    def llr_only(self, M, U, p0, s):
        ok = np.isfinite(p0)
        return [("LLR_binom", s, f, p5.llr_binomial(M[ok], U[ok], p0[ok],
                                                    p5.p_alt(p0, self.T[s], f)[ok]))
                for f in F_GRID]

    def all(self, recs, p0, full=True):
        out = []
        for s in SETS:
            r, q0 = recs[s], p0[s]
            M, U = p5.block_sums(r, self.nb[s])
            out += self.llr_only(M, U, q0, s)
            if not full:
                continue
            ok = np.isfinite(q0)
            rk = r.take(ok[r.block])
            for f in F_GRID:
                q1 = p5.p_alt(q0, self.T[s], f)
                out.append(("LLR_bb", s, f, p5.llr_betabinomial(rk, q0, q1, self.rho[s])))
            beta = np.where(M + U > 0, M / np.maximum(M + U, 1), np.nan)
            for name, mdl in self.transfer.get(s, {}).items():
                b = self.probe_block[s]
                x = np.where(b >= 0, beta[np.maximum(b, 0)], np.nan)
                x = np.where(np.isnan(x), mdl["train_mean"], x)[None]
                out.append((f"T_{name}", s, np.nan,
                            float(np.mean([m.margin(x)[0] for m in mdl["fits"]]))))
        if full:
            out.append(("UXM", "tumour", np.nan, p5.uxm_share(recs["tumour"])))
        return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    ap.add_argument("--samples", type=Path, default=Path("pipeline/samples/phase5_samples.tsv"))
    ap.add_argument("--markers", type=Path, default=Path("pipeline/markers/phase5_markers.tsv"))
    ap.add_argument("--spike", default="SRR28401050")
    ap.add_argument("--out", type=Path, default=Path("results/phase5"))
    ap.add_argument("--no-transfer", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    data = data_dir() / "phase5"
    args.out.mkdir(parents=True, exist_ok=True)
    cache_dir = data / "records"
    cache_dir.mkdir(exist_ok=True)

    smp = pd.read_csv(args.samples, sep="\t")
    blocks = load_blocks(data / "blocks_cpg.tsv")
    off = load_blocks(data / "offmarker_cpg.tsv")
    bsets = {s: blocks[blocks.set == s].reset_index(drop=True) for s in SETS}
    bsets["off"] = off
    nb = {s: len(b) for s, b in bsets.items()}
    ranges = {s: (b.startCpG.to_numpy(), b.endCpG.to_numpy()) for s, b in bsets.items()}
    mk = pd.read_csv(args.markers, sep="\t")
    tumour_mean = {}
    for s in SETS:
        tb = mk[mk.set == s].groupby("block").tumour_mean.mean()
        tumour_mean[s] = tb.reindex(bsets[s].block).to_numpy()
        assert np.isfinite(tumour_mean[s]).all()

    # 1. Records.
    runs = smp.run.tolist() + [args.spike]
    jobs = [(str(data / "pat" / f"{r}.{v}.pat.gz"), str(cache_dir / f"{r}.{v}.npz"), ranges)
            for r in runs for v in ("dedup", "withdup")]
    with ProcessPoolExecutor(args.workers) as ex:
        list(ex.map(build_records, jobs))
    R = {v: {} for v in ("dedup", "withdup")}
    counts = {}
    for r in runs:
        for v, rv in R.items():
            rv[r], c = load_cache(str(cache_dir / f"{r}.{v}.npz"), SETS + ["off"])
            counts[(r, v)] = c
    print(f"records loaded {time.time() - t0:.0f}s", flush=True)

    thr_donors = smp.run[smp.role == "threshold"].tolist()
    eval_donors = smp.run[smp.role == "eval"].tolist()
    cancers = smp.run[smp.role == "cancer"].tolist()
    assert len(thr_donors) == 14 and len(eval_donors) == 16 and len(cancers) == 10
    assert not set(thr_donors) & set(eval_donors)

    # 2. p0 from per-donor block sums (dedup); leave-pool-out versions by subtraction.
    sums = {v: {r: {s: p5.block_sums(R[v][r][s], nb[s]) for s in SETS}
                for r in runs} for v in R}

    def p0_of(donors):
        out = {}
        for s in SETS:
            M = sum(sums["dedup"][d][s][0] for d in donors)
            N = sum(sum(sums["dedup"][d][s]) for d in donors)
            out[s] = np.where(N > 0, M / np.maximum(N, 1), np.nan)
        return out

    p0 = p0_of(thr_donors)
    rho = {s: p5.fit_rho(p5.Records.concat([R["dedup"][d][s] for d in thr_donors]), p0[s])
           for s in SETS}
    h, h_cov = {}, {}  # spike-source block methylation; p0 where the spike has no call
    for s in SETS:
        M, U = sums["dedup"][args.spike][s]
        h_cov[s] = float((M + U > 0).mean())
        h[s] = np.where(M + U > 0, M / np.maximum(M + U, 1), p0[s])

    transfer, probe_block = {}, {}
    if not args.no_transfer:
        tdir = data / "transfer"
        for s in SETS:
            transfer[s] = {}
            blk_idx = {b: i for i, b in enumerate(bsets[s].block)}
            pb = dict(zip(mk.probe[mk.set == s], mk.block[mk.set == s]))
            for name in TRANSFER:
                with open(tdir / f"{s}_{name}.pkl", "rb") as fh:
                    transfer[s][name] = pickle.load(fh)
            probes = transfer[s][TRANSFER[0]]["probes"]
            probe_block[s] = np.array([blk_idx.get(pb.get(p), -1) for p in probes])
    sc = Scorer(nb, tumour_mean, rho, transfer, probe_block)

    rows = []

    def emit(rung, kind, k, pool, f, seed, scores):
        for name, s, fa, val in scores:
            rows.append((rung, kind, k, pool, f, seed, name, s, fa, val))

    def pool_recs(donors, v):
        return {s: p5.Records.concat([R[v][d][s] for d in donors]) for s in SETS + ["off"]}

    def mix(bg, f, rng, v):
        """Background pool records plus spike fragments at share f; returns (recs, q)."""
        if f == 0:
            return bg, 0.0
        sp = R[v][args.spike]
        q = f / (1 - f) * p5.n_fragments(bg["off"]) / p5.n_fragments(sp["off"])
        kept = p5.thin_rows(counts[(args.spike, v)], q, rng)
        return {s: p5.Records.concat([bg[s], p5.keep_rows(sp[s], kept)])
                for s in bg}, q

    # 3. Threshold pools.
    rng0 = np.random.default_rng([SEED, 0])
    tpools = {}
    for k in KS:
        if k == 1:
            tpools[k] = [(d,) for d in thr_donors]
        else:
            seen = set()
            while len(seen) < N_THRESHOLD_POOLS:
                seen.add(tuple(sorted(rng0.choice(thr_donors, size=k, replace=False))))
            tpools[k] = sorted(seen)
    for v, rung in (("dedup", "real"), ("withdup", "real_dup")):
        for k in KS:
            for i, pool in enumerate(tpools[k]):
                rest = [d for d in thr_donors if d not in pool]
                emit(rung, "threshold", k, i, 0.0, 0,
                     sc.all(pool_recs(pool, v), p0_of(rest), full=v == "dedup"))
    print(f"threshold pools scored {time.time() - t0:.0f}s", flush=True)

    # 4. Eval tests; 5. ladder rungs.
    qs = []
    pool_cols = {k: f"pool_k{k}" for k in KS}
    ev = smp[smp.role == "eval"]
    mean_calls = {}
    for k in KS:
        pools = sorted(ev[pool_cols[k]].unique())
        ncalls = []
        for pid in pools:
            donors = ev.run[ev[pool_cols[k]] == pid].tolist()
            assert len(donors) == k
            for v, rung in (("dedup", "real"), ("withdup", "real_dup")):
                bg = pool_recs(donors, v)
                for fi, f in enumerate([0.0, *F_GRID]):
                    for seed in range(args.seeds if f else 1):
                        rng = np.random.default_rng([SEED, 1, k, pid, fi, seed,
                                                     int(v == "withdup")])
                        recs, q = mix(bg, f, rng, v)
                        qs.append((v, k, pid, f, seed, q))
                        if v == "dedup":
                            emit(rung, "eval", k, pid, f, seed, sc.all(recs, p0))
                            nbs = {s: sum(p5.block_sums(recs[s], nb[s])) for s in SETS}
                            if f == 0:
                                ncalls.append({s: nbs[s].mean() for s in SETS})
                            # rung 2: binomial at this test's real calls per block
                            rng2 = np.random.default_rng([SEED, 2, k, pid, fi, seed])
                            sim = []
                            for s in SETS:
                                pm = np.clip((1 - f) * p0[s] + f * h[s], 0, 1)
                                n = nbs[s].astype(np.int64)
                                Ms = rng2.binomial(n, np.nan_to_num(pm))
                                sim += sc.llr_only(Ms, n - Ms, p0[s], s)
                            emit("binom_real", "eval", k, pid, f, seed, sim)
                        else:
                            emit(rung, "eval", k, pid, f, seed, sc.all(recs, p0, full=False))
        mean_calls[k] = {s: float(np.mean([c[s] for c in ncalls])) for s in SETS}
        # rung 2 nulls: threshold pools' real calls per block
        for i, pool in enumerate(tpools[k]):
            rng2 = np.random.default_rng([SEED, 3, k, i])
            sim = []
            for s in SETS:
                n = sum(sum(sums["dedup"][d][s]) for d in pool).astype(np.int64)
                Ms = rng2.binomial(n, np.nan_to_num(p0[s]))
                sim += sc.llr_only(Ms, n - Ms, p0[s], s)
            emit("binom_real", "threshold", k, i, 0.0, 0, sim)
        # rung 1: every block at the mean real calls per block (stochastic rounding)
        for kind, n_tests in (("threshold", N_THRESHOLD_POOLS), ("eval", None)):
            tests = ([(i, 0.0, 0) for i in range(n_tests)] if kind == "threshold" else
                     [(pid, f, seed) for pid in pools for f in [0.0, *F_GRID]
                      for seed in range(args.seeds if f else 1)])
            for j, (pid, f, seed) in enumerate(tests):
                rng1 = np.random.default_rng([SEED, 4, k, int(kind == "eval"), j])
                sim = []
                for s in SETS:
                    nbar = mean_calls[k][s]
                    n = np.floor(nbar) + (rng1.random(nb[s]) < nbar - np.floor(nbar))
                    n = n.astype(np.int64)
                    pm = np.clip((1 - f) * p0[s] + f * h[s], 0, 1)
                    Ms = rng1.binomial(n, np.nan_to_num(pm))
                    sim += sc.llr_only(Ms, n - Ms, p0[s], s)
                emit("binom_mean", kind, k, pid, f, seed, sim)
        print(f"k={k} done {time.time() - t0:.0f}s", flush=True)

    # 6. External check: cancers at k = 1 (eval donors at f = 0 are the negatives above).
    for c in cancers:
        emit("real", "cancer", 1, c, np.nan, 0, sc.all(pool_recs([c], "dedup"), p0))

    cols = ["rung", "kind", "k", "pool", "f", "seed", "score", "set", "f_alt", "value",
            ]
    out = pd.DataFrame(rows, columns=cols)
    out.to_csv(args.out / "scores.tsv.gz", sep="\t", index=False)
    Q = pd.DataFrame(qs, columns=["variant", "k", "pool", "f", "seed", "q"])
    Q.to_csv(args.out / "spike_q.tsv", sep="\t", index=False)
    info = []
    for s in SETS:
        info += [("p0_blocks_covered", s, float(np.isfinite(p0[s]).mean())),
                 ("rho", s, rho[s]),
                 ("spike_blocks_covered", s, h_cov[s])]
        for k in KS:
            info.append((f"mean_calls_per_block_k{k}", s, mean_calls[k][s]))
    info.append(("q_max", "", float(Q.q.max())))
    pd.DataFrame(info, columns=["item", "set", "value"]).to_csv(
        args.out / "run_info.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"wrote {len(out)} rows in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
