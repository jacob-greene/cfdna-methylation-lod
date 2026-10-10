"""Phase 6: optimised scores WZ and FRAG on the Phase 5 tests (pipeline/PREREGISTERED_PHASE6.md).

Stage "threshold" uses the 14 threshold donors only: it fits the null model (src/phase6.py)
and scores the Phase 5 threshold pools leave-pool-out (every fitted value refitted without
the pool's donors). It needs no pre-registration.

Stage "eval" scores the eval tests (same pools, f grid, seeds and spike fragments as
scripts/phase5_score.py) and the 10 cancers with the model fitted on all 14 threshold
donors. It refuses to run unless pipeline/PREREGISTERED_PHASE6.md and the Phase 6 code are
committed and unchanged.

Rungs: "real" (real fragments, duplicates removed) and, for WZ only, "binom_real"
(methylated calls ~ Binomial(N_b, (1 - f) p0_b + f h_b) at each test's real calls per block,
h_b = spike block methylation; same weights as "real").

Writes results/phase6/{threshold_scores.tsv, model_info.tsv} (stage threshold) and
results/phase6/{eval_scores.tsv.gz, effective_share.tsv, confound_calls.tsv} (stage eval).

Usage: uv run python scripts/phase6_score.py --stage threshold|eval
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import phase5 as p5
from src import phase6 as p6

_spec = importlib.util.spec_from_file_location("phase5_score", ROOT / "scripts/phase5_score.py")
p5s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5s)

SEED, KS, F_GRID, SETS = p5s.SEED, p5s.KS, p5s.F_GRID, p5s.SETS
PREREG = ["pipeline/PREREGISTERED_PHASE6.md", "src/phase6.py", "scripts/phase6_score.py",
          "scripts/phase6_summary.py"]


def preregistered_sha() -> str:
    """HEAD SHA if the pre-registration and Phase 6 code are committed and unchanged."""
    git = lambda *a: subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True,
                                    text=True)
    for f in PREREG:
        if git("ls-files", "--error-unmatch", f).returncode:
            sys.exit(f"refusing eval stage: {f} is not committed")
    if git("diff", "--quiet", "HEAD", "--", *PREREG).returncode:
        sys.exit("refusing eval stage: pre-registered files differ from HEAD")
    return git("rev-parse", "HEAD").stdout.strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["threshold", "eval"], required=True)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--samples", type=Path, default=ROOT / "pipeline/samples/phase5_samples.tsv")
    ap.add_argument("--markers", type=Path, default=ROOT / "pipeline/markers/phase5_markers.tsv")
    ap.add_argument("--spike", default="SRR28401050")
    ap.add_argument("--out", type=Path, default=ROOT / "results/phase6")
    args = ap.parse_args()
    sha = preregistered_sha() if args.stage == "eval" else ""
    t0 = time.time()
    data = p5s.data_dir() / "phase5"
    args.out.mkdir(parents=True, exist_ok=True)

    smp = pd.read_csv(args.samples, sep="\t")
    blocks = p5s.load_blocks(data / "blocks_cpg.tsv")
    bsets = {s: blocks[blocks.set == s].reset_index(drop=True) for s in SETS}
    nb = {s: len(b) for s, b in bsets.items()}
    mk = pd.read_csv(args.markers, sep="\t")
    T = {s: mk[mk.set == s].groupby("block").tumour_mean.mean().reindex(bsets[s].block)
         .to_numpy() for s in SETS}
    thr = smp.run[smp.role == "threshold"].tolist()
    assert len(thr) == 14
    if args.stage == "threshold":
        runs = thr  # no eval donor, cancer or spike record is read in this stage
    else:
        runs = smp.run.tolist() + [args.spike]
    R, counts = {}, {}
    for r in runs:
        R[r], counts[r] = p5s.load_cache(str(data / "records" / f"{r}.dedup.npz"),
                                         SETS + ["off"])
    # Phase 5 rho (threshold donors), Phase 5 p0 coverage rule.
    p0_all = {s: p6._p0([p6.summarise(R[d][s], nb[s], 0.0) for d in thr]) for s in SETS}
    rho = {s: p5.fit_rho(p5.Records.concat([R[d][s] for d in thr]), p0_all[s]) for s in SETS}
    S = {r: {s: p6.summarise(R[r][s], nb[s], rho[s]) for s in SETS} for r in runs}
    full = {s: p6.fit([S[d][s] for d in thr], T[s]) for s in SETS}
    print(f"loaded and fitted {time.time() - t0:.0f}s; rho {rho}", flush=True)

    def scores(Sx, model, k):
        return [(name, s, fn(Sx[s], model[s], k)) for s in SETS
                for name, fn in (("WZ", p6.wz), ("FRAG", p6.frag))]

    rows = []
    if args.stage == "threshold":
        rng0 = np.random.default_rng([SEED, 0])  # Phase 5 threshold pools, same draws
        info = []
        for k in KS:
            if k == 1:
                pools = [(d,) for d in thr]
            else:
                seen = set()
                while len(seen) < p5s.N_THRESHOLD_POOLS:
                    seen.add(tuple(sorted(rng0.choice(thr, size=k, replace=False))))
                pools = sorted(seen)
            for i, pool in enumerate(pools):
                rest = [d for d in thr if d not in pool]
                model = {s: p6.fit([S[d][s] for d in rest], T[s]) for s in SETS}
                Sp = {s: sum(S[d][s] for d in pool) for s in SETS}
                for name, s, v in scores(Sp, model, k):
                    rows.append(("real", "threshold", k, i, 0.0, 0, name, s, v))
                rng2 = np.random.default_rng([SEED, 63, k, i])
                for s in SETS:
                    n = p6._get(Sp[s], "N").astype(np.int64)
                    Ms = rng2.binomial(n, np.nan_to_num(full[s].p0))
                    rows.append(("binom_real", "threshold", k, i, 0.0, 0, "WZ", s,
                                 p6.wz(Sp[s], full[s], k, M=Ms)))
                if k == 1:
                    for s in SETS:
                        info.append((s, pool[0], int(model[s].sel.sum()),
                                     model[s].phi_mean, model[s].s2))
            print(f"threshold k={k} done {time.time() - t0:.0f}s", flush=True)
        cols = ["rung", "kind", "k", "pool", "f", "seed", "score", "set", "value"]
        pd.DataFrame(rows, columns=cols).to_csv(args.out / "threshold_scores.tsv", sep="\t",
                                                index=False, float_format="%.6g")
        mi = []
        for s in SETS:
            m = full[s]
            mi += [(s, "blocks", nb[s]), (s, "blocks_p0", int(np.isfinite(m.p0).sum())),
                   (s, "rho", rho[s]), (s, "phi_mean", m.phi_mean),
                   (s, "phi_q90", float(np.quantile(m.phi, 0.9))),
                   (s, "frag_blocks_selected", int(m.sel.sum())),
                   (s, "frag_blocks_hyper", int((m.sel & (m.delta > 0)).sum())),
                   (s, "frag_rate_mean_selected", float(m.rate[m.sel].mean())),
                   (s, "frag_s2", m.s2)]
            loo = [x for x in info if x[0] == s]
            mi += [(s, "loo_frag_blocks_min", min(x[2] for x in loo)),
                   (s, "loo_frag_blocks_max", max(x[2] for x in loo)),
                   (s, "loo_s2_min", min(x[4] for x in loo)),
                   (s, "loo_s2_max", max(x[4] for x in loo))]
        pd.DataFrame(mi, columns=["set", "item", "value"]).to_csv(
            args.out / "model_info.tsv", sep="\t", index=False, float_format="%.6g")
        print(f"wrote {len(rows)} threshold rows {time.time() - t0:.0f}s")
        return

    # Eval stage: same mixes as scripts/phase5_score.py (dedup variant).
    sp = args.spike
    n_sp = p5.n_fragments(R[sp]["off"])
    h = {}
    for s in SETS:
        M, N = p6._get(S[sp][s], "M"), p6._get(S[sp][s], "N")
        h[s] = np.where(N > 0, M / np.maximum(N, 1), full[s].p0)
    ev = smp[smp.role == "eval"]
    eff = []
    for k in KS:
        col = f"pool_k{k}"
        for pid in sorted(ev[col].unique()):
            donors = ev.run[ev[col] == pid].tolist()
            assert len(donors) == k
            Sbg = {s: sum(S[d][s] for d in donors) for s in SETS}
            n_bg = p5.n_fragments_pool([R[d]["off"] for d in donors])
            for fi, f in enumerate([0.0, *F_GRID]):
                for seed in range(args.seeds if f else 1):
                    if f:
                        rng = np.random.default_rng([SEED, 1, k, pid, fi, seed, 0])
                        q = f / (1 - f) * n_bg / n_sp
                        kept = p5.thin_rows(counts[sp], q, rng)
                        Ssp = {s: p6.summarise(p5.keep_rows(R[sp][s], kept), nb[s], rho[s])
                               for s in SETS}
                    else:
                        Ssp = {s: np.zeros_like(Sbg[s]) for s in SETS}
                    Sx = {s: Sbg[s] + Ssp[s] for s in SETS}
                    for name, s, v in scores(Sx, full, k):
                        rows.append(("real", "eval", k, pid, f, seed, name, s, v))
                    rng2 = np.random.default_rng([SEED, 62, k, pid, fi, seed])
                    for s in SETS:
                        m = full[s]
                        n = p6._get(Sx[s], "N").astype(np.int64)
                        pm = np.clip((1 - f) * np.nan_to_num(m.p0) + f * h[s], 0, 1)
                        rows.append(("binom_real", "eval", k, pid, f, seed, "WZ", s,
                                     p6.wz(Sx[s], m, k, M=rng2.binomial(n, pm))))
                        ok = np.isfinite(m.p0)
                        Csp, _ = p6.frag_counts(Ssp[s], m)
                        eff.append((k, pid, f, seed, s,
                                    p6._get(Ssp[s], "N")[ok].sum()
                                    / p6._get(Sx[s], "N")[ok].sum(),
                                    p6._get(Ssp[s], "n3")[m.sel].sum()
                                    / p6._get(Sx[s], "n3")[m.sel].sum(), Csp))
        print(f"eval k={k} done {time.time() - t0:.0f}s", flush=True)
    for c in smp.run[smp.role == "cancer"]:
        for name, s, v in scores(S[c], full, 1):
            rows.append(("real", "cancer", 1, c, np.nan, 0, name, s, v))
    cols = ["rung", "kind", "k", "pool", "f", "seed", "score", "set", "value"]
    pd.DataFrame(rows, columns=cols).to_csv(args.out / "eval_scores.tsv.gz", sep="\t",
                                            index=False, float_format="%.6g")
    pd.DataFrame(eff, columns=["k", "pool", "f", "seed", "set", "share_calls",
                               "share_frag_n3_selected", "spike_tumour_like_selected"]).to_csv(
        args.out / "effective_share.tsv", sep="\t", index=False, float_format="%.6g")
    # Confound check: calls per fragment record, threshold donors vs spike source.
    cf = []
    for s in SETS:
        sel = full[s].sel
        for who, rr in (("threshold_donors", [R[d][s] for d in thr]), ("spike", [R[sp][s]])):
            r = p5.Records.concat(rr)
            for blk, mask in (("all", np.ones(nb[s], bool)), ("frag_selected", sel)):
                kk = mask[r.block]
                n, w = (r.m + r.u)[kk], r.w[kk]
                cf.append((s, who, blk, float((w * n).sum() / w.sum()),
                           float(w[n >= p6.MIN_CALLS].sum() / w.sum())))
    pd.DataFrame(cf, columns=["set", "source", "blocks", "mean_calls_per_record",
                              "share_records_ge3_calls"]).to_csv(
        args.out / "confound_calls.tsv", sep="\t", index=False, float_format="%.4g")
    (args.out / "eval_run_sha.txt").write_text(sha + "\n")
    print(f"wrote {len(rows)} eval rows at {sha[:7]} {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
