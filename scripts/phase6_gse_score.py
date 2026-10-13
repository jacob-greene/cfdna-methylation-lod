"""Phase 6 GSE186007 arm: LABS plasma runs as spikes in the Phase 5 eval pools (amendment 2).

Every GSE186007 run (8 CRC, 4 healthy "HEA"; phase6_gse/runs.tsv) is a spike source on its own.
Each spike is mixed into the same eval pools, f grid, seeds and dosing rule as rung "real" of
scripts/phase6_score.py: the spike's fragments are kept with probability
q = f / (1 - f) * n_bg / n_spike, counted in the off-marker blocks. A test with q > 1 is not
scored (NaN), as in rung "len". Scores are WZ and FRAG (src/phase6.py) with the model fitted
on all 14 threshold donors, so the f = 0 tests equal those of results/phase6 (checked).

Random seeds: default_rng([SEED, 2, k, pool, fi, seed, i]), i = index of the spike run in
phase6_gse/runs.tsv.

Each GSE run is also scored alone ("direct", k = 1, no SPOT-MAS background). These scores mix
protocol and tumour signal; they are descriptive only.

It refuses to run unless pipeline/PREREGISTERED_PHASE6.md and the Phase 6 code are committed
and unchanged (the guard of scripts/phase6_score.py, plus the files of this arm).

Writes results/phase6_gse/{eval_scores.tsv.gz, effective_share.tsv, direct_scores.tsv,
confound_calls.tsv, eval_run_sha.txt}.

Usage: uv run python scripts/phase6_gse_score.py
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
          "scripts/phase6_summary.py", "scripts/phase6_tf.py", "scripts/phase6_gse_score.py",
          "scripts/phase6_gse_summary.py"]


def preregistered_sha() -> str:
    git = lambda *a: subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True,
                                    text=True)
    for f in PREREG:
        if git("ls-files", "--error-unmatch", f).returncode:
            sys.exit(f"refusing: {f} is not committed")
    if git("diff", "--quiet", "HEAD", "--", *PREREG).returncode:
        sys.exit("refusing: pre-registered files differ from HEAD")
    return git("rev-parse", "HEAD").stdout.strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--samples", type=Path, default=ROOT / "pipeline/samples/phase5_samples.tsv")
    ap.add_argument("--markers", type=Path, default=ROOT / "pipeline/markers/phase5_markers.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "results/phase6_gse")
    ap.add_argument("--no-guard", action="store_true",
                    help="skip the commit guard (dry runs on stand-in data only)")
    ap.add_argument("--gse-dir", type=Path, default=None,
                    help="GSE186007 pipeline outdir (default <data>/phase6_gse)")
    args = ap.parse_args()
    sha = "" if args.no_guard else preregistered_sha()
    t0 = time.time()
    data = p5s.data_dir() / "phase5"
    gse = args.gse_dir or p5s.data_dir() / "phase6_gse"
    args.out.mkdir(parents=True, exist_ok=True)

    smp = pd.read_csv(args.samples, sep="\t")
    gruns = pd.read_csv(gse / "runs.tsv", sep="\t")
    blocks = p5s.load_blocks(data / "blocks_cpg.tsv")
    off = p5s.load_blocks(data / "offmarker_cpg.tsv")
    # The GSE pipeline run indexed the same blocks.
    for f in ("blocks_cpg.tsv", "offmarker_cpg.tsv"):
        assert (gse / f).read_bytes() == (data / f).read_bytes(), f
    bsets = {s: blocks[blocks.set == s].reset_index(drop=True) for s in SETS}
    nb = {s: len(b) for s, b in bsets.items()}
    ranges = {s: (b.startCpG.to_numpy(), b.endCpG.to_numpy()) for s, b in bsets.items()}
    ranges["off"] = (off.startCpG.to_numpy(), off.endCpG.to_numpy())
    (gse / "records").mkdir(exist_ok=True)
    for r in gruns.run:
        p5s.build_records((str(gse / "pat" / f"{r}.dedup.pat.gz"),
                           str(gse / "records" / f"{r}.dedup.npz"), ranges))
    mk = pd.read_csv(args.markers, sep="\t")
    T = {s: mk[mk.set == s].groupby("block").tumour_mean.mean().reindex(bsets[s].block)
         .to_numpy() for s in SETS}
    thr = smp.run[smp.role == "threshold"].tolist()
    assert len(thr) == 14
    ev = smp[smp.role == "eval"]
    R, counts = {}, {}
    for r in thr + ev.run.tolist():
        R[r], counts[r] = p5s.load_cache(str(data / "records" / f"{r}.dedup.npz"), SETS + ["off"])
    for r in gruns.run:
        R[r], counts[r] = p5s.load_cache(str(gse / "records" / f"{r}.dedup.npz"), SETS + ["off"])
    # Same model as scripts/phase6_score.py: Phase 5 rho and p0 rule, all 14 threshold donors.
    p0_all = {s: p6._p0([p6.summarise(R[d][s], nb[s], 0.0) for d in thr]) for s in SETS}
    rho = {s: p5.fit_rho(p5.Records.concat([R[d][s] for d in thr]), p0_all[s]) for s in SETS}
    S = {r: {s: p6.summarise(R[r][s], nb[s], rho[s]) for s in SETS} for r in R}
    full = {s: p6.fit([S[d][s] for d in thr], T[s]) for s in SETS}
    print(f"loaded and fitted {time.time() - t0:.0f}s; rho {rho}", flush=True)

    def scores(Sx, k):
        return [(name, s, fn(Sx[s], full[s], k)) for s in SETS
                for name, fn in (("WZ", p6.wz), ("FRAG", p6.frag))]

    rows, eff = [], []
    n_sp = {r: p5.n_fragments(R[r]["off"]) for r in gruns.run}
    for k in KS:
        col = f"pool_k{k}"
        for pid in sorted(ev[col].unique()):
            donors = ev.run[ev[col] == pid].tolist()
            assert len(donors) == k
            Sbg = {s: sum(S[d][s] for d in donors) for s in SETS}
            n_bg = p5.n_fragments_pool([R[d]["off"] for d in donors])
            for name, s, v in scores(Sbg, k):
                rows.append(("none", "none", k, pid, 0.0, 0, name, s, v))
            for i, g in enumerate(gruns.itertuples()):
                for fi, f in enumerate(F_GRID, 1):
                    q = f / (1 - f) * n_bg / n_sp[g.run]
                    for seed in range(args.seeds):
                        if q > 1:
                            rows += [(g.run, g.group, k, pid, f, seed, n, s, np.nan)
                                     for s in SETS for n in ("WZ", "FRAG")]
                            continue
                        rng = np.random.default_rng([SEED, 2, k, pid, fi, seed, i])
                        kept = p5.thin_rows(counts[g.run], q, rng)
                        Ssp = {s: p6.summarise(p5.keep_rows(R[g.run][s], kept), nb[s], rho[s])
                               for s in SETS}
                        Sx = {s: Sbg[s] + Ssp[s] for s in SETS}
                        for name, s, v in scores(Sx, k):
                            rows.append((g.run, g.group, k, pid, f, seed, name, s, v))
                        for s in SETS:
                            m = full[s]
                            ok = np.isfinite(m.p0)
                            Csp, _ = p6.frag_counts(Ssp[s], m)
                            eff.append((g.run, g.group, k, pid, f, seed, s,
                                        p6._get(Ssp[s], "N")[ok].sum()
                                        / p6._get(Sx[s], "N")[ok].sum(),
                                        p6._get(Ssp[s], "n3")[m.sel].sum()
                                        / p6._get(Sx[s], "n3")[m.sel].sum(), Csp))
        print(f"k={k} done {time.time() - t0:.0f}s", flush=True)
    cols = ["spike", "group", "k", "pool", "f", "seed", "score", "set", "value"]
    pd.DataFrame(rows, columns=cols).to_csv(args.out / "eval_scores.tsv.gz", sep="\t",
                                            index=False, float_format="%.6g")
    pd.DataFrame(eff, columns=["spike", "group", "k", "pool", "f", "seed", "set",
                               "share_calls", "share_frag_n3_selected",
                               "spike_tumour_like_selected"]).to_csv(
        args.out / "effective_share.tsv", sep="\t", index=False, float_format="%.6g")
    direct = [(g.run, g.group, name, s, v) for g in gruns.itertuples()
              for name, s, v in scores(S[g.run], 1)]
    pd.DataFrame(direct, columns=["run", "group", "score", "set", "value"]).to_csv(
        args.out / "direct_scores.tsv", sep="\t", index=False, float_format="%.6g")
    # Confound check: calls per fragment record, threshold donors vs each GSE group.
    cf = []
    for s in SETS:
        sel = full[s].sel
        who = [("threshold_donors", thr)] + [(g, gruns.run[gruns.group == g].tolist())
                                             for g in sorted(gruns.group.unique())]
        for name, rr in who:
            r = p5.Records.concat([R[x][s] for x in rr])
            for blk, mask in (("all", np.ones(nb[s], bool)), ("frag_selected", sel)):
                kk = mask[r.block]
                n, w = (r.m + r.u)[kk], r.w[kk]
                cf.append((s, name, blk, float((w * n).sum() / w.sum()),
                           float(w[n >= p6.MIN_CALLS].sum() / w.sum())))
    pd.DataFrame(cf, columns=["set", "source", "blocks", "mean_calls_per_record",
                              "share_records_ge3_calls"]).to_csv(
        args.out / "confound_calls.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame([(r, n_sp[r]) for r in gruns.run], columns=["run", "offmarker_fragments"]
                 ).to_csv(args.out / "spike_sizes.tsv", sep="\t", index=False)
    (args.out / "eval_run_sha.txt").write_text(sha + "\n")
    print(f"wrote {len(rows)} rows at {sha[:7]} {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
