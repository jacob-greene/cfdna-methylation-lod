"""Phase 5 follow-up (skeptic findings 1 and 2): effective spike share and a symmetric p0.

NOT pre-registered. Two additions after the main Phase 5 scoring:

  A. Effective spike share. For every eval test the mix of phase5_score.py is rebuilt with
     the same random seeds (identical spike fragments). The share of marker-block CpG calls
     that come from the spike is reported next to the nominal f (the share of off-marker
     fragments).
     The spike share uses n_fragments_pool (per-donor sum). phase5_score.py counted pooled
     fragments with n_fragments on the concatenated pool, which merges repeated pat row ids
     across donors and under-doses the spike at k >= 2 (--original-dosing measures that).
  B. Symmetric p0 at k = 4 and 8. phase5_score.py scores each threshold pool against p0
     from the other 14 - k threshold donors, but eval pools against p0 from all 14. Here
     each eval test gets p0 from a random 14 - k threshold donors as well (draw d per pool,
     d = 0..seeds-1; the spiked test of seed d uses draw d; the f = 0 test is scored with
     every draw). Rung "real_p014" scores the same corrected-dosing tests with p0 from all
     14 threshold donors (the original eval rule), for a side-by-side comparison.
     Thresholds are unchanged (thresholds.tsv, rung "real"). LLR_binom and
     LLR_bb only; UXM and T_* do not use p0. The original rows are kept as they are.

Writes results/phase5/{effective_share,symmetric_scores,symmetric_detection,
symmetric_lod}.tsv(.gz).

Usage: uv run python scripts/phase5_symmetric.py [--seeds 5]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase5_score import F_GRID, KS, SEED, SETS, Scorer, data_dir, load_blocks, load_cache

from src import phase5 as p5
from src import simulate as sim

SYM_KS = [4, 8]
LLRS = ("LLR_binom", "LLR_bb")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--samples", type=Path, default=Path("pipeline/samples/phase5_samples.tsv"))
    ap.add_argument("--markers", type=Path, default=Path("pipeline/markers/phase5_markers.tsv"))
    ap.add_argument("--spike", default="SRR28401050")
    ap.add_argument("--out", type=Path, default=Path("results/phase5"))
    ap.add_argument("--original-dosing", action="store_true",
                    help="spike q from n_fragments of the concatenated pool (the Phase 5 bug); "
                         "writes effective_share_original_dosing.tsv only")
    args = ap.parse_args()
    t0 = time.time()
    data = data_dir() / "phase5"
    smp = pd.read_csv(args.samples, sep="\t")
    blocks = load_blocks(data / "blocks_cpg.tsv")
    bsets = {s: blocks[blocks.set == s].reset_index(drop=True) for s in SETS}
    nb = {s: len(b) for s, b in bsets.items()}
    mk = pd.read_csv(args.markers, sep="\t")
    tumour_mean = {}
    for s in SETS:
        tb = mk[mk.set == s].groupby("block").tumour_mean.mean()
        tumour_mean[s] = tb.reindex(bsets[s].block).to_numpy()

    thr_donors = smp.run[smp.role == "threshold"].tolist()
    ev = smp[smp.role == "eval"]
    assert len(thr_donors) == 14 and not set(thr_donors) & set(ev.run)
    R, counts = {}, {}
    for r in thr_donors + ev.run.tolist() + [args.spike]:
        R[r], counts[r] = load_cache(str(data / "records" / f"{r}.dedup.npz"), SETS + ["off"])
    sums = {r: {s: p5.block_sums(R[r][s], nb[s]) for s in SETS} for r in thr_donors}

    def p0_of(donors):
        out = {}
        for s in SETS:
            M = sum(sums[d][s][0] for d in donors)
            N = sum(sum(sums[d][s]) for d in donors)
            out[s] = np.where(N > 0, M / np.maximum(N, 1), np.nan)
        return out

    p0 = p0_of(thr_donors)
    rho = {s: p5.fit_rho(p5.Records.concat([R[d][s] for d in thr_donors]), p0[s])
           for s in SETS}
    sc = Scorer(nb, tumour_mean, rho, {}, {})
    sp = R[args.spike]

    share, rows = [], []
    for k in KS:
        col = f"pool_k{k}"
        for pid in sorted(ev[col].unique()):
            donors = ev.run[ev[col] == pid].tolist()
            bg = {s: p5.Records.concat([R[d][s] for d in donors]) for s in SETS + ["off"]}
            draws = []
            if k in SYM_KS:
                for d in range(args.seeds):
                    rngd = np.random.default_rng([SEED, 9, k, pid, d])
                    draws.append(p0_of(rngd.choice(thr_donors, 14 - k, replace=False)))
            for fi, f in enumerate([0.0, *F_GRID]):
                for seed in range(args.seeds if f else 1):
                    # Same seed as phase5_score.py step 4, dedup: identical spike fragments.
                    rng = np.random.default_rng([SEED, 1, k, pid, fi, seed, 0])
                    if f == 0:
                        recs = bg
                    else:
                        n_bg = (p5.n_fragments(bg["off"]) if args.original_dosing else
                                p5.n_fragments_pool([R[d]["off"] for d in donors]))
                        q = f / (1 - f) * n_bg / p5.n_fragments(sp["off"])
                        kept = p5.thin_rows(counts[args.spike], q, rng)
                        spk = {s: p5.keep_rows(sp[s], kept) for s in bg}
                        recs = {s: p5.Records.concat([bg[s], spk[s]]) for s in bg}
                        for s in SETS + ["off"]:
                            n_sp = float(spk[s].w.sum()) if s == "off" else float(
                                sum(x.sum() for x in p5.block_sums(spk[s], nb[s])))
                            n_bg = float(bg[s].w.sum()) if s == "off" else float(
                                sum(x.sum() for x in p5.block_sums(bg[s], nb[s])))
                            share.append((k, pid, f, seed, s, n_sp / (n_sp + n_bg)))
                    if k not in SYM_KS:
                        continue
                    for name, s, fa, val in sc.all(recs, p0):
                        if name in LLRS:
                            rows.append(("real_p014", "eval", k, pid, f, seed, name, s, fa, val))
                    for d in (range(args.seeds) if f == 0 else [seed]):
                        for name, s, fa, val in sc.all(recs, draws[d]):
                            if name in LLRS:
                                rows.append(("real_sym", "eval", k, pid, f, d, name, s, fa, val))
        print(f"k={k} done {time.time() - t0:.0f}s", flush=True)

    out = args.out
    S = pd.DataFrame(share, columns=["k", "pool", "f", "seed", "set", "spike_share"])
    eff = (S.groupby(["set", "k", "f"]).spike_share.agg(["mean", "min", "max"]).reset_index())
    eff.to_csv(out / ("effective_share_original_dosing.tsv" if args.original_dosing else
                      "effective_share.tsv"), sep="\t", index=False, float_format="%.4g")
    print(eff[eff.f == 0.1].to_string(index=False))
    if args.original_dosing:
        return

    cols = ["rung", "kind", "k", "pool", "f", "seed", "score", "set", "f_alt", "value"]
    d = pd.DataFrame(rows, columns=cols)
    d.to_csv(out / "symmetric_scores.tsv.gz", sep="\t", index=False)
    thr = pd.read_csv(out / "thresholds.tsv", sep="\t")
    thr = thr[thr.rung == "real"].drop(columns="rung")
    e = d.merge(thr, on=["score", "set", "f_alt", "k"], how="left")
    KEY = ["rung", "score", "set", "k"]
    assert e.threshold.notna().all()
    e["detected"] = e.value > e.threshold
    fpr = e[e.f == 0].groupby([*KEY, "f_alt"]).detected.agg(["mean", "size"])
    spk = e[(e.f > 0) & np.isclose(e.f_alt, e.f)]
    det_rows = []
    for (rung, score, s, k, f), g in spk.groupby([*KEY, "f"]):
        neg = e[(e.f == 0) & (e.rung == rung) & (e.score == score) & (e.set == s)
                & (e.k == k) & np.isclose(e.f_alt, f)].value
        det_rows.append((rung, score, s, k, f, g.detected.mean(), len(g),
                         p5.auc(g.value, neg)))
    det = pd.DataFrame(det_rows, columns=[*KEY, "f", "rate", "n", "auc"])
    det.to_csv(out / "symmetric_detection.tsv", sep="\t", index=False, float_format="%.4g")
    lod_rows = []
    for (rung, score, s, k), g in det.groupby(KEY):
        g = g.sort_values("f")
        f0 = fpr.loc[(rung, score, s, k)]
        lod_rows.append((rung, score, s, k, sim.lod_interp(g.f, g.rate),
                         float(f0["mean"].max()), int(f0["size"].max()), int(g.n.min())))
    lod = pd.DataFrame(lod_rows, columns=[*KEY, "lod", "fpr_max", "n_null", "n_per_f"])
    lod.to_csv(out / "symmetric_lod.tsv", sep="\t", index=False, float_format="%.4g")
    print(lod.to_string(index=False))
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
