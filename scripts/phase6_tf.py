"""Phase 6 GSE186007 arm: tumour fraction (TF) from arm-level copy number.

A simplified read-depth method in the spirit of ichorCNA (Adalsteinsson et al. 2017,
PMID 29109393): 1 Mb bins, GC correction, a panel of normals, arm-level log2 ratios and a
single-copy model. No methylation call is used.

Rules:
  - Reads: read 1, primary, mapped with mapped mate, not duplicate, MAPQ >= 30, chr1-22.
    Reads that start inside a marker block (blocks_cpg.tsv) are skipped ("off-marker").
  - Bins: 1 Mb; kept if N share < 0.1 and the panel median count > 0.
  - GC correction: count / total, divided by the run's median over 20 equal-count GC bins.
  - Panel of normals: median of the healthy LABS runs; a healthy run uses the other ones.
  - Arms: split at the centromere midpoint (UCSC cytoBand acen); acrocentric p arms
    (13, 14, 15, 21, 22) dropped; an arm needs >= 10 bins.
  - Arm log2 ratio r_a = median of its bins. Noise: s = MAD of neighbouring-bin
    differences / sqrt(2) (robust to copy changes); SE_a = c s / sqrt(n_a). The inflation c
    is the SD of r_a / (s / sqrt(n_a)) over the healthy runs (leave-one-out), all arms.
  - Altered arm: abs(r_a / SE_a) > 4. TF_a = 2 (2^r_a - 1) for gains and 2 (1 - 2^r_a)
    for losses (one copy changed in a diploid genome). TF = median of TF_a over altered
    arms, capped at 1. With no altered arm, TF is reported as below an upper bound
    2 (2^(4 median SE_a) - 1).
  - Interval: 500 bootstrap resamples of bins within arms, rng = default_rng([SEED, 80]);
    2.5 and 97.5 % quantiles of TF (0 for resamples without an altered arm).

Subcommands:
  gc     <genome.fa> <out.tsv>                     GC and N share per 1 Mb bin
  bins   <bam> <out.tsv.gz>                        read counts per bin
  tf     <runs.tsv> <bins_dir> <gc.tsv> <out_dir>  runs.tsv: run, group (CRC or HEA)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import importlib.util

_spec = importlib.util.spec_from_file_location("phase5_score", ROOT / "scripts/phase5_score.py")
p5s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5s)

BIN = 1_000_000
CHROMS = [f"chr{i}" for i in range(1, 23)]
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
Z_ALT, N_BOOT, MIN_BINS = 4.0, 500, 10


def gc(fasta: Path, out: Path) -> None:
    rows = []
    for c in CHROMS:
        seq = subprocess.run(["samtools", "faidx", str(fasta), c], capture_output=True,
                             text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "")
        s = np.frombuffer(seq.upper().encode(), dtype="S1")
        for i in range(0, len(s), BIN):
            w = s[i:i + BIN]
            n = (w == b"N").sum()
            g = ((w == b"G") | (w == b"C")).sum()
            rows.append((c, i // BIN, n / len(w), g / max(len(w) - n, 1)))
        print(c, flush=True)
    pd.DataFrame(rows, columns=["chrom", "bin", "n_share", "gc"]).to_csv(out, sep="\t",
                                                                        index=False)


def bins(bam: Path, out: Path) -> None:
    blk = p5s.load_blocks(p5s.data_dir() / "phase5" / "blocks_cpg.tsv")
    cmd = ["samtools", "view", "-f", "64", "-F", "3852", "-q", "30", str(bam), *CHROMS]
    pr = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    pos = {c: [] for c in CHROMS}
    for line in pr.stdout:
        f = line.split("\t", 4)
        pos[f[2]].append(int(f[3]))
    assert pr.wait() == 0
    rows = []
    for c in CHROMS:
        p = np.asarray(pos[c], dtype=np.int64)
        b = blk[blk.chrom == c].sort_values("start")
        i = np.searchsorted(b.start.to_numpy(), p, side="right") - 1
        inside = (i >= 0) & (p < b.end.to_numpy()[np.maximum(i, 0)])
        cnt = np.bincount(p[~inside] // BIN)
        rows += [(c, j, int(n)) for j, n in enumerate(cnt)]
    pd.DataFrame(rows, columns=["chrom", "bin", "count"]).to_csv(out, sep="\t", index=False)


def _arms(cyto: Path) -> dict:
    cb = pd.read_csv(cyto, sep="\t", header=None, names=["chrom", "s", "e", "band", "stain"])
    cen = cb[cb.stain == "acen"].groupby("chrom").agg(s=("s", "min"), e=("e", "max"))
    return {c: (cen.loc[c, "s"] + cen.loc[c, "e"]) / 2 for c in CHROMS}


def _gc_norm(cnt: np.ndarray, gcv: np.ndarray) -> np.ndarray:
    x = cnt / cnt.sum()
    q = np.quantile(gcv, np.linspace(0, 1, 21))
    gb = np.clip(np.searchsorted(q, gcv, side="right") - 1, 0, 19)
    med = np.array([np.median(x[gb == j]) for j in range(20)])
    return x / np.maximum(med[gb], 1e-12)


def _arm_stats(lr: np.ndarray, arm: np.ndarray, arms: list, idx=None):
    s = np.median(np.abs(np.diff(lr))) * 1.4826 / np.sqrt(2)
    r, n = [], []
    for a in arms:
        v = lr[arm == a] if idx is None else lr[idx[a]]
        r.append(np.median(v))
        n.append(len(v))
    return np.array(r), np.array(n), s


def _tf(r, se):
    z = r / se
    alt = np.abs(z) > Z_ALT
    if not alt.any():
        return 0.0, 0
    t = np.where(r > 0, 2 * (2 ** r - 1), 2 * (1 - 2 ** r))[alt]
    return float(min(np.median(t), 1.0)), int(alt.sum())


def estimate(df: pd.DataFrame, runs: pd.DataFrame, cen: dict, seed=None):
    """TF per run from bin counts. df: chrom, bin, n_share, gc and one count column per run.
    Returns (tf table, arm table, info dict)."""
    hea = runs.run[runs.group == "HEA"].tolist()
    df = df[(df.n_share < 0.1) & (df[hea].median(axis=1) > 0)].reset_index(drop=True)
    df = df.assign(arm=[c[3:] + ("p" if (b + 0.5) * BIN < cen[c] else "q")
                        for c, b in zip(df.chrom, df.bin)])
    df = df[~df.arm.isin([c[3:] + "p" for c in ACRO])]
    n_bins = df.arm.value_counts()
    arms = sorted(n_bins.index[n_bins >= MIN_BINS], key=lambda a: (int(a[:-1]), a[-1]))
    df = df[df.arm.isin(arms)].reset_index(drop=True)
    norm = {r: _gc_norm(df[r].to_numpy(float), df.gc.to_numpy()) for r in runs.run}

    def lr_of(r):
        panel = np.median(np.stack([norm[h] for h in hea if h != r]), axis=0)
        return np.log2(np.maximum(norm[r], 1e-6) / np.maximum(panel, 1e-6))

    LR = {r: lr_of(r) for r in runs.run}
    arm = df.arm.to_numpy()
    # Inflation c from the healthy runs (leave-one-out panels).
    zs = []
    for h in hea:
        r_a, n_a, s = _arm_stats(LR[h], arm, arms)
        zs.append(r_a / (s / np.sqrt(n_a)))
    c_inf = float(np.std(np.concatenate(zs), ddof=1))
    rng = np.random.default_rng([p5s.SEED, 80] if seed is None else seed)
    pos = {a: np.where(arm == a)[0] for a in arms}
    res, arm_rows = [], []
    for r in runs.itertuples():
        r_a, n_a, s = _arm_stats(LR[r.run], arm, arms)
        se = c_inf * s / np.sqrt(n_a)
        t, n_alt = _tf(r_a, se)
        boot = []
        for _ in range(N_BOOT):
            idx = {a: rng.choice(pos[a], size=len(pos[a])) for a in arms}
            rb, nb_, sb = _arm_stats(LR[r.run], arm, arms, idx)
            boot.append(_tf(rb, c_inf * sb / np.sqrt(nb_))[0])
        ub = float(2 * (2 ** (Z_ALT * np.median(se)) - 1))
        res.append((r.run, r.group, t, n_alt, *np.quantile(boot, [0.025, 0.975]), ub, s))
        arm_rows += [(r.run, a, x, x / e) for a, x, e in zip(arms, r_a, se)]
    tf_tab = pd.DataFrame(res, columns=["run", "group", "tf", "n_altered_arms", "tf_lo95",
                                        "tf_hi95", "tf_upper_bound_if_none", "bin_sd"])
    arm_tab = pd.DataFrame(arm_rows, columns=["run", "arm", "log2_ratio", "z"])
    return tf_tab, arm_tab, {"bins": len(df), "arms": len(arms), "inflation_c": c_inf}


def tf(runs_tsv: Path, bins_dir: Path, gc_tsv: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    runs = pd.read_csv(runs_tsv, sep="\t")
    df = pd.read_csv(gc_tsv, sep="\t")
    for r in runs.run:
        c = pd.read_csv(bins_dir / f"{r}.bins.tsv.gz", sep="\t")
        df = df.merge(c.rename(columns={"count": r}), on=["chrom", "bin"], how="left")
    df = df.fillna({r: 0 for r in runs.run})
    cen = _arms(p5s.data_dir() / "ref/annot/hg38.cytoBand.txt.gz")
    tf_tab, arm_tab, info = estimate(df, runs, cen)
    tf_tab.to_csv(out / "tf.tsv", sep="\t", index=False, float_format="%.4g")
    arm_tab.to_csv(out / "tf_arms.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(info.items(), columns=["item", "value"]).to_csv(
        out / "tf_info.tsv", sep="\t", index=False, float_format="%.4g")
    print(tf_tab.to_string(index=False))


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    {"gc": gc, "bins": bins, "tf": tf}[cmd](*map(Path, a))
