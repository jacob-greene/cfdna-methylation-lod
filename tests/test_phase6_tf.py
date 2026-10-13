import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("phase6_tf", ROOT / "scripts/phase6_tf.py")
tfm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tfm)

NB = 60  # bins per chromosome; centromere at bin 30


def _genome(tf_true, seed=0, n_hea=4):
    rng = np.random.default_rng(seed)
    rows = [(c, b) for c in tfm.CHROMS for b in range(NB)]
    df = pd.DataFrame(rows, columns=["chrom", "bin"])
    df["n_share"] = 0.0
    df["gc"] = rng.uniform(0.35, 0.6, len(df))
    lam = 3000 * np.exp(3 * (df.gc - 0.45))  # GC bias shared by all runs
    arm = np.where(df.bin < NB // 2, "p", "q")
    gain = ((df.chrom == "chr8") & (arm == "q")) | ((df.chrom == "chr20") & (arm == "q"))
    loss = ((df.chrom == "chr18") & (arm == "q")) | ((df.chrom == "chr17") & (arm == "p"))
    runs = []
    for i in range(n_hea):
        df[f"H{i}"] = rng.poisson(lam)
        runs.append((f"H{i}", "HEA"))
    mult = np.where(gain, 1 + tf_true / 2, np.where(loss, 1 - tf_true / 2, 1.0))
    df["C0"] = rng.poisson(lam * mult)
    runs.append(("C0", "CRC"))
    cen = {c: NB // 2 * tfm.BIN for c in tfm.CHROMS}
    return df, pd.DataFrame(runs, columns=["run", "group"]), cen


def test_tf_recovers_single_copy_changes():
    df, runs, cen = _genome(0.2)
    t, arms, info = tfm.estimate(df, runs, cen, seed=[1, 2])
    crc = t[t.run == "C0"].iloc[0]
    assert abs(crc.tf - 0.2) < 0.04
    assert crc.n_altered_arms == 4
    assert crc.tf_lo95 <= crc.tf <= crc.tf_hi95
    # Acrocentric p arms are dropped; 22 q arms + 17 p arms remain.
    assert info["arms"] == 39


def test_tf_healthy_runs_and_no_change_give_zero():
    df, runs, cen = _genome(0.0, seed=3)
    t, _, _ = tfm.estimate(df, runs, cen, seed=[1, 2])
    assert (t.tf == 0).all()
    assert (t.tf_upper_bound_if_none > 0).all()


def test_gc_norm_removes_gc_trend():
    rng = np.random.default_rng(0)
    gc = rng.uniform(0.3, 0.6, 5000)
    x = tfm._gc_norm(1000 * np.exp(4 * gc), gc)
    assert abs(np.corrcoef(gc, x)[0, 1]) < 0.3
