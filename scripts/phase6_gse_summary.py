"""Phase 6 GSE186007 arm: detection, LoD per spike run, tumour-fraction scale, LABS control.

Reads results/phase6_gse/{eval_scores.tsv.gz, effective_share.tsv, direct_scores.tsv}
(scripts/phase6_gse_score.py), results/phase6/{thresholds.tsv, eval_scores.tsv.gz} and the
TF table <data>/phase6_gse/tf/tf.tsv (scripts/phase6_tf.py). Rules (amendment 2 of
pipeline/PREREGISTERED_PHASE6.md):
  - thresholds: rung "real" thresholds of results/phase6/thresholds.tsv (bound and model);
  - a (spike, score, set, k, f) cell with any unscored test (NaN) is dropped;
  - LoD per spike run: lowest f detected in >= 95 % of its tests, interpolated on log f;
  - tumour scale: f_tumour = f * TF of the spike run; LoD_tumour = LoD * TF; NaN when the
    run has no altered arm (TF = 0);
  - tumour effective share: LoD mapped through the run's mean share_frag_n3_selected (FRAG)
    or share_calls (WZ) per (k, f), times TF;
  - healthy-LABS control: per (score, set, k, rule, f), the detection rate over all HEA
    spike tests. A CRC cell whose HEA rate exceeds 0.2 is flagged "protocol".

Writes results/phase6_gse/{tf.tsv, fpr.tsv, detection_run.tsv, detection_group.tsv,
lod_run.tsv, lod_summary.tsv, direct.tsv}.

Usage: uv run python scripts/phase6_gse_summary.py [results_dir [tf.tsv]]
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import phase5 as p5
from src import simulate as sim

_spec = importlib.util.spec_from_file_location("phase6_summary", ROOT / "scripts/phase6_summary.py")
p6s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p6s)
_spec = importlib.util.spec_from_file_location("phase5_score", ROOT / "scripts/phase5_score.py")
p5s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5s)

HEA_FLAG = 0.2
SHARE = {"WZ": "share_calls", "FRAG": "share_frag_n3_selected"}
K = ["score", "set", "k", "rule"]


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/phase6_gse"
    ev = pd.read_csv(out / "eval_scores.tsv.gz", sep="\t")
    eff = pd.read_csv(out / "effective_share.tsv", sep="\t")
    tf_path = Path(sys.argv[2]) if len(sys.argv) > 2 else p5s.data_dir() / "phase6_gse/tf/tf.tsv"
    tf = pd.read_csv(tf_path, sep="\t")
    tf.to_csv(out / "tf.tsv", sep="\t", index=False, float_format="%.4g")
    thr = pd.read_csv(ROOT / "results/phase6/thresholds.tsv", sep="\t")
    thr = thr[thr.rung == "real"].drop(columns=["rung", "n_pools"])

    # The f = 0 tests must equal those of the HCT116 arm (same model, same pools).
    p6e = pd.read_csv(ROOT / "results/phase6/eval_scores.tsv.gz", sep="\t")
    p6e = p6e[(p6e.rung == "real") & (p6e.kind == "eval") & (p6e.f == 0)]
    p6e = p6e.assign(pool=p6e.pool.astype(int))  # cancer rows made the column str
    null = ev[ev.spike == "none"]
    chk = null.merge(p6e, on=["k", "pool", "score", "set"], suffixes=("", "_p6"))
    assert len(chk) == len(null) == len(p6e)
    assert np.allclose(chk.value, chk.value_p6, rtol=1e-5, equal_nan=True)

    bad = ev[ev.spike != "none"].groupby(["spike", "score", "set", "k", "f"]).value.apply(
        lambda v: v.isna().any())
    bad = bad[bad].reset_index()[["spike", "score", "set", "k", "f"]]
    ev = ev.merge(bad.assign(_bad=1), how="left", on=["spike", "score", "set", "k", "f"])
    ev = ev[ev._bad.isna()].drop(columns="_bad")
    e = ev.merge(thr, on=["score", "set", "k"])
    e["detected"] = e.value > e.threshold

    n0 = e[e.spike == "none"]
    fp = []
    for key, g in n0.groupby(K):
        x, n = int(g.detected.sum()), len(g)
        fp.append((*key, x / n, x, n, *p6s.clopper_pearson(x, n)))
    fpr = pd.DataFrame(fp, columns=[*K, "fpr", "fp", "n", "fpr_lo95", "fpr_hi95"])
    fpr.to_csv(out / "fpr.tsv", sep="\t", index=False, float_format="%.4g")

    sp = e[e.spike != "none"]
    det = []
    for key, g in sp.groupby(["spike", "group", *K, "f"]):
        neg = n0.merge(pd.DataFrame([key[2:-1]], columns=K), on=K)
        det.append((*key, g.detected.mean(), len(g), g.value.median(),
                    p5.auc(g.value, neg.value)))
    det = pd.DataFrame(det, columns=["spike", "group", *K, "f", "rate", "n", "median_score",
                                     "auc"])
    det.to_csv(out / "detection_run.tsv", sep="\t", index=False, float_format="%.4g")
    grp = (sp.groupby(["group", *K, "f"]).detected.agg(rate="mean", n="size").reset_index())
    hea = grp[grp.group == "HEA"][[*K, "f", "rate"]].rename(columns={"rate": "hea_rate"})
    grp = grp.merge(hea, on=[*K, "f"], how="left")
    grp["protocol_flag"] = grp.hea_rate > HEA_FLAG
    grp.to_csv(out / "detection_group.tsv", sep="\t", index=False, float_format="%.4g")

    mean_eff = eff.groupby(["spike", "set", "k", "f"])[list(SHARE.values())].mean()
    tfm = tf.set_index("run")
    lr = []
    for key, g in det.groupby(["spike", "group", *K]):
        run, group, score, s, k, rule = key
        lod = sim.lod_interp(g.f, g.rate)
        t = float(tfm.loc[run, "tf"])
        e_lod = p6s.lod_effective(lod, mean_eff.loc[(run, s, k)][SHARE[score]])
        tt = t if t > 0 else np.nan
        flag = grp[(grp.group == group) & (grp.score == score) & (grp.set == s) & (grp.k == k)
                   & (grp.rule == rule) & (grp.f <= lod if np.isfinite(lod) else False)]
        lr.append((*key, lod, t, lod * tt, e_lod, e_lod * tt,
                   bool(flag.protocol_flag.any()) if len(flag) else False))
    lod = pd.DataFrame(lr, columns=["spike", "group", *K, "lod", "tf", "lod_tumour",
                                    "lod_effective", "lod_effective_tumour",
                                    "protocol_flag"])
    lod.to_csv(out / "lod_run.tsv", sep="\t", index=False, float_format="%.4g")

    sm = []
    for key, g in lod.groupby(["group", *K]):
        c = g[g.tf > 0] if key[0] == "CRC" else g
        sm.append((*key, len(g), int(np.isfinite(g.lod).sum()), np.nanmedian(g.lod)
                   if np.isfinite(g.lod).any() else np.nan,
                   len(c), int(np.isfinite(c.lod_tumour).sum()),
                   *(np.nanquantile(c.lod_tumour, [0, 0.5, 1])
                     if np.isfinite(c.lod_tumour).any() else [np.nan] * 3),
                   np.nanmedian(c.lod_effective_tumour)
                   if np.isfinite(c.lod_effective_tumour).any() else np.nan))
    sm = pd.DataFrame(sm, columns=["group", *K, "n_runs", "n_lod_finite", "lod_median",
                                   "n_runs_tf", "n_lod_tumour_finite", "lod_tumour_min",
                                   "lod_tumour_median", "lod_tumour_max",
                                   "lod_effective_tumour_median"])
    sm.to_csv(out / "lod_summary.tsv", sep="\t", index=False, float_format="%.4g")

    d = pd.read_csv(out / "direct_scores.tsv", sep="\t").merge(
        thr[thr.k == 1].drop(columns="k"), on=["score", "set"])
    d["detected"] = d.value > d.threshold
    dr = []
    for key, g in d.groupby(["score", "set", "rule"]):
        dr.append((*key, int(g.detected[g.group == "CRC"].sum()), int((g.group == "CRC").sum()),
                   int(g.detected[g.group == "HEA"].sum()), int((g.group == "HEA").sum()),
                   p5.auc(g.value[g.group == "CRC"], g.value[g.group == "HEA"])))
    pd.DataFrame(dr, columns=["score", "set", "rule", "crc_detected", "n_crc", "hea_detected",
                              "n_hea", "auc_crc_vs_hea"]).to_csv(
        out / "direct.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 200)
    print(sm[(sm.rule == "bound") & sm.k.isin([1, 2])].to_string(index=False))


if __name__ == "__main__":
    main()
