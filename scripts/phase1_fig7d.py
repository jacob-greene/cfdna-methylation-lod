"""Phase 1: reproduce Moss et al. 2018 Fig. 7d from GEO processed betas (GSE122126).

Steps
1. Per-sample QC table (missing-beta rate, overall and at atlas CpGs).
2. NNLS deconvolution of every EPIC array against the deduplicated reference atlas.
3. Cross-check against the authors' own `Deconvolve.decon_single_samp` (meth_atlas repo).
4. Healthy baseline: the 12 blanks against Moss Supp. Data 1, Table 4.
5. Fig. 7d: 8 points (5 mixes + the 3 background plasmas at 0 %), as in Supp. Data 1, Table 8.
   The paper reports r = 0.92, p < 1.2e-3 (Results, Fig. 7d paragraph). With n = 8 a Pearson
   r of 0.92 gives p = 1.2e-3, so the paper's r uses all 8 points.

Acceptance tolerance (fixed before the run): |r - 0.92| <= 0.05 for observed colon fraction
vs nominal cancer-plasma percentage, n = 8.

Usage: uv run python scripts/phase1_fig7d.py [--data DIR] [--out results/phase1]
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import (
    deconv,
    geo,
)

SEED = 20261008

BLANKS = {  # GSM -> label; order = rows of Moss Supp. Data 1, Table 4
    "GSM3455776": "Young M1", "GSM3455777": "Young M2", "GSM3455778": "Young F1",
    "GSM3455779": "Young F2", "GSM3455780": "Old M1", "GSM3455781": "Old M2",
    "GSM3455782": "Old F1", "GSM3455783": "Old F2", "GSM3455807": "HI1",
    "GSM3455840": "HI2", "GSM3455846": "HI3", "GSM3455848": "HI4",
}
PL806 = ["GSM3455808", "GSM3455810"]
# Fig. 7d points: GEO GSM, Table 8 mix number, background GSM, cancer-plasma fraction p.
FIG7D = [
    ("GSM3455837", 1, "GSM3455840", 0.10),  # GEO cfDNA_mix_5
    ("GSM3455838", 2, "GSM3455840", 0.05),  # GEO cfDNA_mix_6
    ("GSM3455839", 3, "GSM3455840", 0.03),  # GEO cfDNA_mix_7
    ("GSM3455840", 4, "GSM3455840", 0.00),  # HI2 itself
    ("GSM3455845", 5, "GSM3455846", 0.05),  # GEO cfDNA_mix_1
    ("GSM3455846", 6, "GSM3455846", 0.00),  # HI3 itself
    ("GSM3455847", 7, "GSM3455848", 0.05),  # GEO cfDNA_mix_3
    ("GSM3455848", 8, "GSM3455848", 0.00),  # HI4 itself
]
COLON = "Colon_epithelial_cells"
# Fig. 7d "% Predicted" read from the published panel (Fig7_HTML.png, 2054 x 1171 px): orange
# marker centroids, y axis 0-16 % over pixel rows 892.5-499 (24.59 px per %), x axis checked
# against the 0/3/5/10 % positions (29.7 px per %). Accuracy about +-0.1 percentage point.
# The three 0 % points overlap as one marker; the 5 % points are not labelled by background,
# so they are matched to ours by rank.
PAPER_FIG7D_PCT = {0.10: [13.24], 0.05: [10.78, 7.73, 6.56], 0.03: [2.55], 0.00: [1.27]}


def default_data() -> Path:
    env = os.environ.get("CFDNA_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def load_epic(data: Path) -> pd.DataFrame:
    cache = data / "geo" / "betas_epic.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    return geo.read_betas(data, "GPL21145")


def load_authors_code(data: Path):
    """Import the authors' deconvolve.py and load the atlas with their loader."""
    spec = importlib.util.spec_from_file_location("meth_deconv", data / "meth_atlas" / "deconvolve.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    atl = mod.Deconvolve.load_atlas(str(data / "meth_atlas" / "reference_atlas.csv"))
    return atl, mod


def pearson(x, y) -> dict:
    r, p = stats.pearsonr(x, y)
    fit = stats.linregress(x, y)
    return {"n": len(x), "r": r, "p": p, "slope": fit.slope, "intercept": fit.intercept}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=default_data())
    ap.add_argument("--out", type=Path, default=Path("results/phase1"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    np.random.seed(SEED)

    samples = geo.read_samples(args.data)
    betas = load_epic(args.data)
    raw_atlas = atlas_mod.read_atlas(args.data)
    atl = deconv.dedup_atlas(raw_atlas)
    print(f"atlas rows {len(raw_atlas)}, unique CpGs {len(atl)}, cell types {atl.shape[1]}")

    # 1. QC table
    on_atlas = betas.index.isin(atl.index)
    qc = pd.DataFrame({
        "title": samples.loc[betas.columns, "title"],
        "missing_rate_all": betas.isna().mean(),
        "missing_rate_atlas": betas[on_atlas].isna().mean(),
        "atlas_cpgs_present": betas[on_atlas].notna().sum(),
    })
    qc.index.name = "gsm"

    # 2. Deconvolution
    fracs, info = deconv.deconvolve(betas, atl)
    qc = qc.join(info)
    qc.to_csv(args.out / "qc.tsv", sep="\t", float_format="%.5g")
    fracs.index.name = "gsm"
    fracs.insert(0, "title", qc["title"])
    fracs.to_csv(args.out / "deconv_epic.tsv", sep="\t", float_format="%.6g")
    fracs = fracs.drop(columns="title")

    # 3. Authors' code on the Fig. 7d samples + blanks
    a_atl, mod = load_authors_code(args.data)
    rows = []
    for g in sorted(set(BLANKS) | set(PL806) | {f[0] for f in FIG7D}):
        samp = betas[[g]].reset_index().rename(columns={"probe": "acc"})
        mix, _ = mod.Deconvolve.decon_single_samp(samp, a_atl)
        rows.append({"gsm": g, "max_abs_diff": float(np.max(np.abs(mix - fracs.loc[g].to_numpy())))})
    xc = pd.DataFrame(rows).set_index("gsm")
    xc.to_csv(args.out / "authors_code_crosscheck.tsv", sep="\t", float_format="%.3g")
    print(f"authors' code vs ours: max |diff| over {len(xc)} samples = {xc.max_abs_diff.max():.2e}")

    # 4. Blanks vs Supp. Table 4
    xls = args.data / "moss_supp" / "41467_2018_7466_MOESM3_ESM.xlsx"
    t4 = pd.read_excel(xls, "Table 4", header=0, index_col=0)
    t4 = t4[~t4.index.astype(str).str.contains(r"\(SD\)|Mean|SD \(%\)", regex=True)]
    t4 = t4.drop(columns="GE/ml")
    t4.index = list(BLANKS.values())  # 8 pools, then HI1-HI4 (sexes M, M, F, M match GEO)
    ours = deconv.to_table4_groups(fracs.loc[list(BLANKS)])
    ours.index = list(BLANKS.values())
    rename = {c: c.replace("_", " ") for c in ours.columns}
    ours = ours.rename(columns=rename)
    ours = ours.rename(columns={"Cortical neurons": "Neurons", "Lung cells": "Lung",
                                "Vascular endothelial cells": "Vascular endothelial cells"})
    t4.columns = [str(c).strip() for c in t4.columns]
    common = [c for c in t4.columns if c in ours.columns]
    missing = sorted(set(t4.columns) - set(common))
    if missing:
        print("Table 4 columns not matched:", missing)
    long = []
    for s in t4.index:
        for c in common:
            long.append({"blank": s, "cell_type": c, "paper": float(t4.loc[s, c]),
                         "ours": float(ours.loc[s, c])})
    cmp = pd.DataFrame(long)
    cmp["diff"] = cmp.ours - cmp.paper
    cmp.to_csv(args.out / "blanks_vs_table4_long.tsv", sep="\t", float_format="%.5g")
    summ = cmp.groupby("cell_type").agg(paper_mean=("paper", "mean"), ours_mean=("ours", "mean"),
                                        mean_abs_diff=("diff", lambda d: d.abs().mean()),
                                        max_abs_diff=("diff", lambda d: d.abs().max()))
    summ = summ.sort_values("paper_mean", ascending=False)
    summ.to_csv(args.out / "blanks_vs_table4_summary.tsv", sep="\t", float_format="%.4f")
    print(summ.head(10).to_string(float_format="%.4f"))
    wbc = ours[["Monocytes EPIC", "Lymphocytes", "Granulocytes"]].sum(axis=1)
    base = pd.DataFrame({"WBC": wbc, "Erythrocyte progenitors": ours["Erythrocyte progenitors"],
                         "Vascular endothelium": ours["Vascular endothelial cells"],
                         "Hepatocytes": ours["Hepatocytes"], "Colon": ours["Colon epithelial cells"],
                         "Neurons": ours["Neurons"]})
    base.to_csv(args.out / "blanks_baseline.tsv", sep="\t", float_format="%.4f")
    print(base.describe().loc[["mean", "std", "min", "max"]].to_string(float_format="%.4f"))
    cell_r = stats.pearsonr(cmp.paper, cmp.ours)
    print(f"blanks: all cells x types, r(ours, Table 4) = {cell_r[0]:.4f}, n = {len(cmp)}")

    # 5. Fig. 7d
    c806 = fracs.loc[PL806, COLON]
    c806_mean = float(c806.mean())
    pts = []
    for g, mixno, bg, p in FIG7D:
        c_bg = float(fracs.loc[bg, COLON])
        pts.append({"gsm": g, "title": qc.loc[g, "title"], "table8_mix": mixno, "background": bg,
                    "p_nominal": p, "c_background": c_bg,
                    "expected": p * c806_mean + (1 - p) * c_bg,
                    "observed": float(fracs.loc[g, COLON])})
    pts = pd.DataFrame(pts)
    pts.to_csv(args.out / "fig7d_points.tsv", sep="\t", index=False, float_format="%.5f")
    print(f"c(pl806) = {c806.round(4).to_dict()}, mean {c806_mean:.4f}")
    print(pts.to_string(index=False, float_format="%.4f"))
    res = []
    for label, sub in [("all 8", pts), ("without mix_7", pts[pts.gsm != "GSM3455839"]),
                       ("5 mixes only", pts[pts.p_nominal > 0]),
                       ("5 mixes only, without mix_7", pts[(pts.p_nominal > 0) & (pts.gsm != "GSM3455839")])]:
        for xname in ["p_nominal", "expected"]:
            d = pearson(sub[xname], sub["observed"])
            res.append({"set": label, "x": xname, **d})
    res = pd.DataFrame(res)
    res.to_csv(args.out / "fig7d_stats.tsv", sep="\t", index=False, float_format="%.4g")
    print(res.to_string(index=False, float_format="%.4g"))

    # 6. In-silico control: mix the arrays' betas at the nominal p and deconvolve the result.
    # Beta mixing is linear in molecules, so this isolates the algorithm from the wet-lab mix.
    pl806_beta = betas[PL806].mean(axis=1)
    sil = []
    for g, mixno, bg, p in FIG7D:
        if p == 0:
            continue
        for q in [p, 2 * p]:
            vec = q * pl806_beta + (1 - q) * betas[bg]
            f, _, _ = deconv.nnls_one(vec, atl)
            sil.append({"gsm": g, "title": qc.loc[g, "title"], "p_in_silico": q,
                        "expected": q * c806_mean + (1 - q) * float(fracs.loc[bg, COLON]),
                        "in_silico_colon": float(f[COLON]),
                        "observed_real_mix": float(fracs.loc[g, COLON])})
    sil = pd.DataFrame(sil)
    sil.to_csv(args.out / "fig7d_in_silico_control.tsv", sep="\t", index=False, float_format="%.5f")
    print(sil.to_string(index=False, float_format="%.4f"))
    # Implied cancer-plasma fraction of each real mix: solve observed = q*c806 + (1-q)*c_bg.
    pts["implied_p"] = (pts.observed - pts.c_background) / (c806_mean - pts.c_background)
    pts.to_csv(args.out / "fig7d_points.tsv", sep="\t", index=False, float_format="%.5f")
    print(pts[["title", "p_nominal", "observed", "implied_p"]].to_string(index=False, float_format="%.4f"))

    # Compare with the digitized paper points (5 % points matched by rank).
    cmp7 = []
    for p, vals in PAPER_FIG7D_PCT.items():
        ours_p = sorted(pts.loc[pts.p_nominal == p, "observed"] * 100, reverse=True)
        for i, o in enumerate(ours_p):
            pv = vals[min(i, len(vals) - 1)]
            cmp7.append({"p_nominal": p, "ours_pct": o, "paper_pct": pv, "diff_pct": o - pv})
    cmp7 = pd.DataFrame(cmp7)
    cmp7.to_csv(args.out / "fig7d_vs_paper.tsv", sep="\t", index=False, float_format="%.3f")
    print(cmp7.to_string(index=False, float_format="%.3f"))

    # Figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, xname, xl in [(axes[0], "p_nominal", "Cancer plasma in mix (Supp. Table 8)"),
                          (axes[1], "expected", "Expected colon fraction")]:
        for bg, mk in zip(["GSM3455840", "GSM3455846", "GSM3455848"], ["o", "s", "^"]):
            sub = pts[pts.background == bg]
            ax.scatter(sub[xname], sub.observed, marker=mk, label=f"background {BLANKS[bg]}")
        m7 = pts[pts.gsm == "GSM3455839"]
        ax.scatter(m7[xname], m7.observed, facecolors="none", edgecolors="red", s=160, label="mix_7")
        lim = max(pts[xname].max(), pts.observed.max()) * 1.1
        if xname == "expected":
            ax.plot([0, lim], [0, lim], "k:", lw=1)
        ax.set_xlabel(xl)
        ax.set_ylabel("Observed colon fraction (NNLS)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    figdir = Path("figures")
    figdir.mkdir(exist_ok=True)
    fig.savefig(figdir / "phase1_fig7d.png", dpi=150)


if __name__ == "__main__":
    main()
