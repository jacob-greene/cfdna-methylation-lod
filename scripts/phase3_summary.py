"""Phase 3 summary: detection rates, false-positive rates and LoD from the raw simulation runs.

Reads results/phase3/raw/*.tsv.gz (scripts/phase3_depth_sim.py) and writes:
  detection.tsv     detection rate per (kind, depth, cap, background mode, HI1, markers, score, f)
  fpr.tsv           false-positive rate of the held-out blanks at f = 0, per configuration
  fpr_by_blank.tsv  the same, per held-out blank (which people cause false positives)
  lod.tsv           LoD per configuration: lowest grid f at >= 95 % detection, and log-interpolated
  figures/phase3_depth_lod.png

Threshold rules (column "rule"): "pooled" = EP17 0.95 rank of in-sample training null scores;
"valid" = normal prediction bound over inner-leave-one-out training null scores. The figure
shows the valid rule where it exists.

Rates pool held-out blanks x spike replicates x simulation replicates. An LoD is flagged
"specificity fails" when the held-out false-positive rate exceeds 5 % (EP17 LoB rule).

Usage: uv run python scripts/phase3_summary.py [--raw results/phase3/raw] [--out results/phase3]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import simulate as sim

KEYS = ["kind", "depth", "cap", "bg", "hi1", "markers", "score", "rule"]
LABEL = {"GSM3455776": "Young M1", "GSM3455777": "Young M2", "GSM3455778": "Young F1",
         "GSM3455779": "Young F2", "GSM3455780": "Old M1", "GSM3455781": "Old M2",
         "GSM3455782": "Old F1", "GSM3455783": "Old F2", "GSM3455807": "HI1",
         "GSM3455840": "HI2", "GSM3455846": "HI3", "GSM3455848": "HI4"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=Path("results/phase3/raw"))
    ap.add_argument("--out", type=Path, default=Path("results/phase3"))
    ap.add_argument("--fig", type=Path, default=Path("figures/phase3_depth_lod.png"))
    args = ap.parse_args()
    raw = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(args.raw.glob("*.tsv.gz"))],
                    ignore_index=True)
    raw["markers"] = raw["markers"].astype(str)
    # Runs before the "valid" rule existed carry only the pooled rule.
    raw["rule"] = raw["rule"].fillna("pooled") if "rule" in raw else "pooled"
    raw["blank"] = raw.eval_blank.map(LABEL)

    # False positives: NNLS rows at f = 0; LLR_null rows (one threshold per alternative f).
    null = raw[((raw.score == "NNLS") & (raw.f == 0)) | (raw.score == "LLR_null")].copy()
    null["score"] = null.score.replace({"LLR_null": "LLR"})
    null["alt_f"] = np.where(null.score == "LLR", null.f, np.nan)
    fpr_blank = (null.groupby([*KEYS, "blank"]).agg(n=("n", "sum"), fp=("n_detected", "sum"))
                 .reset_index())
    fpr_blank["fpr"] = fpr_blank.fp / fpr_blank.n
    fpr_blank.to_csv(args.out / "fpr_by_blank.tsv", sep="\t", index=False)
    fpr = null.groupby(KEYS).agg(n=("n", "sum"), fp=("n_detected", "sum")).reset_index()
    fpr["fpr"] = fpr.fp / fpr.n
    worst = (fpr_blank.sort_values("fpr").groupby(KEYS).tail(1)
             .set_index(KEYS)[["blank", "fpr"]].rename(columns={"blank": "worst_blank",
                                                                "fpr": "worst_blank_fpr"}))
    fpr = fpr.join(worst, on=KEYS)
    fpr.to_csv(args.out / "fpr.tsv", sep="\t", index=False)

    det = raw[raw.score.isin(["NNLS", "LLR"]) & (raw.f > 0)]
    det = det.groupby([*KEYS, "f"]).agg(n=("n", "sum"), detected=("n_detected", "sum"),
                                        median_score=("median_score", "median")).reset_index()
    det["rate"] = det.detected / det.n
    det.to_csv(args.out / "detection.tsv", sep="\t", index=False)

    rows = []
    for key, g in det.groupby(KEYS):
        r = dict(zip(KEYS, key))
        r["lod_grid"] = sim.lod_from_rates(g.f, g.rate)
        r["lod_interp"] = sim.lod_interp(g.f, g.rate)
        r["rate_at_1e-3"] = float(g.loc[np.isclose(g.f, 1e-3), "rate"].iloc[0])
        r["rate_at_1e-2"] = float(g.loc[np.isclose(g.f, 1e-2), "rate"].iloc[0])
        rows.append(r)
    lod_tab = pd.DataFrame(rows).merge(fpr[[*KEYS, "fpr", "worst_blank", "worst_blank_fpr"]],
                                       on=KEYS)
    lod_tab["specificity_ok"] = lod_tab.fpr <= 0.05
    lod_tab = lod_tab.sort_values(KEYS)
    lod_tab.to_csv(args.out / "lod.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        print(lod_tab.to_string(index=False))

    # Figure: LoD (interpolated) vs depth. Rows: HI1 with / without. Columns: plasma spike,
    # cell spike, plasma spike with the input cap. Open markers: specificity fails. Dotted /
    # "control": every sample uses the training-blank mean (no between-person variation).
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), sharey=True)
    colors = {"100": "#1f78b4", "500": "#33a02c", "all": "#e31a1c"}
    styles = {"NNLS": "--", "LLR": "-"}
    for i, hi1 in enumerate(["with", "without"]):
        for j, (title, sel) in enumerate([
            ("pl806 plasma spike", {"kind": "plasma", "cap": "none"}),
            ("colon-cell spike", {"kind": "cell", "cap": "none"}),
            ("pl806 spike, LLR all, input cap (pooled rule)", {}),
        ]):
            ax = axes[i, j]
            valid = set(lod_tab[lod_tab.rule == "valid"].kind)
            d = lod_tab[(lod_tab.hi1 == hi1) & (lod_tab.bg == "person")]
            d = d[(d.rule == "valid") | ~d.kind.isin(valid) | (d.cap != "none")]
            for k, v in sel.items():
                d = d[d[k] == v]
            if j < 2:
                for (m, s), g in d.groupby(["markers", "score"]):
                    g = g.sort_values("depth")
                    ax.plot(g.depth, g.lod_interp, styles[s], color=colors[m],
                            label=f"{s}, {m} markers")
                    bad = g[~g.specificity_ok]
                    ax.scatter(bad.depth, bad.lod_interp, facecolors="white",
                               edgecolors=colors[m], zorder=3)
                if j == 0:
                    c = lod_tab[(lod_tab.hi1 == hi1) & (lod_tab.bg == "mean")
                                & (lod_tab.kind == "plasma") & (lod_tab.markers == "all")
                                & (lod_tab.score == "LLR") & (lod_tab.cap == "none")
                                ].sort_values("depth")
                    ax.plot(c.depth, c.lod_interp, ":", color="black",
                            label="LLR, all, no between-person variation (control)")
            else:
                # Capped and control runs carry only the pooled rule; keep one rule per line.
                d = lod_tab[(lod_tab.hi1 == hi1) & (lod_tab.kind == "plasma")
                            & (lod_tab.score == "LLR") & (lod_tab.markers == "all")
                            & (lod_tab.rule == "pooled")]
                for (bg, cap), g in d.groupby(["bg", "cap"]):
                    g = g.sort_values("depth")
                    ax.plot(g.depth, g.lod_interp, "-o" if bg == "mean" else "--x",
                            label=f"{'control' if bg == 'mean' else 'real blanks'}, cap {cap}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim(5e-5, 5e-2)
            ax.set_title(f"{title}; HI1 {hi1}")
            ax.set_xlabel("sequencing depth D (reads per CpG)")
            if j == 0:
                ax.set_ylabel("LoD, tumour fraction f (95 % detection)")
            ax.legend(fontsize=7)
    fig.tight_layout()
    args.fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.fig, dpi=150)
    print("wrote", args.fig)


if __name__ == "__main__":
    main()
