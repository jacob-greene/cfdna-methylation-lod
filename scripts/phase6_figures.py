"""Figures for Phases 5 and 6, drawn from the committed summary tables only.

figures/phase6_hct116_detection.png
    Detection rate against spike fraction f (HCT116 cell-line spike into healthy SPOT-MAS
    cfDNA), tumour marker set, pool sizes k = 1 and 2. Phase 5 scores (UXM, LLR_binom) and
    Phase 6 scores (FRAG, WZ), plus the LLR_binom result under an ideal binomial read model.
figures/phase6_gse_detection.png
    Detection rate against f for each GSE186007 plasma run used as the spike, FRAG score,
    tumour marker set, bound rule, k = 1 and 2. Colour marks the tumour fraction (TF).

Usage: uv run python scripts/phase6_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"

plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 150
for side in ("top", "right"):
    plt.rcParams[f"axes.spines.{side}"] = False
plt.rcParams["axes.grid"] = False
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]

PAIRED = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
          "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"]
KS = (1, 2)
F_TICKS = [0.001, 0.01, 0.1]
F_LABELS = ["0.1 %", "1 %", "10 %"]


def _axes_style(ax, k):
    ax.set_xscale("log")
    ax.set_xticks(F_TICKS, F_LABELS)
    ax.set_yticks([0, 1])
    ax.set_ylim(-0.03, 1.03)
    ax.axhline(0.95, color="0.4", lw=0.8, ls=":")
    ax.set_title(f"k = {k} donor{'s' if k > 1 else ''} per background", fontsize=10)
    ax.set_xlabel("Spike fraction f")


def hct116():
    p5 = pd.read_csv(RES / "phase5/detection.tsv", sep="\t")
    p6 = pd.read_csv(RES / "phase6/detection.tsv", sep="\t")
    p5 = p5[p5.set == "tumour"]
    p6 = p6[(p6.set == "tumour") & (p6["rule"] == "bound")]
    lines = [
        ("FRAG (Phase 6)", p6[(p6.rung == "real") & (p6.score == "FRAG")], PAIRED[5], "-"),
        ("UXM (Phase 5)", p5[(p5.rung == "real") & (p5.score == "UXM")], PAIRED[1], "-"),
        ("WZ (Phase 6)", p6[(p6.rung == "real") & (p6.score == "WZ")], PAIRED[3], "-"),
        ("LLR_binom (Phase 5)", p5[(p5.rung == "real") & (p5.score == "LLR_binom")],
         PAIRED[7], "-"),
        ("LLR_binom, ideal binomial reads", p5[(p5.rung == "binom_real")
                                               & (p5.score == "LLR_binom")], PAIRED[7], "--"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
    for ax, k in zip(axes, KS, strict=True):
        for label, d, col, ls in lines:
            d = d[d.k == k].sort_values("f")
            ax.plot(d.f, d.rate, ls, color=col, marker="o", ms=4, lw=1.5, label=label)
        _axes_style(ax, k)
    axes[0].set_ylabel("Share of tests detected")
    axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(FIG / "phase6_hct116_detection.png", bbox_inches="tight")
    plt.close(fig)


def gse():
    det = pd.read_csv(RES / "phase6_gse/detection_run.tsv", sep="\t")
    tf = pd.read_csv(RES / "phase6_gse/tf.tsv", sep="\t").set_index("run").tf
    det = det[(det.score == "FRAG") & (det.set == "tumour") & (det["rule"] == "bound")]
    tf_colours = {}
    for i, run in enumerate(tf[tf > 0].sort_values(ascending=False).index):
        tf_colours[run] = [PAIRED[5], PAIRED[1], PAIRED[7]][i]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
    for ax, k in zip(axes, KS, strict=True):
        dk = det[det.k == k]
        drawn = set()
        for run, d in dk.groupby("spike"):
            d = d.sort_values("f")
            grp = d.group.iloc[0]
            if run in tf_colours:
                ax.plot(d.f, d.rate, "-o", color=tf_colours[run], ms=4, lw=1.6, zorder=3,
                        label=f"{run}, CRC, TF {tf[run]:.2f}")
            else:
                label = "CRC, TF 0" if grp == "CRC" else "Healthy (HEA)"
                ls = "-" if grp == "CRC" else "--"
                ax.plot(d.f, d.rate, ls, color="0.6" if grp == "CRC" else "0.3", lw=1,
                        zorder=2, label=None if label in drawn else label)
                drawn.add(label)
        _axes_style(ax, k)
    axes[0].set_ylabel("Share of tests detected")
    h, lab = axes[1].get_legend_handles_labels()
    rank = {f"{r}, CRC, TF {tf[r]:.2f}": i for i, r in enumerate(tf_colours)}
    order = sorted(range(len(lab)), key=lambda i: (rank.get(lab[i], len(rank)), lab[i]))
    axes[1].legend([h[i] for i in order], [lab[i] for i in order], frameon=False, fontsize=8,
                   loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(FIG / "phase6_gse_detection.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    hct116()
    gse()
