"""Detection figures for the real-read tests, drawn from the committed summary tables only.

The x-axis is "tumour DNA at the markers": for the HCT116 spike, the effective share (the
pre-registered primary scale: the share of informative fragments in the FRAG-selected blocks
that come from the spike); for a patient-plasma spike, the run's effective share times its
tumour fraction (TF). Bands: 95 % donor-bootstrap intervals (scripts/bootstrap_lod.py). The
point left of the axis break is the false-positive rate (FPR) at f = 0, with its 95 %
Clopper-Pearson interval. Vertical line: the LoD; shaded: the grid points around it. Curves
that fail the 5 % FPR target are grey and get no LoD line.

figures/phase6_hct116_detection.png
    Detection rate, HCT116 cell-line spike into healthy SPOT-MAS cfDNA, tumour marker set,
    k = 1 and 2: FRAG, WZ, UXM and LLR on real reads, and LLR under an ideal binomial model.
figures/phase6_gse_detection.png
    Detection rate per GSE186007 plasma run used as the spike, FRAG, tumour marker set, bound
    rule, k = 1 and 2. Coloured: the 3 cancer runs with TF above 0. Runs with TF 0 have no
    tumour scale; they are drawn at the spike's effective share, an upper bound.

Usage: uv run python scripts/bootstrap_lod.py && uv run python scripts/phase6_figures.py
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
INVALID = "0.6"
FPR_MAX = 0.05
KS = (1, 2)
X_ZERO = 4e-5
F_TICKS = [X_ZERO * 1.3, 1e-3, 1e-2]
F_LABELS = ["0\n(FPR)", "0.1 %", "1 %"]
# Patient runs with TF above 0, named by TF rank (legend names, not accessions).
PATIENTS = {"SRR24725517": "Patient A", "SRR24725519": "Patient B", "SRR24725520": "Patient C"}


def _boot():
    lod = pd.read_csv(RES / "bootstrap/lod.tsv", sep="\t")
    det = pd.read_csv(RES / "bootstrap/detection.tsv", sep="\t")
    return lod[lod["rule"] == "bound"], det[det["rule"] == "bound"]


def _axes_style(ax, k):
    ax.set_xscale("log")
    ax.set_xlim(X_ZERO / 1.5, 0.04)
    ax.set_xticks(F_TICKS, F_LABELS)
    ax.minorticks_off()
    ax.axvline(X_ZERO * 2.6, color="0.3", lw=0.6)
    ax.set_yticks([0, 0.5, 0.95, 1], ["0", "0.5", "0.95", ""])
    ax.set_ylim(-0.03, 1.03)
    ax.axhline(0.95, color="0.4", lw=0.8, ls=":")
    ax.set_title(f"k = {k} donor{'s' if k > 1 else ''} per background", fontsize=10)
    ax.set_xlabel("Tumour DNA at the markers")


def _curve(ax, d, r, col, label, ls="-", mk="o", ms=4, i=0, band=0.15, lw=1.5):
    """Detection curve with band, the f = 0 FPR point, and the LoD with its grid bracket."""
    valid = r.fpr <= FPR_MAX
    c = col if valid else INVALID
    s = d[d.f > 0].sort_values("f")
    ax.fill_between(s.x_marker, s.rate_lo95, s.rate_hi95, color=c, alpha=band, lw=0)
    ax.plot(s.x_marker, s.rate, ls, color=c, marker=mk, ms=ms, lw=lw, label=label)
    z = d[d.f == 0].iloc[0]
    ax.errorbar([X_ZERO * (1 + 0.18 * i)], [z.rate],
                yerr=[[z.rate - z.rate_lo95], [z.rate_hi95 - z.rate]], fmt=mk, color=c, ms=3,
                lw=0.8, capsize=1.5)
    if valid and pd.notna(r.lod_marker):
        f = s.f.to_numpy()
        hi = f[f >= r.lod].min()
        lo = f[f < r.lod].max()
        ax.axvspan(s.x_marker[s.f == lo].iloc[0], s.x_marker[s.f == hi].iloc[0], color=c,
                   alpha=0.07, lw=0)
        ax.axvline(r.lod_marker, color=c, lw=0.9, ls="-.")


def hct116():
    lod, det = _boot()
    lines = [  # label, rung, score, colour, line style, marker, marker size
        ("Fragment count (FRAG)", "real", "FRAG", PAIRED[5], "-", "o", 4),
        ("Weighted Z (WZ)", "real", "WZ", PAIRED[3], "-", "s", 6),
        ("Methylated-fragment share (UXM)", "real", "UXM", PAIRED[1], "-", "^", 4),
        ("Binomial log-likelihood ratio (LLR)", "real", "LLR_binom", PAIRED[7], "-", "v", 4),
        ("LLR, ideal binomial reads", "binom_real", "LLR_binom", PAIRED[7], "--", "v", 4),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=True)
    for ax, k in zip(axes, KS, strict=True):
        dk = det[(det.arm == "hct116") & (det.k == k)]
        lk = lod[(lod.arm == "hct116") & (lod.k == k)]
        for i, (label, rung, score, col, ls, mk, ms) in enumerate(lines):
            d = dk[(dk.rung == rung) & (dk.score == score)]
            r = lk[(lk.rung == rung) & (lk.score == score)].iloc[0]
            _curve(ax, d, r, col, label, ls, mk, ms, i, band=0.12)
        _axes_style(ax, k)
    axes[0].set_ylabel("Share of tests detected")
    _legend(axes)
    fig.tight_layout()
    fig.savefig(FIG / "phase6_hct116_detection.png", bbox_inches="tight")
    plt.close(fig)


def _legend(axes):
    """One legend for both panels; grey entries mark curves that fail the FPR target."""
    from matplotlib.lines import Line2D

    h, lab = axes[0].get_legend_handles_labels()
    h2, lab2 = axes[1].get_legend_handles_labels()
    seen = dict(zip(lab, h, strict=True))
    for la, hh in zip(lab2, h2, strict=True):
        if la not in seen or seen[la].get_color() == INVALID:
            seen[la] = hh
    keys = list(seen)
    hs = [seen[k] for k in keys] + [Line2D([], [], color=INVALID, lw=1.5)]
    axes[1].legend(hs, keys + ["Grey: fails the 5 % FPR target"], frameon=False, fontsize=8,
                   loc="center left", bbox_to_anchor=(1.0, 0.5))


def gse():
    lod, det = _boot()
    det, lod = det[det.arm == "patient"], lod[lod.arm == "patient"]
    tf = pd.read_csv(RES / "phase6_gse/tf.tsv", sep="\t").set_index("run")
    colours = dict(zip(PATIENTS, [PAIRED[5], PAIRED[9], PAIRED[1]], strict=True))
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), sharey=True)
    for ax, k in zip(axes, KS, strict=True):
        dk, lk = det[det.k == k], lod[lod.k == k]
        drawn = set()
        for run in sorted(dk.spike.unique()):
            d = dk[dk.spike == run]
            r = lk[lk.spike == run].iloc[0]
            group = tf.loc[run, "group"]
            if run in PATIENTS:
                _curve(ax, d, r, colours[run],
                       f"{PATIENTS[run]}, colorectal cancer, TF {tf.loc[run, 'tf']:.2f}",
                       i=list(PATIENTS).index(run), lw=1.6)
                continue
            label = "Colorectal cancer, TF 0" if group == "CRC" else "Healthy"
            s = d[d.f > 0].sort_values("f")
            ax.plot(s.x_spike, s.rate, "-" if group == "CRC" else "--",
                    color="0.6" if group == "CRC" else "0.3", lw=0.8, zorder=1,
                    label=None if label in drawn else label)
            drawn.add(label)
        _axes_style(ax, k)
    axes[0].set_ylabel("Share of tests detected")
    h, lab = axes[1].get_legend_handles_labels()
    axes[1].legend(h, lab, frameon=False, fontsize=8, loc="center left",
                   bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(FIG / "phase6_gse_detection.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    hct116()
    gse()
