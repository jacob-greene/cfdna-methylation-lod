"""Key-findings figure for the README, drawn from the committed summary tables only.

figures/key_findings.png
    (a) Limit of detection (LoD) at each step from idealised to realistic data: trained
        classifiers on arrays (Phase 3b), simulated reads (Phase 3), real reads under an
        ideal binomial model (Phase 5), real reads with the HCT116 cell-line spike
        (Phases 5-6), and real reads with a colorectal-cancer plasma spike (Phase 6).
    (b) Detection rate against spike fraction f at k = 2 donors per background, tumour
        marker set: the ideal binomial model, the HCT116 spike and the one patient run
        with a finite LoD.

Usage: uv run python scripts/key_findings_figure.py
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
F_TOP = 0.1          # top of the tested f grid in the read phases
SIM_TOP = 0.03       # top of the f grid in the depth simulation (Phase 3)
FPR_MAX = 0.05       # FPR target
KS = (1, 2)
PATIENT = "SRR24725519"
X_TICKS = [0.001, 0.003, 0.01, 0.03, 0.1]
X_LABELS = ["0.1 %", "0.3 %", "1 %", "3 %", "10 %"]


def ladder_rows() -> list[tuple[str, str, pd.DataFrame]]:
    """Return (label, colour, frame with columns lod, valid) per row, top to bottom."""
    p3b = pd.read_csv(RES / "phase3b/lod.tsv", sep="\t")
    p3b = p3b[(p3b.track == "array") & (p3b.hi1 == "without") & p3b.specificity_ok]
    p3 = pd.read_csv(RES / "phase3/lod.tsv", sep="\t")
    p3 = p3[(p3.depth >= 100) & (p3["rule"] == "valid") & (p3.hi1 == "without")
            & (p3.bg == "person") & (p3.cap == "none") & p3.specificity_ok]
    p5 = pd.read_csv(RES / "phase5/lod.tsv", sep="\t")
    p5 = p5[(p5.set == "tumour") & p5.k.isin(KS)]
    p6 = pd.read_csv(RES / "phase6/lod.tsv", sep="\t")
    p6 = p6[(p6.set == "tumour") & p6.k.isin(KS) & (p6["rule"] == "bound")]
    binom = p5[(p5.rung == "binom_real") & (p5.score == "LLR_binom")]
    uxm = p5[(p5.rung == "real") & (p5.score == "UXM")]
    p6r = p6[(p6.rung == "real") & p6.score.isin(["FRAG", "WZ"])]
    hct = pd.concat([
        pd.DataFrame({"lod": uxm.lod, "valid": uxm.fpr_max <= FPR_MAX}),
        pd.DataFrame({"lod": p6r.lod, "valid": p6r.fpr <= FPR_MAX}),
    ]).dropna()
    return [
        ("Arrays,\ntrained classifiers", PAIRED[3],
         pd.DataFrame({"lod": p3b.lod_interp, "valid": True}).dropna()),
        ("Simulated reads,\n100+ reads per CpG", PAIRED[1],
         pd.DataFrame({"lod": p3.lod_interp, "valid": True})),
        ("Real reads,\nideal binomial model", PAIRED[7],
         pd.DataFrame({"lod": binom.lod, "valid": binom.fpr_max <= FPR_MAX}).dropna()),
        ("Real reads,\nHCT116 spike", PAIRED[5], hct),
    ]


def ladder(ax):
    rows = ladder_rows()
    n = len(rows) + 1
    ax.axvspan(1e-4, 1e-3, color="0.92", zorder=0)
    ax.text(3.2e-4, n - 0.35, "MRD\nneeds", ha="center", va="top", fontsize=7, color="0.35")
    # The 10 % grid top applies to the real-read rows only (rows 0-2 from the bottom).
    ax.vlines(F_TOP, -0.6, 2.4, color="0.5", lw=0.8, ls="--", zorder=1)
    ax.text(F_TOP * 1.06, 2.45, "top of f grid,\nreal reads", ha="left", va="bottom",
            fontsize=7, color="0.35")
    labels = []
    for label, _, d in rows:
        print(f"{label.replace(chr(10), ' ')}: LoD {d.lod.min():.4f}-{d.lod.max():.4f} "
              f"({int(d.lod.notna().sum())} finite, {int(d.lod.isna().sum())} censored; "
              f"{int(d.valid.sum())} of {len(d)} with FPR <= {FPR_MAX})")
    for i, (label, col, d) in enumerate(rows):
        y = n - 1 - i
        labels.append((y, label))
        # Simulation cells not detected at 95 % by the top of their grid (3 %): censored.
        cens = d[d.lod.isna()]
        for j in range(len(cens)):
            yy = y - 0.3 - 0.07 * j
            ax.annotate("", xy=(SIM_TOP * 1.9, yy), xytext=(SIM_TOP, yy),
                        arrowprops={"arrowstyle": "->", "color": col, "lw": 1.2})
        if len(cens):
            ax.text(SIM_TOP / 1.08, y - 0.37, f"{len(cens)} censored at 3 %", ha="right",
                    va="center", fontsize=7, color=col)
        d = d.dropna()
        jit = [(j % 5 - 2) * 0.06 for j in range(len(d))]
        ok = d.valid.to_numpy()
        ys = [y + j for j in jit]
        ax.scatter(d.lod[ok], [v for v, o in zip(ys, ok, strict=True) if o], s=22,
                   color=col, edgecolor=col, zorder=3)
        ax.scatter(d.lod[~ok], [v for v, o in zip(ys, ok, strict=True) if not o], s=22,
                   facecolor="white", edgecolor=col, zorder=3)
    # Patient plasma: one run brackets its LoD between grid points, the others are censored.
    run = pd.read_csv(RES / "phase6_gse/lod_run.tsv", sep="\t")
    run = run[(run.score == "FRAG") & (run.set == "tumour") & (run["rule"] == "bound")
              & run.k.isin(KS) & (run.tf > 0)]
    y = 0
    labels.append((y, "Real reads,\npatient plasma"))
    col = PAIRED[9]
    finite = run[run.lod.notna()]
    for j, (_, r) in enumerate(finite.iterrows()):
        yy = y + 0.12 * j
        ax.plot([0.03, F_TOP], [yy, yy], color=col, lw=2.5, solid_capstyle="butt", zorder=3)
    n_cens = int(run.lod.isna().sum())
    for j in range(n_cens):
        yy = y - 0.12 * (j + 1) + 0.06
        ax.annotate("", xy=(F_TOP * 1.9, yy), xytext=(F_TOP, yy),
                    arrowprops={"arrowstyle": "->", "color": col, "lw": 1.2})
    ax.text(0.03 / 1.08, y + 0.02, f"{len(finite)} bracketed", ha="right", va="center",
            fontsize=7, color=col)
    ax.text(F_TOP * 2.0, y - 0.25, f"{n_cens} censored", ha="left", va="center", fontsize=7,
            color=col)
    ax.set_yticks([v for v, _ in labels], [lab for _, lab in labels], fontsize=8)
    ax.set_ylim(-0.7, n - 0.2)
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 0.4)
    ax.set_xticks([1e-4] + X_TICKS, ["0.01 %"] + X_LABELS, fontsize=8)
    ax.set_xlabel("LoD (spike fraction f)")
    ax.set_title("a   LoD from idealised to real data", loc="left", fontsize=10)


def curves(ax):
    k = 2
    p5 = pd.read_csv(RES / "phase5/detection.tsv", sep="\t")
    p6 = pd.read_csv(RES / "phase6/detection.tsv", sep="\t")
    gse = pd.read_csv(RES / "phase6_gse/detection_run.tsv", sep="\t")
    lines = [
        ("Ideal binomial model (LLR_binom)",
         p5[(p5.rung == "binom_real") & (p5.score == "LLR_binom") & (p5.set == "tumour")],
         PAIRED[7], "--"),
        ("HCT116 spike, real reads (FRAG)",
         p6[(p6.rung == "real") & (p6.score == "FRAG") & (p6.set == "tumour")
            & (p6["rule"] == "bound")], PAIRED[5], "-"),
        (f"Patient plasma {PATIENT} (FRAG)",
         gse[(gse.spike == PATIENT) & (gse.score == "FRAG") & (gse.set == "tumour")
             & (gse["rule"] == "bound")], PAIRED[9], "-"),
    ]
    for label, d, col, ls in lines:
        d = d[d.k == k].sort_values("f")
        ax.plot(d.f, d.rate, ls, color=col, marker="o", ms=4, lw=1.6, label=label)
    ax.axhline(0.95, color="0.4", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(X_TICKS, X_LABELS, fontsize=8)
    ax.set_yticks([0, 0.5, 0.95, 1], ["0", "0.5", "0.95", ""], fontsize=8)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Spike fraction f")
    ax.set_ylabel("Share of tests detected")
    ax.set_title(f"b   Detection at k = {k} donors per background", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.28, 0.8))


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})
    ladder(axes[0])
    curves(axes[1])
    fig.tight_layout()
    fig.savefig(FIG / "key_findings.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
