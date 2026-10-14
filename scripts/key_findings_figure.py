"""Key-findings figure for the README, drawn from the committed summary tables only.

Every x-axis is "tumour DNA at the markers": the share of the marker signal that comes from
tumour DNA. For the HCT116 spike it is the effective share, the pre-registered primary scale
(the share of informative fragments in the FRAG-selected blocks that come from the spike).
For a patient-plasma spike it is the run's effective share times its tumour fraction (TF).
For arrays and simulated reads, f is the tumour (or colon) DNA share of every CpG, so it
needs no mapping. Nominal f stays in the README table.

figures/key_findings.png
    (a) Limit of detection (LoD) at each step. Top: measured on real reads (HCT116 spike,
        patient-plasma spike). Bottom: idealised references (ideal binomial model, simulated
        reads, trained classifiers on arrays). Bars: 95 % donor-bootstrap grid interval
        (scripts/bootstrap_lod.py).
    (b) Detection rate at k = 1 donor per background, tumour marker set, with donor-bootstrap
        bands, the LoD and its grid bracket, and the false-positive rate (FPR) at f = 0.
    (c) Tumour-like spike fragments per test against tumour DNA at the markers (HCT116,
        k = 1), with the count at the LoD.

Usage: uv run python scripts/bootstrap_lod.py && uv run python scripts/key_findings_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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
INVALID = "0.6"      # colour of curves that fail the FPR target
SIM_TOP = 0.03       # top of the f grid in the depth simulation (Phase 3)
FPR_MAX = 0.05       # FPR target
KS = (1, 2)
X_MIN, X_MAX = 1e-4, 0.4
X_TICKS = [1e-4, 0.001, 0.003, 0.01, 0.03, 0.1]
X_LABELS = ["0.01 %", "0.1 %", "0.3 %", "1 %", "3 %", "10 %"]
X_ZERO = 4e-5        # where the f = 0 (FPR) point sits on the log axis of panel b
X_LABEL = "Tumour DNA at the markers"
# Patient runs with a tumour fraction above 0, named by TF rank (legend names, not accessions).
PATIENTS = {"SRR24725517": "Patient A", "SRR24725519": "Patient B", "SRR24725520": "Patient C"}


def boot() -> tuple[pd.DataFrame, pd.DataFrame]:
    lod = pd.read_csv(RES / "bootstrap/lod.tsv", sep="\t")
    det = pd.read_csv(RES / "bootstrap/detection.tsv", sep="\t")
    return lod[lod["rule"] == "bound"], det[det["rule"] == "bound"]


def idealised_rows() -> list[tuple[str, str, pd.DataFrame]]:
    """(label, colour, frame with lod, valid, lo, hi) for the idealised rows, top to bottom."""
    lod, _ = boot()
    binom = lod[(lod.rung == "binom_real") & (lod.score == "LLR_binom")]
    p3b = pd.read_csv(RES / "phase3b/lod.tsv", sep="\t")
    p3b = p3b[(p3b.track == "array") & (p3b.hi1 == "without") & p3b.specificity_ok]
    p3 = pd.read_csv(RES / "phase3/lod.tsv", sep="\t")
    p3 = p3[(p3.depth >= 100) & (p3["rule"] == "valid") & (p3.hi1 == "without")
            & (p3.bg == "person") & (p3.cap == "none") & p3.specificity_ok]
    return [
        ("Ideal binomial\nmodel, HCT116", PAIRED[7],
         pd.DataFrame({"lod": binom.lod_marker, "valid": binom.fpr <= FPR_MAX,
                       "lo": binom.grid_marker_lo95, "hi": binom.grid_marker_hi95})),
        ("Simulated reads,\n100+ reads per CpG", PAIRED[1],
         pd.DataFrame({"lod": p3.lod_interp, "valid": True, "lo": np.nan, "hi": np.nan})),
        ("Arrays, trained\nclassifiers", PAIRED[3],
         pd.DataFrame({"lod": p3b.lod_interp, "valid": True, "lo": np.nan,
                       "hi": np.nan}).dropna(subset=["lod"])),
    ]


def _points(ax, y, d, col, bars=True):
    jit = [(j % 5 - 2) * 0.07 for j in range(len(d))]
    for (_, r), dy in zip(d.iterrows(), jit, strict=True):
        yy = y + dy
        if bars and np.isfinite(r.lo):
            hi = r.hi if np.isfinite(r.hi) else X_MAX
            ax.plot([r.lo, hi], [yy, yy], color=col if r.valid else INVALID, lw=1, zorder=2)
            if not np.isfinite(r.hi):
                ax.annotate("", xy=(X_MAX, yy), xytext=(X_MAX / 1.6, yy),
                            arrowprops={"arrowstyle": "->", "color": col, "lw": 1})
        ax.scatter([r.lod], [yy], s=22, zorder=3, edgecolor=col if r.valid else INVALID,
                   facecolor=col if r.valid else "white")


def _measured(ax, y, r, col, valid=True, label=True):
    """One real-read LoD: point, 95 % grid interval, or a censoring arrow; k = 1 above k = 2."""
    yy = y + (0.17 if r.k == 1 else -0.17)
    c = col if valid else INVALID
    lo = r.grid_marker_lo95
    if np.isfinite(r.lod_marker):
        hi = r.grid_marker_hi95 if np.isfinite(r.grid_marker_hi95) else X_MAX
        ax.plot([lo, hi], [yy, yy], color=c, lw=1, zorder=2)
        if not np.isfinite(r.grid_marker_hi95):
            ax.annotate("", xy=(X_MAX, yy), xytext=(X_MAX / 1.5, yy),
                        arrowprops={"arrowstyle": "->", "color": c, "lw": 1})
        ax.scatter([r.lod_marker], [yy], s=22, zorder=3, edgecolor=c,
                   facecolor=c if valid else "white")
    else:  # censored: the LoD lies above the top of the grid (mapped to the markers)
        ax.annotate("", xy=(lo * 2.2, yy), xytext=(lo, yy),
                    arrowprops={"arrowstyle": "->", "color": c, "lw": 1.2})
    if label and lo > 3 * X_MIN:
        ax.text(lo / 1.08, yy, f"k = {int(r.k)}", fontsize=6.5, ha="right", va="center", color=c)
    elif label:  # no room left of the interval: label right of it
        hi = r.grid_marker_hi95 if np.isfinite(r.grid_marker_hi95) else X_MAX
        ax.text(hi * 1.08, yy, f"k = {int(r.k)}", fontsize=6.5, ha="left", va="center", color=c)


def ladder(ax):
    lod, _ = boot()
    ax.axvspan(1e-4, 1e-3, color="0.92", zorder=0)
    ax.text(3.2e-4, 7.55, "MRD needs", ha="center", va="bottom", fontsize=7, color="0.35")
    labels = []
    # Measured on real reads: HCT116 (FRAG with intervals; UXM and WZ as points only).
    y = 7
    hct = lod[(lod.arm == "hct116") & (lod.rung == "real")]
    for r in hct[hct.score == "FRAG"].itertuples():
        _measured(ax, y, r, PAIRED[5])
    other = hct[hct.score.isin(["UXM", "WZ"]) & np.isfinite(hct.lod_marker)]
    for r in other.itertuples():
        ok = r.fpr <= FPR_MAX
        ax.scatter([r.lod_marker], [y + (0.17 if r.k == 1 else -0.17)], s=16, zorder=3,
                   marker="D", edgecolor=PAIRED[5] if ok else INVALID,
                   facecolor=PAIRED[5] if ok else "white")
    labels.append((y, "HCT116 spike"))
    pat = lod[(lod.arm == "patient") & lod.spike.isin(PATIENTS)]
    tf = pd.read_csv(RES / "phase6_gse/tf.tsv", sep="\t").set_index("run").tf
    for i, (run, name) in enumerate(PATIENTS.items()):
        y = 6 - i
        for r in pat[pat.spike == run].itertuples():
            _measured(ax, y, r, PAIRED[9])
        labels.append((y, f"{name} plasma\n(TF {tf[run]:.2f})"))
    ax.axhline(3.5, color="0.75", lw=0.8)
    ax.text(X_MAX, 7.5, "Measured on real reads", ha="right", fontsize=7.5, style="italic")
    ax.text(X_MAX, 3.25, "Idealised references", ha="right", fontsize=7.5, style="italic")
    for i, (label, c, d) in enumerate(idealised_rows()):
        y = 3 - i * 0.9
        labels.append((y, label))
        d = d.dropna(subset=["lod"])
        if i == 0:  # ideal binomial model: k = 1 and k = 2, like the measured rows
            for (_, r), k in zip(d.iterrows(), (1, 2), strict=True):
                rr = pd.Series({"k": k, "lod_marker": r.lod, "grid_marker_lo95": r.lo,
                                "grid_marker_hi95": r.hi})
                _measured(ax, y, rr, c, valid=bool(r.valid))
        else:
            _points(ax, y, d, c, bars=False)
    # Simulation cells not detected at 95 % by the top of their grid (3 %): censored.
    p3 = pd.read_csv(RES / "phase3/lod.tsv", sep="\t")
    n_cens = int(((p3.depth >= 100) & (p3["rule"] == "valid") & (p3.hi1 == "without")
                  & (p3.bg == "person") & (p3.cap == "none") & p3.specificity_ok
                  & p3.lod_interp.isna()).sum())
    if n_cens:
        ax.annotate("", xy=(SIM_TOP * 1.9, 1.85), xytext=(SIM_TOP, 1.85),
                    arrowprops={"arrowstyle": "->", "color": PAIRED[1], "lw": 1.2})
        ax.text(SIM_TOP * 2.0, 1.85, f"{n_cens} censored", fontsize=6.5, va="center",
                color=PAIRED[1])
    ax.set_yticks([v for v, _ in labels], [lab for _, lab in labels], fontsize=8)
    ax.set_ylim(0.8, 7.9)
    ax.set_xscale("log")
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xticks(X_TICKS, X_LABELS, fontsize=8)
    ax.set_xlabel(X_LABEL)
    ax.set_title("a   LoD, measured and idealised", loc="left", fontsize=10)


def curves(ax):
    k = 1
    lod, det = boot()
    det = det[det.k == k]
    lines = [
        ("HCT116 spike, fragment count (FRAG)",
         ("hct116", "HCT116", "real", "FRAG"), PAIRED[5], "-"),
        ("Patient B plasma (TF 0.41), FRAG",
         ("patient", "SRR24725519", "real", "FRAG"), PAIRED[9], "-"),
        ("Ideal binomial model (LLR)",
         ("hct116", "HCT116", "binom_real", "LLR_binom"), PAIRED[7], "--"),
    ]
    for i, (label, key, col, ls) in enumerate(lines):
        d = det[(det.arm == key[0]) & (det.spike == key[1]) & (det.rung == key[2])
                & (det.score == key[3])].sort_values("f")
        r = lod[(lod.arm == key[0]) & (lod.spike == key[1]) & (lod.rung == key[2])
                & (lod.score == key[3]) & (lod.k == k)].iloc[0]
        valid = r.fpr <= FPR_MAX
        c = col if valid else INVALID
        if not valid:
            label += ", fails FPR target"
        s = d[d.f > 0]
        ax.fill_between(s.x_marker, s.rate_lo95, s.rate_hi95, color=c, alpha=0.18, lw=0)
        ax.plot(s.x_marker, s.rate, ls, color=c, marker="o", ms=3.5, lw=1.5, label=label)
        z = d[d.f == 0].iloc[0]
        xz = X_ZERO * (1 + 0.25 * i)
        ax.errorbar([xz], [z.rate], yerr=[[z.rate - z.rate_lo95], [z.rate_hi95 - z.rate]],
                    fmt="s", color=c, ms=3.5, lw=1, capsize=2)
        if valid and np.isfinite(r.lod_marker):
            # Grid bracket: the grid points (mapped to the markers) around the LoD.
            ax.axvspan(s.x_marker[s.f == 0.03].iloc[0], s.x_marker[s.f == 0.1].iloc[0],
                       color=c, alpha=0.08, lw=0)
            ax.axvline(r.lod_marker, color=c, lw=1, ls="-.")
            ax.text(r.lod_marker * 1.08, 0.3, f"LoD {100 * r.lod_marker:.1f} %", color=c,
                    fontsize=7, rotation=90, va="bottom")
    ax.axhline(0.95, color="0.4", lw=0.8, ls=":")
    ax.axvline(X_ZERO * 2.3, color="0.3", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(X_ZERO / 1.6, 0.05)
    ax.set_xticks([X_ZERO * 1.25] + X_TICKS[1:-1], ["0\n(FPR)"] + X_LABELS[1:-1], fontsize=8)
    ax.minorticks_off()
    ax.set_yticks([0, 0.5, 0.95, 1], ["0", "0.5", "0.95", ""], fontsize=8)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel("Share of tests detected")
    ax.set_title(f"b   Detection at k = {k} donor per background", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="center", bbox_to_anchor=(0.5, 0.58))


def fragments(ax):
    sf = pd.read_csv(RES / "bootstrap/spike_fragments.tsv", sep="\t")
    sf = sf[sf.k == 1]
    lod, _ = boot()
    r = lod[(lod.arm == "hct116") & (lod.rung == "real") & (lod.score == "FRAG")
            & (lod.k == 1)].iloc[0]
    x, y = np.log10(sf.x_marker.to_numpy()), np.log10(sf["mean"].to_numpy())
    at_lod = 10 ** np.interp(np.log10(r.lod_marker), x, y)
    slope = float(np.exp(np.mean(np.log(sf["mean"] / sf.x_marker))))  # fragments per share
    ax.axvspan(1e-4, 1e-3, color="0.92", zorder=0)
    ax.plot(sf.x_marker, sf["mean"], "-o", color=PAIRED[5], ms=4, lw=1.5,
            label="Mean over 80 tests (16 donors)")
    ax.axhline(at_lod, color="0.4", lw=0.8, ls="--")
    ax.text(X_MIN * 1.2, at_lod * 1.15, f"about {at_lod:.0f} at the LoD", fontsize=7,
            color="0.3")
    ax.text(3.2e-4, 60, "MRD", ha="center", fontsize=7, color="0.35")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(X_MIN, 0.05)
    ax.set_ylim(0.005, 100)
    ax.set_xticks(X_TICKS[:-1], X_LABELS[:-1], fontsize=8)
    ax.set_yticks([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"], fontsize=8)
    ax.minorticks_off()
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel("Tumour-like spike fragments\nper test")
    ax.set_title("c   Tumour fragments per test, k = 1", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    print(f"tumour-like spike fragments at the LoD (k = 1): {at_lod:.1f}; at 0.01 % and "
          f"0.1 % (proportional fit): {slope * 1e-4:.3f} and {slope * 1e-3:.2f}")


def main():
    fig = plt.figure(figsize=(10, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1], hspace=0.3, wspace=0.3)
    ladder(fig.add_subplot(gs[0, :]))
    curves(fig.add_subplot(gs[1, 0]))
    fragments(fig.add_subplot(gs[1, 1]))
    fig.savefig(FIG / "key_findings.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
