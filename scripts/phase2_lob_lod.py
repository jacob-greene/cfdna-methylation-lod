"""Phase 2: limit of blank (LoB) and limit of detection (LoD) of the array colon fraction.

Measurand: the NNLS colon-epithelium fraction of one EPIC array (Moss et al. 2018 atlas).

LoB (CLSI EP17-A2, nonparametric): the 0.95 rank of the blank colon fractions, with a
percentile bootstrap CI. Blank sets:
  - all 12 healthy blanks (8 pools + HI1-HI4), and the 8-unit version with rep1/rep2 pool
    pairs averaged ("collapsed");
  - default split: LoB from the 8 pools + HI1; LoD evaluated on mixtures over HI2-HI4;
  - reverse split: LoB from HI2-HI4; LoD evaluated on mixtures over the 8 pools + HI1.
LoB-setting blanks are never mixture backgrounds in the same split.

LoD:
  (a) parametric, EP17: LoB + c_p * pooled SD of low-level mixtures (f = 0.1-3 % colon);
      translated from colon-fraction units to spike fraction f by a linear calibration
      fitted on the same split's mixtures;
  (b) empirical: the lowest grid level detected (colon fraction > LoB) in >= 95 % of
      replicates, and at every higher level.

In-silico mixes at f in {0, 0.1, 0.3, 1, 3, 10} % colon, plus 1.5 and 2 % for resolution:
  Option 1, cell-type spike: a colon-epithelium replicate (GSM3455852/63/64) mixed into a
    blank's betas. The atlas colon column is replaced by the mean of the other two
    replicates (leave-one-replicate-out). The atlas colon column is the mean of the three
    replicates (checked below), so this removes the spiked replicate from the reference.
    LoB is recomputed with the same leave-one-out atlas. Unlike the paper's simulations,
    the CpG set is not re-selected (Methods, "Simulations"), so this is slightly optimistic.
  Option 2, real-plasma spike: the pl806 colon-cancer plasma (250 ng and 50 ng arrays) mixed
    into a blank at plasma fraction p = f / c(pl806). The paper's 3 % plasma point
    (p = 0.03, f = 0.03 * c(pl806)) is added to this grid.
Replicates = backgrounds x spike replicates. A sensitivity analysis adds technical
re-measurement noise: (rep1 - rep2) / sqrt(2) beta vectors from the three sepsis 250 ng
duplicate pairs, with both signs, applied to blanks and mixtures alike (7 noise states
including none). This doubles the technical noise variance of the background part, so it
is conservative.

Also: the real GEO wet-lab mixes against the default-split LoB (mix_7 is shown but is not
an LoD anchor by design decision), and an off-target analysis of the
~2 % neuronal signal (Moss et al. 2018, Results "Tissue origins of healthy cfDNA", Fig. 4).

Usage: uv run python scripts/phase2_lob_lod.py [--data DIR] [--out results/phase2]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import deconv, lod

SEED = 20261008
N_BOOT = 10_000
COLON, NEURON = "Colon_epithelial_cells", "Cortical_neurons"

POOLS = {
    "GSM3455776": "Young M1", "GSM3455777": "Young M2", "GSM3455778": "Young F1",
    "GSM3455779": "Young F2", "GSM3455780": "Old M1", "GSM3455781": "Old M2",
    "GSM3455782": "Old F1", "GSM3455783": "Old F2",
}
POOL_PAIRS = [("GSM3455776", "GSM3455777"), ("GSM3455778", "GSM3455779"),
              ("GSM3455780", "GSM3455781"), ("GSM3455782", "GSM3455783")]
INDIV = {"GSM3455807": "HI1", "GSM3455840": "HI2", "GSM3455846": "HI3", "GSM3455848": "HI4"}
BLANKS = {**POOLS, **INDIV}
LEUK = "GSM3455862"  # leukocytes_1, HI2's matched leukocyte DNA (pure-blood blank)
COLON_REPS = ["GSM3455852", "GSM3455863", "GSM3455864"]
NEURON_REPS = [("450K", "GSM3455773"), ("EPIC", "GSM3455834"), ("EPIC", "GSM3455835")]
PL806 = ["GSM3455808", "GSM3455810"]
SEPSIS_DUPS = [("GSM3455786", "GSM3455787"), ("GSM3455790", "GSM3455791"),
               ("GSM3455793", "GSM3455794")]
REAL_MIXES = [  # GSM, GEO title, background GSM, nominal plasma fraction p, LoD anchor?
    ("GSM3455837", "mix_5", "GSM3455840", 0.10, True),
    ("GSM3455838", "mix_6", "GSM3455840", 0.05, True),
    ("GSM3455839", "mix_7", "GSM3455840", 0.03, False),
    ("GSM3455845", "mix_1", "GSM3455846", 0.05, True),
    ("GSM3455847", "mix_3", "GSM3455848", 0.05, True),
]
F_GRID = [0.0, 0.001, 0.003, 0.01, 0.03, 0.10]  # study f grid
F_FINE = [0.015, 0.02]  # added between 1 and 3 % to locate the empirical LoD
LOW_LEVELS = [0.001, 0.003, 0.01, 0.03]  # EP17 low-level samples for the pooled SD

SPLITS = {
    "default": {"lob": list(POOLS) + ["GSM3455807"],
                "eval": ["GSM3455840", "GSM3455846", "GSM3455848"]},
    "reverse": {"lob": ["GSM3455840", "GSM3455846", "GSM3455848"],
                "eval": list(POOLS) + ["GSM3455807"]},
}


def default_data() -> Path:
    env = os.environ.get("CFDNA_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


class Decon:
    """NNLS against one fixed atlas, on numpy vectors aligned to the atlas rows."""

    def __init__(self, atl: pd.DataFrame):
        self.cols = list(atl.columns)
        self.a = atl.to_numpy(dtype=float)
        self.i_colon = self.cols.index(COLON)
        self.i_neuron = self.cols.index(NEURON) if NEURON in self.cols else None

    def fracs(self, y: np.ndarray) -> np.ndarray:
        return lod.nnls_fractions(y, self.a)

    def colon(self, y: np.ndarray) -> float:
        return float(self.fracs(y)[self.i_colon])


def collapse_pairs(colon: pd.Series) -> pd.Series:
    """Average each rep1/rep2 pool pair into one unit; keep other samples as they are."""
    out = {}
    paired = {g for p in POOL_PAIRS for g in p}
    for a, b in POOL_PAIRS:
        if a in colon.index and b in colon.index:
            out[f"{a}+{b}"] = (colon[a] + colon[b]) / 2
        else:
            out.update({g: colon[g] for g in (a, b) if g in colon.index})
    for g in colon.index:
        if g not in paired:
            out[g] = colon[g]
    return pd.Series(out)


def lob_row(name: str, colon: pd.Series, rng: np.random.Generator, n_units_note: str) -> dict:
    v = colon.to_numpy()
    lo, hi = lod.bootstrap_ci(v, lod.lob_nonparametric, N_BOOT, rng)
    return {"blank_set": name, "n": len(v), "units": n_units_note,
            "mean": v.mean(), "sd": v.std(ddof=1), "min": v.min(), "max": v.max(),
            "lob": lod.lob_nonparametric(v), "lob_ci_lo": lo, "lob_ci_hi": hi,
            "lob_parametric": v.mean() + 1.645 * v.std(ddof=1)}


def run_option(option: str, dec_for_spike, spikes: dict, f_to_mix, betas: dict,
               noise: list[np.ndarray], f_grid: list[float]) -> pd.DataFrame:
    """Deconvolve every (split, background, spike replicate, f, noise state) combination.

    dec_for_spike(spike_id) -> Decon used for that spike (leave-one-out atlas for option 1).
    f_to_mix(f) -> mixing weight of the spike profile.
    Blanks (f = 0) are deconvolved once per (spike atlas, background, noise state).
    """
    rows = []
    for sid, spike in spikes.items():
        dec = dec_for_spike(sid)
        for bg in BLANKS:
            for k, eps in enumerate(noise):
                for f in f_grid:
                    y = lod.mix(betas[bg], spike, f_to_mix(f))
                    y = np.clip(y + eps, 0, 1)
                    rows.append({"option": option, "spike": sid, "background": bg,
                                 "noise": k, "f": f, "colon": dec.colon(y)})
    return pd.DataFrame(rows)


def evaluate(mixes: pd.DataFrame, option: str, noise_states: list[int],
             grid: list[float]) -> tuple[list[dict], list[dict]]:
    """LoB, detection rates and LoDs for one option, both splits, one noise setting."""
    lod_rows, det_rows = [], []
    label = "none" if noise_states == [0] else "re-measurement"
    for split, sets in SPLITS.items():
        d = mixes[mixes.noise.isin(noise_states)]
        # LoB per spike replicate (option 1 atlases differ per replicate), then detection.
        lob_by_spike = lod.lob_by_spike(d, sets["lob"])
        ev = d[d.background.isin(sets["eval"])].copy()
        ev["lob"] = ev.spike.map(lob_by_spike)
        ev["detected"] = ev.colon > ev.lob
        rates = ev.groupby("f").detected.mean()
        n_rep = ev.groupby("f").size()
        for f in grid:
            det_rows.append({"option": option, "split": split, "noise": label, "f": f,
                             "n_replicates": int(n_rep[f]), "detected": int(
                                 ev[ev.f == f].detected.sum()),
                             "rate": rates[f], "mean_colon": ev[ev.f == f].colon.mean(),
                             "sd_colon": ev[ev.f == f].colon.std(ddof=1)})
        lob = float(lob_by_spike.mean())
        groups = [ev[ev.f == f].colon.to_numpy() for f in LOW_LEVELS]
        par = lod.lod_parametric(lob, groups)
        cal = stats.linregress(ev[ev.f > 0].f, ev[ev.f > 0].colon)
        emp = lod.empirical_lod(list(rates.index), list(rates.values))
        lod_rows.append({
            "option": option, "split": split, "noise": label,
            "lob_blanks": len(sets["lob"]), "eval_backgrounds": len(sets["eval"]),
            "spike_replicates": ev.spike.nunique(),
            "replicates_per_level": int(n_rep.iloc[0]),
            "lob_colon": lob, "lob_range_over_spikes": float(lob_by_spike.max()
                                                              - lob_by_spike.min()),
            "sd_low_pooled": par["sd_pooled"], "c_p": par["c_p"],
            "lod_param_colon": par["lod"],
            "cal_intercept": cal.intercept, "cal_slope": cal.slope,
            "lod_param_f": (par["lod"] - cal.intercept) / cal.slope,
            "lod_empirical_f": emp,
            "fpr_at_f0": rates[0.0],
        })
    return lod_rows, det_rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=default_data())
    ap.add_argument("--out", type=Path, default=Path("results/phase2"))
    ap.add_argument("--fig", type=Path, default=Path("figures/phase2_lob_lod.png"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    atl = deconv.dedup_atlas(atlas_mod.read_atlas(args.data))
    epic = pd.read_parquet(args.data / "geo" / "betas_epic.parquet").reindex(atl.index)
    k450 = pd.read_parquet(args.data / "geo" / "betas_450k.parquet").reindex(atl.index)
    betas = {g: epic[g].to_numpy(dtype=float) for g in epic.columns}
    for g in k450.columns:
        betas.setdefault(g, k450[g].to_numpy(dtype=float))
    print(f"atlas {atl.shape}; EPIC arrays {epic.shape[1]}")
    checks = []

    # -- 0. Reference checks -------------------------------------------------------------
    colon_reps = epic[COLON_REPS]
    diff = (atl[COLON] - colon_reps.mean(axis=1)).abs()
    checks.append(("atlas colon column = mean of GSM3455852/63/64, max abs diff",
                   f"{diff.max():.2e}"))
    neu = pd.concat([k450["GSM3455773"], epic["GSM3455834"], epic["GSM3455835"]], axis=1)
    diff = (atl[NEURON] - neu.mean(axis=1)).abs()
    checks.append(("atlas neuron column = mean of GSM3455773 (450K)/834/835, max abs diff",
                   f"{diff.max():.2e}"))

    def loo_atlas(col: str, reps: pd.DataFrame, drop: str) -> pd.DataFrame:
        a = atl.copy()
        a[col] = reps.drop(columns=drop).mean(axis=1).fillna(atl[col])
        return a

    dec_std = Decon(atl)
    dec_loo = {r: Decon(loo_atlas(COLON, colon_reps, r)) for r in COLON_REPS}

    # -- 1. Blank colon fractions and LoB -------------------------------------------------
    blank_fr = pd.DataFrame({g: dec_std.fracs(betas[g]) for g in [*BLANKS, LEUK]},
                            index=atl.columns).T
    colon = blank_fr[COLON]
    lob_rows = []
    sets = {
        "all 12 blanks": list(BLANKS),
        "default-split LoB set (8 pools + HI1)": SPLITS["default"]["lob"],
        "reverse-split LoB set (HI2-HI4)": SPLITS["reverse"]["lob"],
    }
    for name, gs in sets.items():
        c = colon[gs]
        lob_rows.append(lob_row(name, c, rng, "arrays"))
        cc = collapse_pairs(c)
        if len(cc) < len(c):
            lob_rows.append(lob_row(name + ", pool pairs collapsed", cc, rng,
                                    "pool-pair means + individuals"))
    lob_tab = pd.DataFrame(lob_rows)
    lob_tab.to_csv(args.out / "lob.tsv", sep="\t", index=False)
    blank_tab = pd.DataFrame({"label": {**BLANKS, LEUK: "HI2 leukocytes"}, "colon": colon,
                              "neuron": blank_fr[NEURON]})
    blank_tab.index.name = "gsm"
    blank_tab.to_csv(args.out / "blank_colon.tsv", sep="\t")
    print(lob_tab.to_string())

    # Are rep1/rep2 pool pairs closer than unrelated pools? Compare with the sepsis duplicates.
    pair_rows = []
    for kind, pairs in [("pool rep1/rep2", POOL_PAIRS), ("sepsis 250 ng duplicates",
                                                         SEPSIS_DUPS)]:
        for a, b in pairs:
            fa, fb = dec_std.colon(betas[a]), dec_std.colon(betas[b])
            d = betas[a] - betas[b]
            pair_rows.append({"kind": kind, "a": a, "b": b, "median_abs_dbeta":
                              float(np.nanmedian(np.abs(d))), "colon_a": fa, "colon_b": fb,
                              "abs_dcolon": abs(fa - fb)})
    unrel = []
    pools = list(POOLS)
    for i in range(len(pools)):
        for j in range(i + 1, len(pools)):
            if (pools[i], pools[j]) in POOL_PAIRS:
                continue
            d = betas[pools[i]] - betas[pools[j]]
            unrel.append({"kind": "unrelated pools", "a": pools[i], "b": pools[j],
                          "median_abs_dbeta": float(np.nanmedian(np.abs(d))),
                          "colon_a": colon[pools[i]], "colon_b": colon[pools[j]],
                          "abs_dcolon": abs(colon[pools[i]] - colon[pools[j]])})
    pair_tab = pd.DataFrame(pair_rows + unrel)
    pair_tab.to_csv(args.out / "replicate_pairs.tsv", sep="\t", index=False)

    # -- 2. Noise vectors from the sepsis duplicates -------------------------------------
    noise = [np.zeros(len(atl))]
    for a, b in SEPSIS_DUPS:
        e = np.nan_to_num((betas[a] - betas[b]) / np.sqrt(2))
        noise += [e, -e]
    noise_all = list(range(len(noise)))
    # Positive control for the noise model: colon-fraction spread it induces on blanks vs the
    # spread between the real duplicates.
    sim_sd = np.std([dec_std.colon(np.clip(betas[g] + e, 0, 1)) - colon[g]
                     for g in BLANKS for e in noise[1:]], ddof=1)
    real = pair_tab[pair_tab.kind == "sepsis 250 ng duplicates"].abs_dcolon.to_numpy()
    checks.append(("colon SD added by one noise draw (12 blanks x 6 draws)",
                   f"{sim_sd * 100:.3f} points"))
    checks.append((("per-array technical colon SD from the 3 sepsis duplicates, "
                   "sqrt(mean(d^2)/2)"), f"{np.sqrt(np.mean(real ** 2) / 2) * 100:.3f} points"))

    # -- 3. In-silico mixtures ----------------------------------------------------------
    c806 = float(np.mean([dec_std.colon(betas[g]) for g in PL806]))
    checks.append(("c(pl806), mean colon fraction of the two pl806 arrays", f"{c806:.4f}"))
    grid1 = sorted(F_GRID + F_FINE)
    grid2 = sorted(grid1 + [0.03 * c806])
    opt1 = run_option("1: colon-cell spike", lambda s: dec_loo[s],
                      {r: betas[r] for r in COLON_REPS}, lambda f: f, betas, noise, grid1)
    opt2 = run_option("2: pl806 plasma spike", lambda s: dec_std,
                      {g: betas[g] for g in PL806}, lambda f: f / c806, betas, noise, grid2)
    mixes = pd.concat([opt1, opt2], ignore_index=True)
    mixes.to_csv(args.out / "mixes_long.tsv.gz", sep="\t", index=False)

    # Known-answer positive control: noise-free recovery of the added colon fraction.
    rec_rows = []
    for opt, d in mixes[mixes.noise == 0].groupby("option"):
        base = d[d.f == 0].set_index(["spike", "background"]).colon
        d = d[d.f > 0].copy()
        d["added"] = d.colon - d.set_index(["spike", "background"]).index.map(base)
        for f, g in d.groupby("f"):
            rec_rows.append({"option": opt, "f": f, "n": len(g),
                             "mean_added": g.added.mean(), "min_added": g.added.min(),
                             "max_added": g.added.max(), "recovery": g.added.mean() / f})
    rec_tab = pd.DataFrame(rec_rows)
    rec_tab.to_csv(args.out / "recovery.tsv", sep="\t", index=False)
    print(rec_tab.to_string())

    # -- 4. LoB / LoD per option, split and noise setting --------------------------------
    lod_rows, det_rows = [], []
    for opt, grid in [(opt1, grid1), (opt2, grid2)]:
        name = opt.option.iloc[0]
        for ns in ([0], noise_all):
            a, b = evaluate(opt, name, ns, grid)
            lod_rows += a
            det_rows += b
    lod_tab = pd.DataFrame(lod_rows)
    det_tab = pd.DataFrame(det_rows)
    lod_tab.to_csv(args.out / "lod.tsv", sep="\t", index=False)
    det_tab.to_csv(args.out / "detection.tsv", sep="\t", index=False)
    print(lod_tab.to_string())

    # -- 5. Real wet-lab mixes against the default-split LoB -----------------------------
    lob_default = lod.lob_nonparametric(colon[SPLITS["default"]["lob"]])
    real_rows = []
    for gsm, title, bg, p, anchor in REAL_MIXES:
        c = dec_std.colon(betas[gsm])
        real_rows.append({"gsm": gsm, "title": title, "background": INDIV[bg],
                          "nominal_p": p, "expected_colon": p * c806 + (1 - p) * colon[bg],
                          "observed_colon": c, "lob_default": lob_default,
                          "detected": c > lob_default, "lod_anchor": anchor,
                          "margin_points": (c - lob_default) * 100})
        hits = [dec_std.colon(np.clip(betas[gsm] + e, 0, 1)) > lob_default for e in noise]
        real_rows[-1]["detected_under_noise"] = f"{sum(hits)}/{len(hits)}"
    real_tab = pd.DataFrame(real_rows)
    real_tab.to_csv(args.out / "real_mixes.tsv", sep="\t", index=False)
    print(real_tab.to_string())

    # Leave-one-blank-out: LoB from the other 11 blanks, then the spike fraction this blank
    # needs to exceed it, f_needed = (LoB_-b - colon_b) / recovery. Pool pairs are left out
    # together so a pool's twin never sets its LoB.
    rec1 = rec_tab[rec_tab.option.str.startswith("1")].recovery.mean()
    rec2 = rec_tab[rec_tab.option.str.startswith("2")].recovery.mean()
    twin = {a: b for a, b in POOL_PAIRS} | {b: a for a, b in POOL_PAIRS}
    lobo_rows = []
    for g, label in BLANKS.items():
        rest = [h for h in BLANKS if h != g and h != twin.get(g)]
        lob_g = lod.lob_nonparametric(colon[rest])
        need = lob_g - colon[g]
        lobo_rows.append({"gsm": g, "label": label, "n_lob_blanks": len(rest),
                          "colon": colon[g], "lob_without": lob_g,
                          "f_needed_opt1": max(need, 0) / rec1,
                          "f_needed_opt2": max(need, 0) / rec2,
                          "false_positive": colon[g] > lob_g})
    lobo = pd.DataFrame(lobo_rows)
    lobo.to_csv(args.out / "leave_one_blank_out.tsv", sep="\t", index=False)
    print(lobo.to_string())
    for c in ("f_needed_opt1", "f_needed_opt2"):
        checks.append((f"leave-one-blank-out {c}: median / 0.95 rank / max",
                       (f"{lobo[c].median() * 100:.2f} / "
                       f"{lod.lob_nonparametric(lobo[c]) * 100:.2f} / "
                       f"{lobo[c].max() * 100:.2f} %")))

    # -- 6. Off-target signal: neurons and colon ----------------------------------------
    off_rows = []
    others = [c for c in atl.columns if c != COLON]
    base_lob = {k: lod.lob_nonparametric(colon[v]) for k, v in sets.items()}
    for drop in others:
        dec = Decon(atl.drop(columns=drop))
        cd = pd.Series({g: dec.colon(betas[g]) for g in BLANKS})
        row = {"dropped": drop, "mean_blank_colon": cd.mean(),
               "delta_mean_points": (cd.mean() - colon[list(BLANKS)].mean()) * 100}
        for k, v in sets.items():
            row[f"delta_lob_points[{k}]"] = (lod.lob_nonparametric(cd[v]) - base_lob[k]) * 100
        off_rows.append(row)
    off_tab = pd.DataFrame(off_rows).sort_values("delta_mean_points", key=abs,
                                                 ascending=False)
    off_tab.to_csv(args.out / "atlas_dropout.tsv", sep="\t", index=False)

    cor_rows = []
    for c in others:
        if blank_fr.loc[list(BLANKS), c].std() == 0:
            continue
        r, p = stats.spearmanr(blank_fr.loc[list(BLANKS), COLON], blank_fr.loc[list(BLANKS), c])
        cor_rows.append({"cell_type": c, "spearman_rho_vs_colon": r, "p": p,
                         "mean_blank_fraction": blank_fr.loc[list(BLANKS), c].mean()})
    cor_tab = pd.DataFrame(cor_rows).sort_values("p")
    cor_tab.to_csv(args.out / "blank_colon_correlations.tsv", sep="\t", index=False)

    # Crosstalk: spike neurons (leave-one-out neuron column) and read colon; and the reverse.
    neu_reps = neu.copy()
    neu_reps.columns = [g for _, g in NEURON_REPS]
    xt_rows = []
    for _, r in NEURON_REPS:
        dec = Decon(loo_atlas(NEURON, neu_reps, r))
        for f in (0.0, 0.02, 0.05):
            for bg in BLANKS:
                fr = dec.fracs(lod.mix(betas[bg], betas[r], f))
                xt_rows.append({"spike": "neuron", "replicate": r, "background": bg, "f": f,
                                "colon": fr[dec.i_colon], "neuron": fr[dec.i_neuron]})
    for r in COLON_REPS:
        dec = dec_loo[r]
        for f in (0.0, 0.02, 0.05):
            for bg in BLANKS:
                fr = dec.fracs(lod.mix(betas[bg], betas[r], f))
                xt_rows.append({"spike": "colon", "replicate": r, "background": bg, "f": f,
                                "colon": fr[dec.i_colon], "neuron": fr[dec.i_neuron]})
    xt = pd.DataFrame(xt_rows)
    base = xt[xt.f == 0].set_index(["spike", "replicate", "background"])
    xt = xt[xt.f > 0].join(base[["colon", "neuron"]], on=["spike", "replicate", "background"],
                           rsuffix="_0")
    xt["d_colon_points"] = (xt.colon - xt.colon_0) * 100
    xt["d_neuron_points"] = (xt.neuron - xt.neuron_0) * 100
    xt_sum = xt.groupby(["spike", "f"]).agg(
        n=("d_colon_points", "size"),
        d_colon_mean=("d_colon_points", "mean"), d_colon_min=("d_colon_points", "min"),
        d_colon_max=("d_colon_points", "max"), d_neuron_mean=("d_neuron_points", "mean"),
        d_neuron_min=("d_neuron_points", "min"), d_neuron_max=("d_neuron_points", "max"),
    ).reset_index()
    xt_sum.to_csv(args.out / "crosstalk.tsv", sep="\t", index=False)
    print(xt_sum.to_string())
    print(off_tab.head(8).to_string())
    print(cor_tab.head(6).to_string())

    pd.DataFrame(checks, columns=["check", "value"]).to_csv(args.out / "checks.tsv", sep="\t",
                                                            index=False)
    for c, v in checks:
        print(f"{c}: {v}")

    plot(det_tab, lod_tab, colon, real_tab, lob_default, args.fig)


def plot(det: pd.DataFrame, lods: pd.DataFrame, colon: pd.Series, real: pd.DataFrame,
         lob_default: float, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    x0 = 3e-4  # where f = 0 is drawn on the log axis
    styles = {("default", "none"): {"color": "#1f78b4", "ls": "-", "marker": "o",
                                         "label": "default split"},
              ("default", "re-measurement"): {"color": "#a6cee3", "ls": "-", "marker": "s",
                                                  "label": "default split + noise"},
              ("reverse", "none"): {"color": "#e31a1c", "ls": "--", "marker": "^",
                                        "label": "reverse split"}}
    for ax, opt in zip(axes[:2], det.option.unique()):
        d = det[det.option == opt]
        for (split, noise), st in styles.items():
            s = d[(d.split == split) & (d.noise == noise)]
            x = np.where(s.f == 0, x0, s.f)
            ax.plot(x * 100, s.rate * 100, ms=5, **st)
            row = lods[(lods.option == opt) & (lods.split == split) & (lods.noise == noise)]
            emp = row.lod_empirical_f.iloc[0]
            if np.isfinite(emp):
                ax.axvline(emp * 100, color=st["color"], ls=":", lw=1)
        row = lods[(lods.option == opt) & (lods.split == "default") & (lods.noise == "none")]
        ax.axvline(row.lod_param_f.iloc[0] * 100, color="k", lw=1)
        ax.axhline(95, color="grey", lw=0.8)
        ax.set_xscale("log")
        ax.set_xticks([x0 * 100, 0.1, 0.3, 1, 2, 3, 10])
        ax.set_xticklabels(["0", "0.1", "0.3", "1", "2", "3", "10"])
        ax.set_xlabel("spiked colon fraction f (%)")
        ax.set_ylabel("detection rate (% of replicates > LoB)")
        ax.set_title(f"Option {opt}", fontsize=10)
        ax.set_ylim(-3, 103)
    axes[0].legend(fontsize=8, loc="upper left")
    axes[1].text(0.98, 0.03, "black line: parametric LoD (default)\n"
                 "dotted: empirical LoD", transform=axes[1].transAxes, ha="right",
                 fontsize=7)

    ax = axes[2]
    lob_set = SPLITS["default"]["lob"]
    ev_set = SPLITS["default"]["eval"]
    ax.scatter(np.zeros(len(lob_set)) - 0.1, colon[lob_set] * 100, color="#1f78b4",
               label="LoB blanks (8 pools + HI1)")
    ax.scatter(np.zeros(len(ev_set)) + 0.1, colon[ev_set] * 100, color="#33a02c",
               label="HI2-HI4 (mix backgrounds)")
    for _, r in real.iterrows():
        x = r.nominal_p * 100
        ax.scatter(x, r.observed_colon * 100, color="#ff7f00" if r.lod_anchor else "white",
                   edgecolor="#ff7f00", zorder=3)
        ax.annotate(r.title, (x, r.observed_colon * 100), fontsize=7, xytext=(4, 2),
                    textcoords="offset points")
    ax.axhline(lob_default * 100, color="k", lw=1, label=f"LoB {lob_default * 100:.2f} %")
    ax.set_xlabel("nominal pl806 plasma fraction p (%)")
    ax.set_ylabel("NNLS colon fraction (%)")
    ax.set_title("Real GEO mixes vs default-split LoB\n(open marker: mix_7, not an anchor)",
                 fontsize=10)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"figure -> {path}")


if __name__ == "__main__":
    main()
