"""Phase 5 fix pass: compare original and corrected spike dosing; recompute the gap attribution.

The original scores (results/phase5/original_dosing/) dosed the spike with a pooled fragment
count that merged donors' pat row ids, so k >= 2 pools were under-dosed (fixed in
src/phase5.n_fragments_pool). The corrected scores are results/phase5/scores.tsv.gz.

Outputs:
  - results/phase5/dosing_comparison.tsv: per rung x score x set x k, detection rate at
    f = 0.1, AUC at f = 0.1 and LoD, original vs corrected dosing;
  - results/phase5/gap_attribution.tsv: a decomposition of the real vs binom_real
    gap (LLR_binom, eval tests, f_alt = f), for both dosings. shift = mean paired
    (spiked - f = 0) score of the same pool and f_alt; z = shift / SD of the f = 0 eval
    scores; deficit share = log(1 / shift ratio) / log(z_binom_real / z_real).

Usage: uv run python scripts/phase5_dosing_compare.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("results/phase5")
VERSIONS = {"original": OUT / "original_dosing", "corrected": OUT}


def comparison() -> pd.DataFrame:
    parts = []
    for ver, p in VERSIONS.items():
        det = pd.read_csv(p / "detection.tsv", sep="\t")
        lod = pd.read_csv(p / "lod.tsv", sep="\t")
        d = det[np.isclose(det.f, 0.1)][["rung", "score", "set", "k", "rate", "auc"]]
        d = d.merge(lod[["rung", "score", "set", "k", "lod"]], on=["rung", "score", "set", "k"])
        parts.append(d.rename(columns={c: f"{c}_{ver}" for c in ("rate", "auc", "lod")}))
    c = parts[0].merge(parts[1], on=["rung", "score", "set", "k"], validate="1:1")
    return c[["rung", "score", "set", "k", "rate_original", "rate_corrected", "auc_original",
              "auc_corrected", "lod_original", "lod_corrected"]]


def shift_z(d: pd.DataFrame, rung: str, s: str, k: int, f: float) -> tuple[float, float]:
    e = d[(d.kind == "eval") & (d.rung == rung) & (d.score == "LLR_binom") & (d.set == s)
          & (d.k == k) & np.isclose(d.f_alt, f)]
    null = e[e.f == 0].set_index("pool").value
    sp = e[np.isclose(e.f, f)]
    shift = float((sp.value - sp.pool.map(null)).mean())
    return shift, shift / float(null.std(ddof=1))


def gap() -> pd.DataFrame:
    rows = []
    for ver, p in VERSIONS.items():
        d = pd.read_csv(p / "scores.tsv.gz", sep="\t")
        for s in ("atlas", "tumour"):
            for k in (1, 2, 4, 8):
                for f in (0.01, 0.1):
                    sb, zb = shift_z(d, "binom_real", s, k, f)
                    sr, zr = shift_z(d, "real", s, k, f)
                    ratio = sr / sb
                    rows.append((ver, s, k, f, ratio, zb, zr, zr / ratio,
                                 np.log(1 / ratio) / np.log(zb / zr)))
    return pd.DataFrame(rows, columns=["dosing", "set", "k", "f", "shift_ratio", "z_binom_real",
                                       "z_real", "z_real_no_deficit", "deficit_share"])


def main() -> None:
    c = comparison()
    c.to_csv(OUT / "dosing_comparison.tsv", sep="\t", index=False, float_format="%.4g")
    g = gap()
    g.to_csv(OUT / "gap_attribution.tsv", sep="\t", index=False, float_format="%.4g")
    print(c[c.rung == "real"].to_string(index=False))
    print(g.to_string(index=False))


if __name__ == "__main__":
    main()
