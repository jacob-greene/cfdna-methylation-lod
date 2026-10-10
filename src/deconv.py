"""NNLS deconvolution against the Moss et al. 2018 reference atlas.

Mirrors `deconvolve.py` in github.com/nloyfer/meth_atlas: the atlas is deduplicated by CpG ID,
each sample is joined to the atlas, CpGs with a missing value are dropped per sample, and the
NNLS coefficients are rescaled to sum to 1.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize


def dedup_atlas(atlas: pd.DataFrame) -> pd.DataFrame:
    """Drop repeated CpG rows (keep the first), as deconvolve.py does."""
    return atlas[~atlas.index.duplicated(keep="first")].sort_index()


def nnls_one(sample: pd.Series, atlas: pd.DataFrame) -> tuple[pd.Series, int, float]:
    """Deconvolve one beta vector. Returns (fractions summing to 1, n CpGs used, NNLS residual)."""
    joined = atlas.join(sample.rename("_s"), how="inner").dropna()
    a = joined[atlas.columns].to_numpy(dtype=float)
    y = joined["_s"].to_numpy(dtype=float)
    coef, resid = optimize.nnls(a, y)
    total = coef.sum()
    frac = coef / total if total > 0 else coef
    return pd.Series(frac, index=atlas.columns), len(joined), float(resid)


def deconvolve(betas: pd.DataFrame, atlas: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Deconvolve every column of a probe x sample beta matrix.

    Returns (fractions: sample x cell type, info: sample x [n_cpgs, residual]).
    """
    fracs, info = {}, {}
    for s in betas.columns:
        f, n, r = nnls_one(betas[s], atlas)
        fracs[s] = f
        info[s] = {"n_cpgs": n, "residual": r}
    return pd.DataFrame(fracs).T, pd.DataFrame(info).T


# Groups used in Moss Supp. Data 1, Table 4 (lymphocytes = B + CD4 T + CD8 T + NK).
TABLE4_GROUPS = {
    "Monocytes EPIC": ["Monocytes_EPIC"],
    "Lymphocytes": ["B-cells_EPIC", "CD4T-cells_EPIC", "CD8T-cells_EPIC", "NK-cells_EPIC"],
    "Granulocytes": ["Neutrophils_EPIC"],
}


def to_table4_groups(fracs: pd.DataFrame) -> pd.DataFrame:
    """Collapse atlas cell types to the column groups of Moss Supp. Table 4."""
    out = pd.DataFrame(index=fracs.index)
    used = set()
    for name, cols in TABLE4_GROUPS.items():
        out[name] = fracs[cols].sum(axis=1)
        used.update(cols)
    for c in fracs.columns:
        if c not in used:
            out[c] = fracs[c]
    return out


def synthetic_mix(atlas: pd.DataFrame, props: np.ndarray, noise_sd: float,
                  rng: np.random.Generator) -> pd.Series:
    """Atlas-based synthetic sample: atlas @ props plus Gaussian noise, clipped to [0, 1]."""
    y = atlas.to_numpy(dtype=float) @ props + rng.normal(0, noise_sd, len(atlas))
    return pd.Series(np.clip(y, 0, 1), index=atlas.index)
