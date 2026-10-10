"""TCGA level-3 SeSAMe beta files (GDC) and Aran et al. 2015 consensus purity."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_manifest(data_dir: Path, project: str = "TCGA-COAD") -> pd.DataFrame:
    return pd.read_csv(Path(data_dir) / "tcga" / project / "manifest.tsv", sep="\t")


def read_betas(data_dir: Path, project: str = "TCGA-COAD", probes=None) -> pd.DataFrame:
    """probe x file-id beta matrix. Restrict to `probes` to save memory."""
    d = Path(data_dir) / "tcga" / project
    man = read_manifest(data_dir, project)
    cols = {}
    for fid, fn in zip(man["id"], man["filename"]):
        s = pd.read_csv(d / "files" / fn, sep="\t", header=None, index_col=0,
                        names=["probe", fid], na_values=["NA"]).iloc[:, 0]
        cols[fid] = s if probes is None else s.reindex(probes)
    return pd.DataFrame(cols).astype("float32")


def read_purity(data_dir: Path, cancer: str = "COAD") -> pd.DataFrame:
    """Aran et al. 2015 Supplementary Data 1: ESTIMATE, ABSOLUTE, LUMP, IHC and CPE per sample."""
    path = Path(data_dir) / "purity" / "aran2015_supplementary_data_1.xlsx"
    p = pd.read_excel(path, skiprows=3).iloc[:, :7]
    p.columns = ["sample", "cancer", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"]
    for c in p.columns[2:]:
        p[c] = pd.to_numeric(p[c], errors="coerce")
    return p[p["cancer"] == cancer].reset_index(drop=True)
