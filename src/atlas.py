"""Load the Moss et al. 2018 reference atlas and the hg38 probe manifests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_atlas(data_dir: Path, full: bool = False) -> pd.DataFrame:
    """CpG x cell-type mean-beta atlas (7,890 x 25, or ~388.6K x 25 when full=True)."""
    d = Path(data_dir) / "meth_atlas"
    path = d / ("full_atlas.csv.gz" if full else "reference_atlas.csv")
    df = pd.read_csv(path, index_col=0)
    df.index.name = "probe"
    return df


def read_manifest(data_dir: Path, platform: str) -> pd.DataFrame:
    """hg38 probe coordinates. CpG_beg is 0-based; the CpG spans [CpG_beg, CpG_end)."""
    name = {"EPIC": "EPIC", "450K": "HM450"}[platform]
    path = Path(data_dir) / "manifests" / f"{name}.hg38.manifest.tsv.gz"
    cols = ["Probe_ID", "CpG_chrm", "CpG_beg", "CpG_end", "type", "mapQ_A"]
    df = pd.read_csv(path, sep="\t", usecols=cols, low_memory=False)
    return df.set_index("Probe_ID")
