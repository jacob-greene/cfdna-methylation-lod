"""Parse GEO series-matrix files for GSE122126 into a sample table and a beta matrix."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

PLATFORMS = {"GPL13534": "450K", "GPL21145": "EPIC"}


def _matrix_path(data_dir: Path, gpl: str) -> Path:
    return Path(data_dir) / "geo" / f"GSE122126-{gpl}_series_matrix.txt.gz"


def read_samples(data_dir: Path) -> pd.DataFrame:
    """Return one row per GSM with title, description, characteristics and platform."""
    frames = []
    for gpl, label in PLATFORMS.items():
        rows: dict[str, list[list[str]]] = {}
        with gzip.open(_matrix_path(data_dir, gpl), "rt") as fh:
            for line in fh:
                if line.startswith("!series_matrix_table_begin"):
                    break
                if not line.startswith("!Sample_"):
                    continue
                key, *vals = line.rstrip("\n").split("\t")
                vals = [v.strip('"') for v in vals]
                rows.setdefault(key[len("!Sample_"):], []).append(vals)
        df = pd.DataFrame({"gsm": rows["geo_accession"][0]})
        for key, lists in rows.items():
            for i, vals in enumerate(lists):
                col = key if len(lists) == 1 else f"{key}_{i + 1}"
                df[col] = vals
        df["platform"] = label
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    keep = ["gsm", "title", "platform", "source_name_ch1"]
    keep += sorted(c for c in out.columns if c.startswith(("characteristics", "description")))
    keep += ["data_processing", "data_row_count"]
    return out[keep].set_index("gsm")


def read_betas(data_dir: Path, gpl: str) -> pd.DataFrame:
    """Return the probe x GSM beta matrix of one platform (NaN where GEO left blanks)."""
    path = _matrix_path(data_dir, gpl)
    with gzip.open(path, "rt") as fh:
        for n, line in enumerate(fh):
            if line.startswith("!series_matrix_table_begin"):
                skip = n + 1
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, index_col=0, comment="!",
                     na_values=["", "null", "NA"], low_memory=False)
    df.index.name = "probe"
    return df.astype("float32")
