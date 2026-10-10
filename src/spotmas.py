"""SPOT-MAS (PRJNA929650) run metadata joined to the paper's Table S1 labels."""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import pandas as pd


def read_table_s1(docx: Path) -> pd.DataFrame:
    """Parse Table S1 of Nguyen et al. 2023 (one row per subject) from the .docx."""
    xml = zipfile.ZipFile(docx).read("word/document.xml").decode()
    table = re.search(r"<w:tbl>.*?</w:tbl>", xml, flags=re.DOTALL).group(0)  # Table S1 is first
    rows = []
    for tr in re.findall(r"<w:tr[ >].*?</w:tr>", table, flags=re.DOTALL):
        cells = [re.sub(r"<[^>]+>", "", tc).strip()
                 for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.DOTALL)]
        rows.append(cells)
    hi = next(i for i, r in enumerate(rows) if "SampleID" in r)
    # Two header rows; "Background lesions" spans three sub-columns.
    header = ["No", "Group", "Label", "SampleID", "Age", "Gender", "Tumor Stage",
              "Tumor type", "Metastasis", "Number of tumors", "Tumor diameter", "Subtype"]
    body = [r for r in rows[hi + 1:] if r and r[0].isdigit()]
    if any(len(r) != len(header) for r in body):
        raise ValueError("Table S1 row width differs from header width")
    return pd.DataFrame(body, columns=header)


def read_runs(data_dir: Path) -> pd.DataFrame:
    """SRA runinfo joined to BioSample sample names and Table S1 labels."""
    d = Path(data_dir) / "sra"
    runs = pd.read_csv(d / "PRJNA929650_runs_labeled.tsv", sep="\t")
    s1 = read_table_s1(d / "spotmas_supp" / "Table_1.docx")
    return runs.merge(s1, left_on="sample_name", right_on="SampleID", how="left",
                      validate="one_to_one")
