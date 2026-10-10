"""Summarise the Bismark QC outputs of SPOT-MAS runs (Phase 0 item 4).

For each run directory written by `phase0_spotmas_qc.sbatch` this reports
- read pairs, mapping efficiency, duplicate rate (Bismark reports);
- non-CpG conversion = 100 - CHH methylation % (also CHG);
- M-bias: R1/R2 read positions whose CHH methylation exceeds 2 % (Adaptase tail);
- mean depth and the fraction of atlas CpGs (6,105 unique, hg38 from the HM450 manifest)
  with >= 1 call, plus the same for all 450K cg probes as a reference.

Usage: uv run python scripts/phase0_spotmas_qc_summary.py <run_dir> [...] > out.tsv
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import deconv


def grab(text: str, pattern: str) -> float:
    m = re.search(pattern, text)
    return float(m.group(1).replace(",", "")) if m else float("nan")


def mbias(path: Path) -> pd.DataFrame:
    rows, ctx = [], None
    for line in path.read_text().splitlines():
        m = re.match(r"(CpG|CHG|CHH) context \((R1|R2)\)", line)
        if m:
            ctx = m.groups()
            continue
        f = line.split("\t")
        if ctx and len(f) == 5 and f[0].isdigit():
            rows.append({"context": ctx[0], "read": ctx[1], "pos": int(f[0]),
                         "pct": float(f[3]) if f[3] else float("nan"), "cov": int(f[4])})
    return pd.DataFrame(rows)


def high_positions(mb: pd.DataFrame, read: str, thresh: float = 2.0) -> str:
    sub = mb[(mb.context == "CHH") & (mb.read == read) & (mb.pct > thresh) & (mb["cov"] > 100)]
    return ",".join(map(str, sub.pos)) or "none"


def main() -> None:
    data = Path(os.environ.get("CFDNA_DATA",
                               Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"))
    atl = deconv.dedup_atlas(atlas_mod.read_atlas(data))
    man = atlas_mod.read_manifest(data, "450K").dropna(subset=["CpG_chrm", "CpG_beg"])
    man = man[man.index.str.startswith("cg")].copy()
    man["CpG_beg"] = man["CpG_beg"].astype(int)
    is_atlas = man.index.isin(atl.index)
    out = []
    for d in map(Path, sys.argv[1:]):
        rep = next(d.glob("*_PE_report.txt")).read_text()
        dd = next(d.glob("*deduplication_report.txt")).read_text()
        mb = mbias(next(d.glob("*M-bias.txt")))
        stats_txt = next(d.glob("*.samtools_stats.txt")).read_text()
        cov = pd.read_csv(next(d.glob("*.bismark.cov.gz")), sep="\t", header=None,
                          usecols=[0, 1, 4, 5], names=["chr", "pos", "m", "u"])
        calls = (cov.m + cov.u).to_numpy()
        keys = pd.Series(calls, index=pd.MultiIndex.from_arrays([cov.chr, cov.pos]))
        n = (keys.reindex(pd.MultiIndex.from_arrays([man.CpG_chrm, man.CpG_beg + 1])).fillna(0)
             .to_numpy()
             + keys.reindex(pd.MultiIndex.from_arrays([man.CpG_chrm, man.CpG_beg + 2])).fillna(0)
             .to_numpy())
        chh = grab(rep, r"C methylated in CHH context:\s+([\d.]+)%")
        bases = grab(stats_txt, r"SN\tbases mapped \(cigar\):\t(\d+)")
        out.append({
            "run": d.name,
            "read_pairs": int(grab(rep, r"Sequence pairs analysed in total:\s+(\d+)")),
            "mapping_pct": grab(rep, r"Mapping efficiency:\s+([\d.]+)%"),
            "duplicate_pct": grab(dd, r"duplicated alignments removed:\s+\d+ \(([\d.]+)%\)"),
            "CpG_meth_pct": grab(rep, r"C methylated in CpG context:\s+([\d.]+)%"),
            "CHG_meth_pct": grab(rep, r"C methylated in CHG context:\s+([\d.]+)%"),
            "CHH_meth_pct": chh,
            "conversion_pct": 100 - chh,
            "R1_pos_CHH_gt2pct": high_positions(mb, "R1"),
            "R2_pos_CHH_gt2pct": high_positions(mb, "R2"),
            "mean_depth_x": bases / 3.1e9,
            "atlas_cpgs_ge1_call": int((n[is_atlas] >= 1).sum()),
            "atlas_frac_ge1": float((n[is_atlas] >= 1).mean()),
            "all450k_frac_ge1": float((n >= 1).mean()),
            "atlas_mean_calls": float(n[is_atlas].mean()),
        })
    pd.DataFrame(out).to_csv(sys.stdout, sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
