#!/usr/bin/env python3
"""Per-sample QC table for the Phase 5 pipeline.

Reads the per-run outputs of the pipeline from one directory and writes one row per run:
read pairs, duplication, mapping, insert size, bisulfite conversion (CHH calls on one
chromosome), M-bias range, genome-wide CpG calls, and CpG calls in the marker blocks
(duplicates removed). Needs numpy and pandas only.
"""

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def fastp_pairs(path):
    s = json.loads(Path(path).read_text())["summary"]
    return s["before_filtering"]["total_reads"] // 2, s["after_filtering"]["total_reads"] // 2


def markdup_rate(path):
    kv = dict(re.findall(r"^([A-Z_ ]+): (\d+)", Path(path).read_text(), flags=re.MULTILINE))
    reads = int(kv["READ"]) - int(kv.get("EXCLUDED", 0))
    return int(kv["DUPLICATE TOTAL"]) / reads if reads else np.nan


def samtools_sn(path):
    sn = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("SN\t"):
            _, key, val = line.split("\t")[:3]
            sn[key.rstrip(":")] = float(val)
    return sn


def mbias_range(path):
    """Max minus min CpG methylation (%) over read positions, per read (strands summed)."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df = df[df["Context"] == "CpG"] if "Context" in df.columns else df
    out = {}
    for read, g in df.groupby("Read"):
        pos = g.groupby("Position")[["nMethylated", "nUnmethylated"]].sum()
        tot = pos.sum(axis=1)
        pct = 100 * pos["nMethylated"] / tot.where(tot >= 100)
        out[f"mbias_r{read}_range_pct"] = float(pct.max() - pct.min())
    return out


def bedgraph_calls(path):
    calls = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.split("\t")
            calls += int(f[4]) + int(f[5])
    return calls


def chh_conversion(path):
    kv = dict(line.split("\t") for line in Path(path).read_text().split("\n") if line)
    m, u = int(kv["chh_meth"]), int(kv["chh_unmeth"])
    return u / (m + u) if m + u else np.nan


def block_calls(pat_path, blocks):
    """CpG calls and fragments per block of one marker set from a wgbstools pat file."""
    starts = blocks["startCpG"].to_numpy()
    ends = blocks["endCpG"].to_numpy()
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]
    calls = np.zeros(len(blocks), dtype=np.int64)
    frags = np.zeros(len(blocks), dtype=np.int64)
    with gzip.open(pat_path, "rt") as fh:
        for line in fh:
            _, s, pat, n = line.rstrip("\n").split("\t")[:4]
            s, n = int(s), int(n)
            # Blocks of one set are disjoint (CpG indices are genome-wide); a fragment can span two.
            i = np.searchsorted(starts, s + len(pat), side="left") - 1
            while i >= 0 and ends[i] > s:
                lo, hi = max(s, starts[i]), min(s + len(pat), ends[i])
                k = sum(c in "CT" for c in pat[lo - s : hi - s])
                if k:
                    calls[i] += k * n
                    frags[i] += n
                i -= 1
    res = pd.DataFrame({"calls": calls, "frags": frags})
    res.index = blocks.index[order]
    return res.sort_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("indir")
    ap.add_argument("--blocks", required=True, help="blocks_cpg.tsv from wgbstools convert")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = Path(a.indir)
    blocks = pd.read_csv(
        a.blocks, sep="\t", header=None,
        names=["chrom", "start", "end", "startCpG", "endCpG", "block", "set", "n_probes"],
    )
    rows = []
    for pat in sorted(d.glob("*.dedup.pat.gz")):
        run = pat.name.removesuffix(".dedup.pat.gz")
        r = {"run": run}
        if (d / f"{run}.fastp.json").exists():
            r["pairs_raw"], r["pairs_trimmed"] = fastp_pairs(d / f"{run}.fastp.json")
        if (d / f"{run}.markdup.txt").exists():
            r["dup_rate"] = markdup_rate(d / f"{run}.markdup.txt")
        if (d / f"{run}.stats.txt").exists():
            sn = samtools_sn(d / f"{run}.stats.txt")
            r["mapped_frac"] = sn["reads mapped"] / sn["raw total sequences"]
            r["insert_mean"] = sn["insert size average"]
            r["insert_sd"] = sn["insert size standard deviation"]
        r["chh_conversion"] = chh_conversion(d / f"{run}.chh.txt")
        r.update(mbias_range(d / f"{run}.mbias.txt"))
        r["cpg_calls_genome"] = bedgraph_calls(d / f"{run}_CpG.bedGraph.gz")
        for s, b in blocks.groupby("set"):
            b = block_calls(pat, b)
            r[f"{s}_calls_per_block"] = b["calls"].mean()
            r[f"{s}_frags_per_block"] = b["frags"].mean()
            r[f"{s}_blocks_covered"] = (b["calls"] > 0).mean()
        rows.append(r)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
