"""Phase 0: atlas CpGs present (non-NA) in the TCGA-COAD SeSAMe level-3 betas.

Positive control: cg00000029 (not an atlas CpG) is read by the same code path.
Usage: uv run python scripts/phase0_tcga_atlas_check.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import (
                           deconv,
                           tcga,
)

data = Path(os.environ.get("CFDNA_DATA",
                           Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"))
atl = deconv.dedup_atlas(atlas_mod.read_atlas(data))
man = tcga.read_manifest(data)
probes = list(atl.index) + ["cg00000029"]
b = tcga.read_betas(data, probes=probes)
print(f"files {b.shape[1]} (manifest {len(man)}), probe rows read {b.shape[0]}")
ab = b.reindex(atl.index)
print(f"atlas CpGs as rows in TCGA files: {ab.notna().any(axis=1).sum()} of {len(atl)}")
print(f"atlas CpGs non-NA in all files: {ab.notna().all(axis=1).sum()}")
print(f"atlas CpGs non-NA in >= 95% of files: {(ab.notna().mean(axis=1) >= 0.95).sum()}")
per = ab.notna().sum()
print(f"non-NA atlas CpGs per file: min {per.min()}, median {int(per.median())}, max {per.max()}")
pos = b.loc["cg00000029"] if "cg00000029" in b.index else None
print("positive control cg00000029 non-NA files:", None if pos is None else int(pos.notna().sum()))
