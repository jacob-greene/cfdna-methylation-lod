"""Check the Phase 3b rerun against the first version's raw scores.

NNLS and LLR are not refit, so their scores and thresholds must be equal to the first
version (the rerun only adds blank-into-blank rows, drawn from their own generator). The
"<model>_refit" rows reproduce the first version's trained-model rule; LR_L2 and XGB must
match it, LR_L1 may not (it was unseeded before).

Usage: uv run python scripts/phase3b_check_v1.py --old results/phase3b/raw_v1 \
           [--new results/phase3b/raw]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ["fold", "kind", "person", "source", "f", "f_alt"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", type=Path, required=True)
    ap.add_argument("--new", type=Path, default=Path("results/phase3b/raw"))
    args = ap.parse_args()
    rows, bad = [], 0
    for p in sorted(args.new.glob("scores_*.tsv.gz")):
        # A test run on a fold subset (scores_<tag>_f0-1) is checked on those folds only.
        old_p = args.old / re.sub(r"_f[\d-]+(?=\.tsv)", "", p.name)
        if not old_p.exists():
            continue
        new = pd.read_csv(p, sep="\t", dtype={"track": str})
        old = pd.read_csv(old_p, sep="\t", dtype={"track": str})
        new = new[new.kind != "b2b"]
        old = old[old.fold.isin(new.fold.unique())]
        for m_old, g_old in old.groupby("model"):
            m_new = m_old if m_old in ("NNLS", "LLR") else f"{m_old}_refit"
            g_new = new[new.model == m_new]
            if len(g_new) != len(g_old):
                rows.append((p.name, m_old, m_new, len(g_old), len(g_new), np.nan, np.nan))
                bad += m_old in ("NNLS", "LLR")
                continue
            a = g_old[KEY].reset_index(drop=True)
            b = g_new[KEY].reset_index(drop=True)
            same_keys = a.astype(str).equals(b.astype(str))
            ds = np.abs(g_old.score.to_numpy() - g_new.score.to_numpy()).max()
            dt = np.abs(g_old.threshold.to_numpy() - g_new.threshold.to_numpy()).max()
            rows.append((p.name, m_old, m_new, len(g_old), len(g_new), ds, dt))
            if m_old in ("NNLS", "LLR") and not (same_keys and ds == 0 and dt == 0):
                bad += 1
    r = pd.DataFrame(rows, columns=["file", "model_v1", "model_v2", "n_v1", "n_v2",
                                    "max_abs_score_diff", "max_abs_threshold_diff"])
    pd.set_option("display.width", 200)
    print(r.to_string(index=False))
    print(f"NNLS/LLR mismatches: {bad}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
