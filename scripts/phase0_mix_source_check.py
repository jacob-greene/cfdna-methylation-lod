"""Phase 0 check: which cancer plasma is in each GEO cfDNA mix?

For mix M over healthy background H, M - H = p * (C - H) if C is the spiked plasma.
We report corr(M - H, C - H) and the least-squares p for each candidate C.
GEO says the colon plasma is pl806; Moss et al. 2018 Supp. Data 1 Table 8 says CC2.
"""

from pathlib import Path

import numpy as np
import pandas as pd

D = Path(__file__).resolve().parents[1].parent / "cfdna-methylation-lod-data"
CANDS = {"pl806_250": "GSM3455808", "pl806_50": "GSM3455810", "pl684": "GSM3455809",
         "pl2010": "GSM3455842", "lung_1478": "GSM3455841", "breast_1040": "GSM3455843"}
MIXES = {"mix_5 10% HI2": ("GSM3455837", "GSM3455840"), "mix_6 5% HI2": ("GSM3455838", "GSM3455840"),
         "mix_7 3% HI2": ("GSM3455839", "GSM3455840"), "mix_1 5% HI3": ("GSM3455845", "GSM3455846"),
         "mix_3 5% HI4": ("GSM3455847", "GSM3455848")}

b = pd.read_parquet(D / "geo/betas_epic.parquet")
rows = []
for mn, (m, h) in MIXES.items():
    for cn, c in CANDS.items():
        d = b[[m, h, c]].dropna()
        y, x = d[m] - d[h], d[c] - d[h]
        rows.append((mn, cn, float((x * y).sum() / (x * x).sum()), float(np.corrcoef(x, y)[0, 1])))
t = pd.DataFrame(rows, columns=["mix", "candidate", "fitted_p", "corr"])
out = Path(__file__).resolve().parents[1] / "results/phase0/mix_source_check.tsv"
t.to_csv(out, sep="\t", index=False, float_format="%.4f")
print(t.pivot(index="mix", columns="candidate", values="corr").round(3).to_string())
