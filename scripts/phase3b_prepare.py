"""Phase 3b, step 1: build the TCGA-COAD and GEO inputs of the classifier comparison.

Writes to $CFDNA_DATA/phase3b/ (data, not committed) and results/phase3b/prepare.log:
  geo.parquet       candidate CpG x GEO array betas (blanks, colon replicates, colon plasmas)
  tumours.parquet   candidate CpG x TCGA-COAD primary tumour betas (one file per case,
                    CPE purity >= 0.6), after normalization (see below)
  normals.parquet   candidate CpG x TCGA-COAD solid-tissue normal betas (raw)
  blood.parquet     candidate CpG x Moss full-atlas blood cell types
  tumours.tsv       case, sample, file id, CPE purity
  atlas_cpgs.txt    atlas CpGs (feature set a) inside the candidate set

Candidate CpGs: non-missing in every needed GEO EPIC array, in the Moss full-atlas blood
columns, and in >= 90 % of the kept TCGA tumour files. Remaining TCGA NAs stay NA; the spike
then adds no tumour signal at that CpG (src/classify.spike).

Normalization check (Phase 3b): TCGA betas are SeSAMe; GEO betas are minfi
preprocessIllumina. Compare the TCGA solid-tissue-normal mean with the mean of the three Moss
colon-epithelium replicates at the atlas CpGs. If the median absolute offset exceeds 0.05,
subtract the per-probe offset (estimated on normals only) from every tumour, clipped to [0, 1].

Usage: uv run python scripts/phase3b_prepare.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import atlas as atlas_mod
from src import deconv, tcga

BLANKS = ["GSM3455776", "GSM3455777", "GSM3455778", "GSM3455779", "GSM3455780",
          "GSM3455781", "GSM3455782", "GSM3455783", "GSM3455807", "GSM3455840",
          "GSM3455846", "GSM3455848"]
COLON_REPS = ["GSM3455852", "GSM3455863", "GSM3455864"]
PLASMAS = ["GSM3455808", "GSM3455810", "GSM3455809", "GSM3455842"]  # 806, 806-50ng, 684, 2010
BLOOD = ["Monocytes_EPIC", "B-cells_EPIC", "CD4T-cells_EPIC", "NK-cells_EPIC",
         "CD8T-cells_EPIC", "Neutrophils_EPIC", "Erythrocyte_progenitors"]
PURITY_MIN = 0.6
TCGA_MIN_PRESENT = 0.9
OFFSET_LIMIT = 0.05


def data_dir() -> Path:
    env = os.environ.get("CFDNA_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "cfdna-methylation-lod-data"


def main() -> None:
    data = data_dir()
    out = data / "phase3b"
    out.mkdir(exist_ok=True)
    log_path = Path("results/phase3b/prepare.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    def log(s: str) -> None:
        print(s, flush=True)
        lines.append(s)

    # TCGA samples: primary tumours with CPE >= 0.6, one file per case; all normals.
    man = tcga.read_manifest(data)
    pur = tcga.read_purity(data)
    pur = pur.dropna(subset=["CPE"]).groupby("sample", as_index=False)["CPE"].mean()
    tum = man[man.sample_type == "Primary Tumor"].merge(pur, on="sample", how="left")
    log(f"primary tumour files {len(tum)}, with CPE {tum.CPE.notna().sum()}")
    tum = tum[tum.CPE >= PURITY_MIN].sort_values(["case", "sample", "id"])
    tum = tum.drop_duplicates("case", keep="first")
    log(f"kept tumours (CPE >= {PURITY_MIN}, one per case) {len(tum)}; "
        f"CPE median {tum.CPE.median():.3f}, range {tum.CPE.min():.3f}-{tum.CPE.max():.3f}")
    nor = man[man.sample_type == "Solid Tissue Normal"].sort_values(["case", "sample", "id"])
    nor = nor.drop_duplicates("case", keep="first")
    log(f"solid-tissue normals (one per case) {len(nor)}")

    # GEO and atlas candidates.
    epic = pd.read_parquet(data / "geo" / "betas_epic.parquet")
    need = BLANKS + COLON_REPS + PLASMAS
    geo_ok = epic[need].notna().all(axis=1)
    full = atlas_mod.read_atlas(data, full=True)
    full = full[~full.index.duplicated(keep="first")]
    blood_ok = full[BLOOD].notna().all(axis=1)
    cand = geo_ok[geo_ok].index.intersection(blood_ok[blood_ok].index)
    log(f"CpGs non-missing in all needed GEO arrays {int(geo_ok.sum())}; "
        f"and in the full-atlas blood columns {len(cand)}")

    keep = pd.concat([tum, nor])
    sub = man[man.id.isin(keep.id)]
    d = data / "tcga" / "TCGA-COAD" / "files"
    cols = {}
    for i, (fid, fn) in enumerate(zip(sub["id"], sub["filename"])):
        s = pd.read_csv(d / fn, sep="\t", header=None, index_col=0, names=["probe", fid],
                        na_values=["NA"]).iloc[:, 0]
        cols[fid] = s.reindex(cand).astype("float32")
        if i % 50 == 0:
            print(f"  read {i + 1}/{len(sub)}", flush=True)
    b = pd.DataFrame(cols)
    T = b[tum.id.tolist()]
    N = b[nor.id.tolist()]
    present = T.notna().mean(axis=1)
    cand = present[present >= TCGA_MIN_PRESENT].index
    log(f"candidate CpGs (also non-missing in >= {TCGA_MIN_PRESENT:.0%} of kept tumours) "
        f"{len(cand)}; tumour NA share {float(T.loc[cand].isna().mean().mean()):.4f}")
    T, N = T.loc[cand], N.loc[cand]

    # Normalization check at the atlas CpGs (feature set a).
    atl = deconv.dedup_atlas(atlas_mod.read_atlas(data))
    atl = atl[atl.notna().all(axis=1)]
    a_cpgs = atl.index.intersection(cand)
    moss_colon = epic.loc[cand, COLON_REPS].mean(axis=1)
    offset = (N.mean(axis=1) - moss_colon).astype(float)
    med_atlas = float(offset.loc[a_cpgs].abs().median())
    med_all = float(offset.abs().median())
    log(f"atlas CpGs in candidate set {len(a_cpgs)} of {len(atl)}")
    log(f"median |TCGA normal mean - Moss colon mean|: atlas CpGs {med_atlas:.4f}; "
        f"all candidates {med_all:.4f}; signed median (atlas) "
        f"{float(offset.loc[a_cpgs].median()):+.4f}")
    # Context: the same offset between the three Moss replicates (leave-one-out).
    loo = [float((epic.loc[a_cpgs, r] - epic.loc[a_cpgs, [x for x in COLON_REPS if x != r]]
                  .mean(axis=1)).abs().median()) for r in COLON_REPS]
    log("median |Moss replicate - mean of other two| (atlas CpGs): "
        + ", ".join(f"{v:.4f}" for v in loo))
    if med_atlas > OFFSET_LIMIT:
        T = (T.sub(offset, axis=0)).clip(0, 1)
        log(f"offset > {OFFSET_LIMIT}: per-probe mean-matching applied to tumours "
            "(offset estimated on normals only)")
        after = float((N.sub(offset, axis=0).mean(axis=1) - moss_colon).loc[a_cpgs].abs()
                      .median())
        log(f"after matching, median |normal mean - Moss colon| (atlas CpGs) {after:.4f}")
        normalized = True
    else:
        log(f"offset <= {OFFSET_LIMIT}: no normalization")
        normalized = False

    tum = tum.assign(normalized=normalized)
    T.astype("float32").to_parquet(out / "tumours.parquet")
    N.astype("float32").to_parquet(out / "normals.parquet")
    epic.loc[cand, need].astype("float32").to_parquet(out / "geo.parquet")
    full.loc[cand, BLOOD].astype("float32").to_parquet(out / "blood.parquet")
    offset.rename("offset").to_frame().to_parquet(out / "offset.parquet")
    tum[["case", "sample", "id", "CPE", "normalized"]].to_csv(out / "tumours.tsv", sep="\t",
                                                               index=False)
    (out / "atlas_cpgs.txt").write_text("\n".join(a_cpgs) + "\n")
    log_path.write_text("\n".join(lines) + "\n")
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
