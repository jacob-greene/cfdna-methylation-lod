# cfdna-methylation-lod

How low a tumour fraction can a cfDNA methylation test detect? This repository restates the
plasma-dilution result of Moss et al. 2018 (Nat Commun 9:5068, Fig. 7d) as a clinical-style
limit of blank and limit of detection. It then simulates sequencing depth and compares
untrained scores with logistic regression and XGBoost.

All results are in silico. They use public array data only. Validation on real sequencing
reads is not done yet (see [Status](#status)).

## Terms

| Term | Meaning |
|---|---|
| Blank | A healthy cfDNA sample with no tumour DNA. There are 12 EPIC arrays: 8 pools and 4 individuals, HI1–HI4. |
| f | The tumour-cell DNA fraction of a simulated sample. |
| D, sequencing depth | The number of reads that cover one CpG. |
| LoB, limit of blank | The 95th percentile of a score over blanks (CLSI EP17-A2). |
| LoD, limit of detection | The lowest f detected in at least 95 % of tests (CLSI EP17), interpolated on log f. |
| FPR, false-positive rate | The share of held-out blank results at f = 0 scored above the threshold. EP17 needs 5 % or less. |
| NNLS colon fraction | The colon-epithelium share from non-negative least-squares deconvolution against the Moss et al. 2018 atlas. |
| LLR, log-likelihood ratio | A score that sums, over marker CpGs, how much better "tumour at fraction f" explains the read counts than "background only". |
| Valid rule | Each training blank is scored by a model built without its group. The threshold is a one-sided 95 % normal prediction bound over people. |
| LR_L2, LR_L1 | Logistic regression with an L2 or L1 penalty, on standardized features. |
| XGB | XGBoost gradient-boosted trees. |
| Atlas features | The 5,258 Moss atlas CpGs present in all inputs. |
| Seed | The random seed. It sets the fold assignment, the training subsample and the LR_L1 solver. |

## Question

Moss et al. report that the colon signal of a colon-cancer plasma is detected "down to 3 %"
after dilution into healthy plasma (r = 0.92). This repository asks four questions:

1. Do the public GEO arrays reproduce Fig. 7d?
2. What are the LoB and the LoD of that measurement under CLSI EP17?
3. Does deeper sequencing at the same marker CpGs move the LoD toward MRD-level fractions
   (0.01–0.1 %)?
4. Do trained classifiers (logistic regression, XGBoost) beat the untrained NNLS and LLR
   scores?

## Data provenance

Only openly downloadable data is used. Downloads go to `$CFDNA_DATA` (default:
`../cfdna-methylation-lod-data`), never into the repository. `data/` holds the download
scripts and `data/checksums.sha256`.

| Source | What | Script |
|---|---|---|
| GEO GSE122126 (Moss et al. 2018) | Processed betas (series matrices): 12 healthy cfDNA blanks, the colon-cancer plasmas 806, 684 and 2010, the Fig. 7d mixes, cell-type replicates | `data/download_geo.sh` |
| `nloyfer/meth_atlas` (GitHub) | Reference atlas, 7,890 CpGs × 25 cell types; NNLS code | `data/download_atlas.sh` |
| Moss et al. 2018 supplement | Supplementary tables | `data/download_moss_supp.sh` |
| Illumina manifests (hg38) | 450K and EPIC probe annotation | `data/download_manifests.sh` |
| GDC, TCGA-COAD | 352 open-access 450K beta files (SeSAMe level 3) | `data/download_tcga.sh` |
| Aran et al. 2015 | Consensus purity estimates (CPE) for TCGA tumours | `data/download_purity.sh` |
| SRA PRJNA929650 (SPOT-MAS) | Run metadata and 2 runs for bisulfite QC (Phase 0) | `data/download_sra_metadata.sh`, `data/download_spotmas_*.sh` |

Run the scripts from the repository root, for example `bash data/download_geo.sh`. Each
script skips files that already exist. Verify with
`cd "$CFDNA_DATA" && sha256sum -c <repo>/data/checksums.sha256`.

## Results in plain language

### Phases 1–3: the array result and the depth simulation

| Phase | Result |
|---|---|
| 1, reproduce Fig. 7d | All 8 points match the published figure within 0.15 percentage points. Our r is 0.940; the paper states 0.92. |
| 2, array LoB and LoD | 95 % detection needs about 2.0–2.3 % colon DNA (about 3.0–3.4 % pl806 plasma). The paper's "down to 3 %" sits just below the LoD. |
| 3, depth simulation | More reads do not lower the LoD over real blanks. It stays at 1.2–2.7 % from D = 100 up. Between-person variation sets the floor. |

![LoD against sequencing depth, by score and marker count](figures/phase3_depth_lod.png)

Figure 1. LoD (95 % detection) against depth D. Top row: HI1 included. Bottom row: HI1
removed. Left: plasma spike. Middle: colon-cell spike. Right: LLR on all markers with an
input cap, pooled rule (superseded; the capped runs exist only for that rule). Open circles:
FPR above 5 %, so that LoD is not valid. Dotted and "control" lines: no between-person
variation. Without between-person variation, the LoD falls with depth. With real blanks, it
stays flat.

### Phase 3b: trained classifiers against untrained scores

Tumour DNA comes from 276 TCGA-COAD primary tumours (CPE ≥ 0.6), spiked into held-out
blanks. No person and no tumour is in both a training and a test fold. Thresholds use the
valid rule. HI1 is excluded from the headline; with HI1 included, no model meets the 5 % FPR
rule.

| Model, track | LoD, % tumour-cell DNA | Held-out false positives | Label |
|---|---|---|---|
| LR_L1, array, atlas features | 1.6 at the fixed seed; 0.9–1.7 across 7 seeds (median 1.5)¹ | 0 of 11 people | unverified |
| NNLS, array and reads from D = 100 up | 1.9–2.0 | within the 5 % rule (point estimate) | array: 1.88–1.92 across 7 seeds¹ |
| Trained models, reads from D = 100 up, atlas features | 2.2–3.4 | per row in `results/phase3b/lod.tsv` | — |

¹ The 7-seed ranges come from an independent re-run with a modified seed. The scripts here
use one fixed seed (`SEED` in `scripts/phase3b_classify.py`), so they reproduce only the
fixed-seed values.

- LR_L1 is below NNLS in 7 of 7 seeds, so the ranking holds. The value does not: it moves
  with the seed. It is not a floor.
- XGBoost does not beat logistic regression consistently. Within folds, the largest AUC gap
  is +0.07.
- Where XGBoost gains on tumour-specific features, part of the gain is colon-tissue signal,
  not cancer signal (TCGA normal-colon spikes are called too).

![Sensitivity against tumour fraction per model, feature set and track](figures/phase3b_classifiers.png)

Figure 2. Sensitivity at the valid-rule threshold against f, HI1 excluded. Top row: atlas
features. Bottom row: tumour features selected inside each training fold. Columns: array
betas, then simulated reads at D = 30, 100 and 300. Bands: 95 % group-bootstrap CI. Dashed
lines with open markers: held-out FPR above 5 %. Dotted line: 95 %.

The full tables are in `results/phase3b/lod.tsv`, `fpr.tsv`, `auc.tsv` and the control
tables beside them.

## Limitations

- **Every FPR rests on 11 people in 7 groups.** "Passes the 5 % rule" is a point estimate,
  not proof. On the array track, 0 of 11 false positives has a one-sided 95 % upper bound of
  24 %.
- **The LR_L1 LoD depends on the seed** (0.9–1.7 % over 7 seeds, from an independent re-run
  not scripted here) and is unverified. Its
  tuned penalty varies from fold to fold.
- **HI1 is an outlier blank.** It is 4.4 SD above the other 11 blanks in Phase 3. With it,
  no model meets the 5 % FPR rule. Results are given with and without it.
- **The read model is simple.** Reads are Binomial(D, beta): independent CpGs, no
  fragment-level co-methylation, no bisulfite loss or PCR duplicates, array betas taken as
  molecular rates. Real reads and a beta-binomial model are Phase 5 work.
- **f is not the same quantity in every phase.** Phases 2–3 spike the pl806 colon fraction or
  colon cells; Phase 3b spikes the tumour-cell fraction of TCGA tumours.
- The f grid is coarse near the LoD (1, 2, 3 %), so an interpolated LoD can be off by about
  0.3 points.
- TCGA betas (SeSAMe) and GEO betas (minfi) are processed differently. The median absolute
  offset at the atlas CpGs is 0.028, below the 0.05 limit set in advance.

## Layout

| Path | Contents |
|---|---|
| `data/` | Download scripts and checksums only |
| `src/deconv.py` | Deconvolution: NNLS against the atlas |
| `src/lod.py` | LoB and LoD (CLSI EP17-A2) |
| `src/simulate.py` | Read-count simulation, LLR, prediction bounds |
| `src/classify.py` | Tumour spikes, grouped folds, fold-internal feature selection, AUC |
| `src/atlas.py`, `src/geo.py`, `src/tcga.py`, `src/spotmas.py` | Data readers |
| `scripts/` | One script per analysis step (`phase*_*.py`) and Slurm wrappers (`*.sbatch`) |
| `results/` | Summary tables per phase (raw per-run scores are not tracked) |
| `figures/` | Figures per phase |
| `notebooks/01_walkthrough.ipynb` | Walk-through of the results from the summary tables |
| `tests/` | pytest suite |
| `pipeline/` | Reserved for the Phase 5 Nextflow pipeline |

Quality control has no module of its own. It lives in `scripts/phase0_spotmas_qc*` (bisulfite
conversion, M-bias) and `scripts/phase1_fig7d.py` (array missing rates).

## Reproduce

```bash
uv sync                      # Python 3.11, pinned by uv.lock
bash data/download_geo.sh    # and the other data/download_*.sh scripts
uv run pytest -q
uv run ruff check .
uv run python scripts/phase1_fig7d.py
uv run python scripts/phase2_lob_lod.py
REPO_DIR=$PWD CFDNA_DATA=<data dir> sbatch --array=0-42 scripts/phase3_depth_sim.sbatch
uv run python scripts/phase3_summary.py
REPO_DIR=$PWD CFDNA_DATA=<data dir> STEP=prepare sbatch scripts/phase3b.sbatch
REPO_DIR=$PWD CFDNA_DATA=<data dir> sbatch --array=0-15 scripts/phase3b.sbatch
uv run python scripts/phase3b_summary.py
uv run --with jupyter jupyter lab notebooks/01_walkthrough.ipynb   # reads results/ only
```

Every simulation and model uses a fixed seed (see the `SEED` constant of each script).

## Tests

| Test | File |
|---|---|
| NNLS output is non-negative and sums to 1 | `tests/test_deconv.py` |
| NNLS recovers known proportions on synthetic mixes | `tests/test_deconv.py`, `tests/test_lod.py` |
| Simulations are reproducible under a fixed seed | `tests/test_simulate.py` |
| LoB uses only blank rows of the LoB set | `tests/test_lod.py` |
| No person or tumour is in both train and test folds, at the outer and inner level | `tests/test_leakage.py`, `tests/test_classify.py` |
| Fold-internal feature selection never reads a test tumour | `tests/test_leakage.py`, `tests/test_classify.py` |

The leakage tests were checked by mutation: each deliberate leak in the fold code makes them
fail.

## Status

| Phase | State |
|---|---|
| 0, data audit | done |
| 1, reproduce Fig. 7d | done |
| 2, array LoB and LoD | done |
| 3, depth simulation | done |
| 3b, classifier comparison | done |
| 4, repository polish | this version |
| 5, real sequencing reads (read-level in-silico mixing, Nextflow) | not started |

## Reference

Moss J, Magenheim J, Neiman D, et al. Comprehensive human cell-type methylation atlas reveals
origins of circulating cell-free DNA in health and disease. Nat Commun 9, 5068 (2018).
