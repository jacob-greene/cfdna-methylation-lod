# Phase 6 pre-registration: optimised detectors

This file, `src/phase6.py`, `scripts/phase6_score.py` and `scripts/phase6_summary.py` are
committed before any eval donor, cancer or spike-in test is scored with a Phase 6 method.
`scripts/phase6_score.py --stage eval` refuses to run if any of them is uncommitted or
changed. Later changes must be listed in the results with the reason.

## Terms

| Term | Meaning |
|---|---|
| f | Spike fraction: the share of fragments from the spike source (HCT116 WGBS), counted in 2,000 off-marker windows (Phase 5 rule) |
| Effective share | The share of marker calls (WZ) or of ≥ 3-call fragments in the selected blocks (FRAG) that come from the spike |
| k | Pool size: the number of healthy donors pooled into one test sample |
| Threshold donor | One of the 14 SPOT-MAS healthy donors of Phase 5. The only data used to fit or choose anything in Phase 6 |
| Eval donor | One of the 16 other healthy donors; backgrounds of the spike-in tests |
| Block | A 500 bp window around a marker CpG (Phase 5 atlas and tumour sets) |
| Δ_b | TCGA-COAD tumour mean of block b minus p0_b, the pooled threshold-donor methylation |
| Tumour-like fragment | A fragment with ≥ 3 CpG calls in a block, all methylated if Δ_b > 0, all unmethylated if Δ_b < 0 |
| Leave-pool-out | Every fitted value is refitted without the donors of the threshold pool that is scored |
| LoD | Lowest f detected in ≥ 95 % of eval tests, interpolated on log f (Phase 5 rule) |
| FPR | Share of f = 0 eval tests above the threshold |

## Data and tests (unchanged from Phase 5)

Same samples, marker sets, block records (duplicates removed), eval pools, f grid
{0, 0.1 %, 0.3 %, 1 %, 3 %, 10 %}, 5 seeds and spike-mixing random seeds as
`scripts/phase5_score.py`. Each spiked test therefore holds the same spike fragments as in
Phase 5. Only methylation calls at marker blocks are used. No insert size, fragment length,
coverage or other technical feature enters any score.

## Scores (`src/phase6.py`)

| Score | Rule | Fixed parameters |
|---|---|---|
| WZ (O1) | Σ_b Δ_b (x_b − p0_b) / v_b ÷ √(Σ_b Δ_b² / v_b) over blocks with calls. x_b = block methylation of the test. v_b = q(1 − q) (G_b / N_b² + φ_b / k), q = p0 clipped to [0.01, 0.99], G_b = Σ_fragments n (1 + (n − 1) ρ), ρ = Phase 5 intra-fragment correlation | φ_b = per-block between-donor excess variance (method of moments, leave-one-donor-out p0), shrunk to the set mean with weight 14 donors, floored at 0; set mean for blocks with < 3 covered donors |
| FRAG (O2, O3) | Z = (C − E) / √(E + s² E² / k). C = tumour-like fragments in the selected blocks. E = Σ_b n3_b r_b, n3_b = ≥ 3-call fragments, r_b = background tumour-like rate | Selection (O3): abs(Δ_b) ≥ 0.3 and r_b ≤ 0.02. r_b = (c_b + 10 r̄) / (n3_b + 10), r̄ = set rate over blocks with abs(Δ_b) ≥ 0.3. s² = method-of-moments overdispersion, each donor scored with rates and selection refitted without it |

Both scores run on both marker sets. For a test, the model is fitted on all 14 threshold
donors. For a threshold pool, it is fitted leave-pool-out.

`binom_real` rung (positive control, WZ only): methylated calls ~ Binomial(N_b,
(1 − f) p0_b + f h_b) at each test's real calls per block, h_b = spike block methylation,
scored with the same weights. Thresholds from the same draws over the threshold pools.

## Thresholds

| Rule | Scores | Definition |
|---|---|---|
| bound (primary) | WZ, FRAG | Phase 5 rule: one-sided 95 % normal prediction bound over the threshold pools of the same k (14 single donors at k = 1, 50 random k-subsets otherwise), df = max(floor(14 / k) − 1, 1) |
| model (O4, secondary) | FRAG | Z > 1.645 |

Target FPR: 5 %.

## Endpoints

| Endpoint | Definition |
|---|---|
| Primary | LoD (nominal f) of FRAG on the tumour set, bound rule, at k = 1 and k = 2. Baseline: Phase 5 UXM (9.1 % at k = 1, 7.4 % at k = 2) |
| Secondary | FRAG atlas; WZ atlas and tumour (baseline Phase 5 LLR_binom); LoD on effective share; O4 FPR and LoD |
| FPR | Per score and k, with a Clopper-Pearson 95 % interval. Resolution 1/16 at k = 1, 1/8 at k = 2 |
| Matched FPR (descriptive) | For every Phase 5 `real` score and every Phase 6 score: detection with the threshold at the highest eval f = 0 score of the same k, the resulting LoD, and AUC |
| k = 4 and 8 | Reported as direction only: 4 and 2 independent eval backgrounds |
| External check | The 10 breast cancers vs the 16 eval donors at k = 1, bound thresholds. The markers are colon markers, so no signal is expected |
| Confound check | Mean calls per fragment record and share of records with ≥ 3 calls, spike vs threshold donors, all blocks and FRAG blocks |

A `real` WZ detection rate far above its `binom_real` rate is treated as a bug until shown
otherwise.

## Design-stage results (threshold donors only; no eval donor, cancer or spike was read)

`scripts/phase6_score.py --stage threshold` wrote `results/phase6/threshold_scores.tsv`
and `results/phase6/model_info.tsv` before this file was committed.

| Item | atlas | tumour |
|---|---|---|
| Blocks with p0 | 3,549 | 358 |
| ρ (Phase 5 value, refitted) | 0.174 | 0.082 |
| φ set mean | 0.130 | 0.066 |
| FRAG blocks selected (leave-one-out range) | 412 (409–423) | 348 (347–349) |
| FRAG blocks with Δ_b > 0 | 332 | 348 |
| FRAG mean background rate, selected blocks | 0.0060 | 0.0010 |
| FRAG s² (leave-one-out range) | 0.18 (0.14–0.25) | 3.2 (1.2–5.9) |

Design check of O1 (unverified prediction for the spike): on leave-one-out threshold
donors, the null SD of three weightings was compared with the shift that a TCGA-direction
signal (p0 + f Δ) would give. Smaller f is better.

| Weighting | f for a 2.5-SD shift, atlas | Same, tumour |
|---|---|---|
| Binomial (Phase 5 LLR_binom at small f) | 0.020 | 0.014 |
| Binomial + co-methylation | 0.021 | 0.014 |
| WZ (+ between-donor φ_b) | 0.023 | 0.016 |

So per-block weighting does not reduce the null spread: the excess spread is mostly a
donor-level shift shared by many blocks, not a few noisy blocks. WZ is still scored once,
as pre-registered, as a test of O1. The prediction is no gain over LLR_binom.

Other design-stage observations:

- Leave-pool-out WZ null scores are not centred on 0 (atlas mean 3.2 at k = 1, 17.5 at
  k = 8). A plug-in p0 with a p0-dependent variance gives a positive bias, and the bias
  grows when fewer donors fit the model (14 − k). Eval tests use a model from 14 donors, so
  the threshold pools may give conservative thresholds at large k.
- The FRAG tumour s² comes mostly from one threshold donor (9 tumour-like fragments vs 1.2
  expected). The same donor has the highest WZ and FRAG atlas null scores. It stays in.
- Shares of threshold pools above Z = 1.645 (O4): FRAG atlas 0.07, 0.02, 0, 0.02 and FRAG
  tumour 0.07, 0.16, 0.16, 0.32 at k = 1, 2, 4, 8. So O4 is expected to exceed the 5 % target
  for FRAG tumour at k ≥ 2.

## Caveats stated in advance

- HCT116 is CIMP-high, so it is more methylated at these markers than a typical colon tumour.
  Detection may look better than for a typical patient.
- The spike is a cell-line library with longer inserts (Phase 5). Fragments with more calls
  pass the ≥ 3-call rule more often. The confound check reports this.
- The effective share is below nominal f (Phase 5: atlas 0.039, tumour 0.018 at f = 0.1).
  LoD is reported on both scales.

## Amendment 1 (2026-10-10): effective share primary, length-matched spike rung

Decision of the project lead, after a survey of public spike sources. This
amendment is committed before any eval donor, cancer or spike-in test is scored with a
Phase 6 method. Nothing above is removed; where this section differs, this section applies.
The GSE186007 plasma arm is pre-registered in a later, separate amendment, before its first
methylation scoring.

### A1. Effective share is the primary scale of f

| Endpoint | Definition |
|---|---|
| Primary | LoD of FRAG on the tumour set, bound rule, rung `real`, k = 1 and k = 2, **on the effective share** (mean share of ≥ 3-call fragments in the FRAG-selected blocks that come from the spike, per k and f). LoD on nominal f is reported next to it |
| Scale per score | FRAG: `share_frag_n3_selected`. WZ: `share_calls` (spike share of calls in blocks with p0) |
| Mapping | Nominal-f LoD mapped through the mean effective share per (rung, set, k, f), interpolated on log-log axes (`lod_effective` in `scripts/phase6_summary.py`) |
| Phase 5 baselines on the same scale | UXM through `share_frag_n3_selected`, LLR_binom through `share_calls`, of the same set and k, rung `real`. The tests are the same, so the ranking of scores is the same on both scales |

### A2. Length-matched spike rung (`len`)

The spike library has longer inserts than SPOT-MAS cfDNA. In-silico trimming of the reads
cannot add the unsequenced middle of a long fragment, so a trimmed fragment would carry fewer
calls than a real short fragment. This rung therefore **selects** HCT116 fragments so that
their insert sizes match the SPOT-MAS distribution. It changes the spike only. No score uses
insert size or any other technical feature.

| Item | Rule (`scripts/phase6_lenmatch.py`, `scripts/phase6_lenmatch.sbatch`) |
|---|---|
| Reference distribution | Insert sizes of the 14 threshold donors (read 1, off-marker blocks, MAPQ ≥ 10, duplicates excluded) |
| Bins | 25 equal-count bins of the reference; outer edges 0 and 1000 bp |
| Keep probability | (1/25 ÷ HCT116 share in the bin), divided by its maximum over bins |
| HCT116 fragments | Read names with a read in a marker or off-marker block, 0 < abs(TLEN) < 1000 bp; others dropped |
| Random draws | `default_rng([SEED, 70])` over read names in sorted order |
| Then | Filtered BAM → `bam2pat` (wgbstools 0.3.0, flags 1796, as in `pipeline/main.nf`) → records cache `SRR28401050_len.dedup.npz` |
| Mixing | Same pools, f grid and 5 seeds as rung `real`; random seeds `[SEED, 1, k, pid, fi, seed, 1]` (rung `real` ends in 0, unchanged) |
| Dosing | Phase 5 rule with the length-matched spike's own off-marker fragment count |
| Unscored tests | If the keep probability q = f/(1 − f) × n_bg / n_spike exceeds 1, the test is not scored (NaN). A (k, f) cell with any unscored test is dropped. LoD is NaN when no scored f reaches 95 % |
| Null tests and thresholds | Shared with rung `real` (the f = 0 tests hold no spike) |
| Scores | WZ and FRAG, both sets. No `binom_real` for this rung |

Design data of the selection (spike and threshold donors only; no score was computed):

| Source | Fragments | Insert q5 | q25 | Median | q75 | q95 |
|---|---|---|---|---|---|---|
| Threshold donors (reference) | 65,519 | 120 | 146 | 159 | 173 | 207 |
| HCT116, all (in blocks) | 20,402 | 127 | 174 | 213 | 260 | 350 |
| HCT116, kept | 2,426 | 120 | 145 | 159 | 172 | 204 |

| Count | Value |
|---|---|
| Off-marker fragments, HCT116 all / kept | 5,098 / 604 |
| Off-marker fragments per threshold donor, median (max) | 4,136 (5,393) |
| Expected unscored cells (q > 1), from threshold-donor counts | f = 10 % at k ≥ 2; f = 3 % at k = 8; possibly f = 3 % at k = 4 |

Caveats stated in advance for this rung:

- Only 2,426 spike fragments remain. At high f the 5 seeds reuse most of them, so the
  spread over seeds understates the true sampling spread.
- Selection on length can also select on methylation if fragment length depends on it in the
  HCT116 library. The confound table reports calls per record for `spike_len` as well.
- The rung is a sensitivity check of the length confound, not a new primary endpoint.

### A3. Files under the eval guard

`scripts/phase6_lenmatch.py` and `scripts/phase6_lenmatch.sbatch` join the files that
`scripts/phase6_score.py --stage eval` requires to be committed and unchanged.

### A4. Bug fix before the first eval scores (2026-10-10)

The first eval run (Slurm job 10356043) stopped at its first test with an error and wrote
no score. In the `binom_real` rung, a block with no p0 and no spike calls got a NaN spike
methylation h_b, and the binomial draw refused it. Such blocks are excluded from WZ, so h_b
is now 0 there. No rule changes. A dry run with threshold donors as stand-in backgrounds
checks the eval code path before the one real eval run.


## Amendment 2 (2026-10-12): GSE186007 plasma arm

Decisions of the project lead: add the plasma arm, and never keep raw FASTQ after a
successful alignment (disk budget). This
amendment is committed before any GSE186007 run is scored with a methylation score. Up to this
commit, the GSE186007 data were used only for alignment QC (`phase6_gse/qc_table.tsv`) and
for the tumour fraction (B2), which uses read depth only. Nothing above is removed.

### Terms added

| Term | Meaning |
|---|---|
| LABS | The library protocol of GSE186007 (linear-amplification bisulfite sequencing of plasma cfDNA; doi:10.1186/s13059-024-03262-2). A different protocol from SPOT-MAS |
| CRC run, HEA run | A GSE186007 plasma run from a colorectal cancer patient or a healthy person |
| TF | Tumour fraction: the share of a run's cfDNA that comes from the tumour, estimated from copy number (B2) |
| f_tumour | f × TF: the share of fragments in a test that come from tumour DNA |
| Healthy-LABS control | The same mixing with an HEA run as the spike. It measures how much of the detection comes from the protocol and not from tumour DNA |

### B1. Data (fixed from metadata before any download)

| Item | Rule |
|---|---|
| Runs | The first 8 CRC and the first 4 HEA paired-end NovaSeq runs by run accession in `cfdna-methylation-lod-data/sra/gse186007_runs.tsv`: CRC SRR24725513–SRR24725520, HEA SRR24725563–SRR24725566 |
| Reads | The first 20 M read pairs of each run, streamed from ENA (`--stream_runs`) |
| Processing | `pipeline/main.nf` with the Phase 5 parameters unchanged (the fastp clips and the 85 bp cap included), outdir `cfdna-methylation-lod-data/phase6_gse` |
| FASTQ | Downloaded and trimmed FASTQ deleted after the BAM check passes; Nextflow work files deleted after each batch |

### B2. Tumour fraction (`scripts/phase6_tf.py`)

Read depth only; no methylation call is used. The method is a simplified ichorCNA
(Adalsteinsson et al. 2017, PMID 29109393). The full rules are in the script docstring:

| Step | Rule |
|---|---|
| Reads | Read 1, primary, mapped with mapped mate, not duplicate, MAPQ ≥ 30, chr1–22; reads that start in a marker block are skipped |
| Bins | 1 Mb; kept if N share < 0.1 and the median HEA count > 0 |
| GC correction | Count share divided by the run's median over 20 equal-count GC bins |
| Panel of normals | Median of the 4 HEA runs; an HEA run uses the other 3 |
| Arms | Split at the centromere midpoint (cytoBand acen); acrocentric p arms dropped; ≥ 10 bins |
| Altered arm | abs(r_a / SE_a) > 4, with SE_a inflated by c, the SD of the HEA arm z-scores (leave-one-out) |
| TF | Median over altered arms of 2 (2^r − 1) (gains) or 2 (1 − 2^r) (losses), capped at 1. TF = 0 when no arm is altered |
| Interval | 500 bootstrap resamples of bins within arms, `default_rng([SEED, 80])` |

### B3. Mixing design

| Item | Rule (`scripts/phase6_gse_score.py`) |
|---|---|
| Spike | Each of the 12 GSE186007 runs alone, duplicates removed (`dedup` pat) |
| Backgrounds, f grid, seeds | The eval pools, f grid and 5 seeds of rung `real` |
| Dosing | Phase 5 rule with the spike run's own off-marker fragment count. A test with keep probability q > 1 is not scored, as in rung `len` |
| Random seeds | `default_rng([SEED, 2, k, pool, fi, seed, i])`, i = index of the run in `phase6_gse/runs.tsv` |
| f = 0 tests | Recomputed; they must equal the rung `real` f = 0 scores of `results/phase6` (the summary checks this) |
| Tumour scale | f_tumour = f × TF of the spike run |
| Direct scores (descriptive) | Each GSE186007 run scored alone at k = 1 |

### B4. Scores and thresholds

Unchanged: WZ and FRAG on both sets, with the model fitted on the 14 threshold donors. The
thresholds are the rung `real` thresholds in `results/phase6/thresholds.tsv` (bound rule
primary, model rule for FRAG secondary).

### B5. Endpoints (`scripts/phase6_gse_summary.py`)

| Endpoint | Definition |
|---|---|
| Primary | FRAG, tumour set, bound rule, k = 1 and k = 2: the median over CRC runs with TF > 0 of LoD_tumour = LoD (nominal f of that run) × TF. Reported with its range and the number of runs with a finite LoD |
| Healthy-LABS control (reported next to every CRC result) | The detection rate over all HEA spike tests per (score, set, k, rule, f). If it exceeds 0.2 at a (k, f), CRC detection there is flagged "protocol" and is not attributed to tumour DNA. With no protocol effect, the rate is near the FPR |
| Secondary | WZ and the atlas set; LoD on nominal f per run; LoD on the tumour effective share (the run's effective share at the LoD × TF); pooled CRC and HEA detection rates per f; direct scores (CRC vs HEA, k = 1) |
| Comparison to HCT116 (descriptive) | LoD_tumour and the tumour effective share next to the HCT116 FRAG tumour LoD (2.15 % and 2.01 % effective share at k = 1 and 2) |
| CRC runs with TF = 0 | LoD on nominal f only; their TF upper bound is reported |

### B6. Caveats stated in advance

- TF comes from 20 M read pairs and a panel of only 4 normals. Its interval is wide, and a
  run with a small TF may show no altered arm.
- The single-copy model gives a biased TF for arms with more than one copy changed or with
  subclonal changes. The median over arms limits this.
- An LABS spike also carries the patient's non-tumour cfDNA and the LABS protocol signature.
  The healthy-LABS control measures this part.
- At high f the 5 seeds reuse most spike fragments, so the spread over seeds understates the
  sampling spread.
- 8 CRC runs, with no stage or TF in the GEO metadata.

### B7. Files under the guard

`scripts/phase6_gse_score.py` refuses to run unless this file, `src/phase6.py`,
`scripts/phase6_score.py`, `scripts/phase6_summary.py`, `scripts/phase6_tf.py`,
`scripts/phase6_gse_score.py` and `scripts/phase6_gse_summary.py` are committed and unchanged.

### B8. Tumour fraction and QC, fixed before scoring

The TF table below was computed with the B2 rules before any GSE186007 run was scored
(`phase6_gse/tf/tf.tsv`; 2,711 bins, 39 arms, inflation c = 1.742). TF depends on read depth
only. It is now fixed for the endpoints in B5.

| Run | Group | TF | Altered arms | 95 % interval | Upper bound if no arm altered |
|---|---|---|---|---|---|
| SRR24725513 | CRC | 0 | 0 | 0–0.120 | 0.065 |
| SRR24725514 | CRC | 0 | 0 | 0–0.103 | 0.055 |
| SRR24725515 | CRC | 0 | 0 | 0–0.048 | 0.062 |
| SRR24725516 | CRC | 0 | 0 | 0–0.060 | 0.059 |
| SRR24725517 | CRC | 0.599 | 18 | 0.566–0.617 | – |
| SRR24725518 | CRC | 0 | 0 | 0–0.072 | 0.052 |
| SRR24725519 | CRC | 0.409 | 22 | 0.369–0.427 | – |
| SRR24725520 | CRC | 0.099 | 1 (19q gain) | 0–0.135 | – |
| SRR24725563 | HEA | 0 | 0 | 0–0 | 0.058 |
| SRR24725564 | HEA | 0 | 0 | 0–0 | 0.060 |
| SRR24725565 | HEA | 0 | 0 | 0–0 | 0.060 |
| SRR24725566 | HEA | 0 | 0 | 0–0.123 | 0.095 |

So the primary endpoint (B5, CRC runs with TF > 0) uses 3 runs: SRR24725517, SRR24725519
and SRR24725520. The TF of SRR24725520 rests on one arm, and its interval includes 0.

Processing record (no rule change):

| Item | Record |
|---|---|
| Alignment QC | `phase6_gse/qc_table_all12.tsv`, one pipeline pass over all 12 runs. All runs: 20 M pairs, mapped share ≥ 0.979, CHH conversion ≥ 0.985 |
| M-bias | R1 CpG methylation is higher at positions 1–4 and 82–85 (about 78–85 %) than mid-read (about 74 %). The R1 range is 16–27 %, against about 8 % in Phase 5. R2 is flat. The Phase 5 clips stay unchanged (B1) |
| Pipeline failures | Two transient failures, both rerun: a bam2pat temp-dir cleanup error on NFS (Errno 39), and a broken ENA stream for SRR24725565 caught by the line-count check. `BAM2PAT` and `FETCH_STREAM` now retry twice. A rerun starts each run from its published BAM; no run was fetched twice after alignment |
| FASTQ | 48 FASTQ files (12 runs × 4) deleted after the BAM check (`phase6_gse/fastq_cleanup.tsv`) |
