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

Operator decision on issue 115 (comment 6100586554, a reply to comment 6100316928). This
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
