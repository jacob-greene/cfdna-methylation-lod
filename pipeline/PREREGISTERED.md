# Phase 5 pre-registration: real-read validation

This file, `samples/phase5_samples.tsv` and `markers/` are committed before any SPOT-MAS
sample is scored. Later changes to them must be listed in the results with the reason.

## Terms

| Term | Meaning |
|---|---|
| Spike source | HCT116 colon cancer cell line WGBS, GSM8156741 (run SRR28401050), about 20M read pairs streamed from the start of the FASTQ pair |
| Threshold donor | One of 14 SPOT-MAS healthy donors that set every threshold and the background profile p0 |
| Eval donor | One of 16 other SPOT-MAS healthy donors; backgrounds of the spike-in tests. Never used to set a threshold |
| Pool of k | Reads of k donors combined; k = 1, 2, 4, 8 makes the depth axis |
| f | Share of fragments in a mixed sample that come from the spike source |
| Block | A 500 bp window centred on a marker CpG; overlapping windows of one set are merged |

## Samples (`samples/phase5_samples.tsv`, `scripts/phase5_select.py`, seed 20261010)

- 30 healthy controls (all donors female) and 10 early-stage breast cancers (stage I–IIB).
- The 16 eval donors form disjoint pools: 16 of k = 1, 8 of k = 2, 4 of k = 4, 2 of k = 8.
- Threshold pools of k: 50 random k-subsets of the 14 threshold donors (overlapping).

## Markers (`markers/`, `scripts/phase5_markers.py`)

| Set | CpGs | Blocks | Selection |
|---|---|---|---|
| atlas | 5,258 | 3,550 | Moss atlas CpGs of Phase 3b feature set (a) |
| tumour | 500 | 358 | Top 500 of TCGA-COAD tumour minus max(normal colon, Moss blood), max blood beta ≤ 0.2, on all 276 kept tumours (Phase 3b rule) |

HCT116 is CIMP-high, so it is more methylated at these markers than a typical colon tumour.
Detection may look better than for a typical patient.

## Processing (identical for every source, from FASTQ)

fastp: clip 15 nt from the R2 5′ end (Adaptase tail) and 15 nt from the R1 3′ end, no R2 3′
clip, then cap both reads at 85 nt so that the PE150 spike reads match the PE100 SPOT-MAS
reads. bwa-meth to hg38; duplicates marked, not removed; MethylDackel per-CpG counts and
M-bias; wgbstools `bam2pat` fragment patterns. Mixing is done on fragment patterns inside
the blocks, after duplicate removal.

## Scores (features are methylation calls at marker blocks only)

| Score | Rule |
|---|---|
| LLR, binomial | Sum over blocks of CpG-call LLR: p0 = pooled block methylation of the threshold donors; p1 = p0 + f × (TCGA tumour mean − p0) at the tested f. Both marker sets, reported separately |
| LLR, beta-binomial | Same, with each fragment's calls in a block as one beta-binomial draw; the intra-fragment correlation is estimated on threshold donors only |
| UXM-style | Tumour blocks: share of fragments with ≥ 3 CpG calls in the block that are all methylated |
| Transfer models | Phase 3b LR_L2, LR_L1 and XGB, trained on the Phase 3b D = 30 simulation (non-HI1 blanks, all 276 tumours, Phase 3b tuning). Real feature = block methylation of the probe's block; no call = training mean. Applied unchanged with the Phase 3b valid-rule threshold. AUC also reported |

## Thresholds and evaluation

- Untrained scores: one-sided 95 % normal prediction bound over threshold pools of the same
  k, df = max(floor(14 / k) − 1, 1).
- Spike-in grid: f ∈ {0, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1} × k ∈ {1, 2, 4, 8} × 5 seeds per
  eval pool. f = 0 rows give the held-out FPR.
- LoD: lowest f detected in ≥ 95 % of eval tests, interpolated on log f (Phase 3 rule).
- Binomial gap ladder at matched k: (1) Phase 3 binomial at the mean real calls per block;
  (2) binomial at each pool's real per-block calls (coverage dispersion); (3) real
  fragments (co-methylation and between-person variation); (4) real fragments with
  duplicates kept.
- Batch control: cross-validated logistic regression on fragment-level technical features
  (insert size, mapping quality, mismatches, read length) and on off-marker methylation,
  spike fragments vs threshold-donor fragments. AUC well above 0.5 is reported as a batch
  confound.
- External check: 10 cancers vs the 16 eval donors at k = 1, with threshold-donor
  thresholds. No tuning on these samples.
