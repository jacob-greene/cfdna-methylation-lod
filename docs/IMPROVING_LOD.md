# How the LoD could improve

This page explains what limits the limit of detection (LoD) in this project, and which data
would improve it. Results and all terms: [`RESULTS.md`](RESULTS.md#terms).

- f: the spike-in fraction, the share of all fragments that come from the spike-in.
- k: the pool size, the number of healthy donors whose reads form one test background
  (k = 1 is one person's plasma).
- Tumour DNA at the markers: the share of informative fragments at the marker blocks that come
  from tumour DNA. All LoD values use this scale.
- FPR: false-positive rate, the share of spike-free tests called positive.

## What limits the LoD

At low spike-in levels, the number of tumour-like spike-in fragments
in a test varies widely from test to test (Figure 1c in the [README](../README.md), error bars). For one HCT116 spike-in in the
16 test donors at k = 1, with 80 tests per level (SD: standard deviation over the 80 tests):

| Tumour DNA at the markers | Mean tumour-like spike-in fragments per test | SD | Tests with none |
|---|---|---|---|
| 0.01 % | 0.09 | 0.28 | 91 % |
| 0.07 % | 0.53 | 0.95 | 65 % |
| 0.21 % | 1.85 | 1.30 | 12.5 % |
| 0.63 % | 5.6 | 3.1 | 2.5 % |
| 2.2 % | 20.2 | 7.8 | 0 % |

Each one-person test holds about 1,600 informative fragments at the tumour markers. At the
LoD, about 20 are tumour-like spike-in fragments. Background differences between people dominate
the score's noise at this depth. At the levels needed to monitor minimal residual disease (MRD, 0.01–0.1 %), a test would hold about one tumour fragment or
fewer, and most tests hold none, so MRD needs many more informative fragments, not only a
better score. To improve the LoD, model gains must come from more healthy donors to set
thresholds. Such data can also train a score that models co-methylation within fragments and
variation between people. It can also find markers that are consistent across tumour types,
for panel design.

## Data sets and what would help

All data are public. Each source fills one role in the spike-in design.

| Source | Role | Why this source | Limit for the model |
|---|---|---|---|
| SPOT-MAS cfDNA (SRA PRJNA929650) | Healthy backgrounds | 30 healthy donors with whole-genome bisulfite reads: enough to keep threshold donors (14) apart from test donors (16) | About 15 million read pairs per donor (median; range 8–26 million). Only 16 spike-free tests at k = 1 and 8 at k = 2: an FPR of 0 of 16 has a 95 % upper bound of 21 % |
| HCT116 cell line (GSM8156741) | Tumour spike-in | A colon-cancer genome with whole-genome bisulfite reads and no normal-cell DNA | One cell line. HCT116 is CIMP-high, so it is more methylated at the markers than a typical tumour |
| Colorectal-cancer plasma (GSE186007) | Patient spike-in | Real patient cfDNA, with healthy runs from the same protocol as a control | 8 cancer runs; only 3 have tumour DNA measurable by copy number. A different library protocol from the backgrounds |
| TCGA-COAD arrays (GDC) | Marker selection | 276 tumours with purity of 0.6 or more | Markers come from tumour tissue. The patient run with the highest tumour fraction is rarely detected |
| Moss et al. 2018 arrays and atlas (GSE122126) | Idealised reference | A published plasma-dilution series to reproduce | 12 healthy blanks; array, not sequencing, data |

A better data set for this model would have:

| Change | Reason, from the results |
|---|---|
| More healthy donors, not more depth per donor | Between-person variation sets the floor; depth beyond 100 reads per CpG does not move the LoD in simulation |
| More cancer plasmas with an independent tumour fraction (for example tumour-informed) | Copy number finds tumour DNA in 3 of 8 runs only. Runs without copy-number changes cannot be placed on the tumour-DNA scale |
| Backgrounds and patient spikes from the same library protocol | Scored on their own, all 4 healthy runs of the other protocol are called on the atlas marker set, so that set does not transfer |
| More f grid points between 3 % and 10 % | Every finite LoD on real reads lies in that interval, so the grid cannot rank scores within it |
| More than one cell line | One CIMP-high line may overstate detection of typical tumours |

