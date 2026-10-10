# Phase 5 pipeline: real-read validation

This Nextflow pipeline turns FASTQ files into fragment-level methylation patterns at the
marker blocks. Every source goes through the same steps: the SPOT-MAS cfDNA runs and the
spike source (HCT116 WGBS). The processing rules are fixed in `PREREGISTERED.md`.

## Terms

| Term | Meaning |
|---|---|
| Spike source | HCT116 colon cancer cell line WGBS, GSM8156741 (run SRR28401050); the first 20M read pairs are streamed |
| Block | A 500 bp window centred on a marker CpG; overlapping windows of one set are merged |
| Off-marker block | One of 2,000 random 500 bp autosomal windows at least 10 kb from every marker block |
| pat file | wgbstools fragment pattern file: one row per fragment, with its first CpG index and its C/T calls |

## Steps

| Process | Does |
|---|---|
| `FETCH_ENA` | Downloads a SPOT-MAS FASTQ pair from ENA and checks the MD5 sums |
| `FETCH_STREAM` | Streams the first N read pairs of the spike source |
| `FASTP` | Clips 15 nt from the R2 5′ end and the R1 3′ end, then caps reads at 85 nt |
| `BWAMETH_INDEX` | Builds the bwa-meth hg38 index once (`storeDir`) |
| `BWAMETH_ALIGN` | Aligns, fixes mate tags, sorts, marks duplicates (not removed) |
| `METHYLDACKEL` | Per-CpG counts, M-bias table, CHH calls on chr21 (conversion rate) |
| `BAM2PAT` | Fragment patterns in marker and off-marker blocks, duplicates removed and kept |
| `FRAG_FEATURES` | Insert size, MAPQ, mismatches and read length in off-marker blocks |
| `BLOCK_INDEX` | Block coordinates to CpG index ranges |
| `QC_TABLE` | One QC row per run |

## Profiles

| Profile | Use |
|---|---|
| `slurm` | Slurm executor; partition from `NXF_SLURM_QUEUE`; at most 8 jobs at a time |
| `conda` | Pinned micromamba environment (`envs/phase5.explicit.txt`) with the container tool versions |
| `apptainer` | Biocontainers through Apptainer; the wgbstools image from `containers/wgbstools/` |
| `docker` | The same containers through Docker (used by the CI smoke test) |
| `awsbatch` | AWS Batch executor; written for portability, not run in this project |
| `test` | Smoke test on a small bundled BAM (`tests/data/`) |

The HPC run used `slurm,conda`. Apptainer cannot start containers inside the agent sandbox
on this cluster, also not in Slurm jobs (no setuid starter, no mount propagation). The
`conda` environment pins the same tool versions as the containers:

| Tool | Version |
|---|---|
| fastp | 0.23.4 |
| bwa-meth | 0.2.10 |
| bwa | 0.7.18 |
| MethylDackel | 0.6.1 |
| samtools | 1.21 |
| wgbstools | 0.3.0 (built from the GitHub release tarball) |

## Run

```bash
export CFDNA_DATA=/path/to/cfdna-methylation-lod-data NXF_SLURM_QUEUE=campus-new
cd $CFDNA_DATA/phase5/nxf
nextflow run /path/to/cfdna-methylation-lod/pipeline/main.nf -profile slurm,conda -resume
```

- `--only_runs SRR1,SRR2` restricts the SPOT-MAS runs (tests).
- `--spike_run false` (or `''`) skips the spike source.
- Run the Nextflow head process inside a Slurm job; it runs for hours.

Smoke test (Docker, as in CI):

```bash
docker build -t cfdna-wgbstools:0.3.0 pipeline/containers/wgbstools
nextflow run pipeline/main.nf -profile test,docker
```

## After the pipeline

The scoring scripts read the pipeline outputs under `$CFDNA_DATA/phase5`:

| Script | Does |
|---|---|
| `scripts/phase5_transfer_train.py` | Trains the Phase 3b models once on the D = 30 simulation |
| `scripts/phase5_score.py` | Read-level mixing, all scores, the binomial gap ladder, external check |
| `scripts/phase5_summary.py` | Thresholds, detection rates, LoD, held-out FPR, AUC |
| `scripts/phase5_batch.py` | Batch control on off-marker features |
