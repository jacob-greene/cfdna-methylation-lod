# Pipeline

This Nextflow pipeline runs the whole real-read analysis in one run:

1. `READS_COHORT` turns the SPOT-MAS cfDNA runs and the spike source (HCT116 WGBS) into
   fragment-level methylation patterns at the marker blocks (`<data>/phase5`).
2. `READS_PLASMA` does the same for the GSE186007 colorectal-cancer and healthy plasma runs
   (`samples/phase6_gse_runs.tsv`; `<data>/phase6_gse`).
3. `ANALYSIS` runs the scoring, LoD, tumour-fraction, bootstrap and figure scripts. They write
   `results/` and `figures/` in the repository.

Both read steps use the same processes (`modules/reads.nf`). The processing rules are fixed in
`PREREGISTERED.md`, the scoring rules in `PREREGISTERED_PHASE6.md`.

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

Analysis (`modules/analysis.nf`; each process runs one or two scripts of `scripts/`):

| Process | Script | Output |
|---|---|---|
| `TRANSFER_TRAIN` | `phase5_transfer_train.py` | Classifiers trained on the array simulation (`<data>/phase5/transfer`) |
| `PHASE5_SCORE` | `phase5_score.py`, `phase5_summary.py` | `results/phase5/` |
| `LENMATCH` | `phase6_lenmatch.sbatch` | Length-matched HCT116 spike (`<data>/phase5/pat`) |
| `PHASE6_SCORE` | `phase6_score.py` (threshold, eval), `phase6_summary.py` | `results/phase6/` |
| `TF_GC`, `TF_BINS`, `TF_CALL` | `phase6_tf.py` | Tumour fraction per patient run (`<data>/phase6_gse/tf`) |
| `PLASMA_SCORE` | `phase6_gse_score.py`, `phase6_gse_summary.py` | `results/phase6_gse/` |
| `BOOTSTRAP` | `bootstrap_lod.py` | `results/bootstrap/` |
| `FIGURES` | `key_findings_figure.py`, `phase6_figures.py` | `figures/` |

The analysis processes are not cached, so a rerun recomputes every table. The scoring scripts
refuse to run if the pre-registration or the scoring code are uncommitted or changed.

## Inputs

All inputs live in the data directory `$CFDNA_DATA` (default `../cfdna-methylation-lod-data`
next to the repository):

| Path in `$CFDNA_DATA` | Contents | Source |
|---|---|---|
| `sra/PRJNA929650_ena_filereport.tsv` | SPOT-MAS FASTQ URLs and MD5 sums | `data/download_sra_metadata.sh` |
| `sra/PRJNA771691_ena_filereport.tsv` | GSE186007 FASTQ URLs and MD5 sums | ENA file report of PRJNA771691 (same fields) |
| `ref/hg38/genome.fa` (+ `.fai`) | hg38 FASTA; the bwa-meth index is built on the first run | UCSC hg38 |
| `ref/annot/hg38.cytoBand.txt.gz` | Centromere positions for the tumour fraction | UCSC hg38 `cytoBand.txt.gz` |
| `tools/wgbs_tools-0.3.0` | wgbstools 0.3.0 with an initialised `hg38` genome | as in `containers/wgbstools/Dockerfile` |
| `tools/p5env` | Tool environment (`conda` profile) | `envs/phase5.explicit.txt` |
| `phase3b/` | Array simulation inputs for `TRANSFER_TRAIN` | `scripts/phase3b.sbatch` with `STEP=prepare` |

## Profiles

| Profile | Use |
|---|---|
| `slurm` | Slurm executor; partition from `NXF_SLURM_QUEUE`; at most 8 jobs at a time |
| `conda` | Pinned micromamba environment (`envs/phase5.explicit.txt`) with the container tool versions |
| `apptainer` | Biocontainers through Apptainer; the wgbstools image from `containers/wgbstools/` |
| `docker` | The same containers through Docker (used by the CI smoke test) |
| `awsbatch` | AWS Batch executor; written for portability, not run in this project |
| `test` | Smoke test on a small bundled BAM (`tests/data/`) |
| `test_fastq` | FASTQ-deletion test: FASTQ made from the bundled BAM, aligned, then deleted (run by `tests/fastq_cleanup_test.sh`) |

The HPC run used `slurm,conda`. Apptainer could not start containers in the restricted
environment of that cluster, also not in Slurm jobs (no setuid starter, no mount propagation). The
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
sbatch -p campus-new --export=ALL,REPO_DIR=$PWD,CFDNA_DATA=/path/to/data scripts/run_pipeline.sbatch
```

The script runs `nextflow run pipeline/main.nf -profile slurm,conda -resume` in
`$CFDNA_DATA/nxf`.

- `--gse false --analysis false` processes the cohort reads only.
- `--analysis false` stops after the reads.

- `--only_runs SRR1,SRR2` restricts the SPOT-MAS runs (tests). Add `--analysis false`: the
  analysis needs every run.
- `--spike_run false` (or `''`) skips the spike source.
- FASTQ deletion is the default. After a run's BAM passes `BAM_CHECK` (`samtools quickcheck`
  and more than 0 mapped reads in `samtools idxstats`), and the published copy in
  `<outdir>/bam` has the same size, the pipeline deletes the run's raw FASTQ (`fastq_dir`) and
  its trimmed FASTQ (the `FASTP` work directory). A failed check keeps the files and logs a
  warning. Every decision goes to `<outdir>/fastq_cleanup.tsv`. `--keep_fastq` keeps all FASTQ.
- `--reuse_bams` (default true): a run whose BAM and index are already in `<outdir>/bam` starts
  from that BAM. A rerun therefore does not fetch and align the runs whose FASTQ are gone.
  Their downstream tasks rerun once, because the BAM path differs from the cached one.
- `--input_fastq 'dir/*_{1,2}.fastq.gz'` runs local paired FASTQ instead of ENA downloads.
  These input files may be the only copy, so the pipeline keeps them and deletes only the
  trimmed FASTQ. `--delete_input_fastq` deletes them too.
- Run the Nextflow head process inside a Slurm job; it runs for hours.
  `scripts/run_pipeline.sbatch` does this. Pass `REPO_DIR` and `CFDNA_DATA` with
  `sbatch --export=ALL,REPO_DIR=...,CFDNA_DATA=...`.

Smoke test (Docker, as in CI):

```bash
docker build -t cfdna-wgbstools:0.3.0 pipeline/containers/wgbstools
nextflow run pipeline/main.nf -profile test,docker
pipeline/tests/fastq_cleanup_test.sh docker   # needs samtools on PATH (or SAMTOOLS=...)
```
