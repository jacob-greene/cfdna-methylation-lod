#!/usr/bin/env nextflow
// cfDNA methylation LoD pipeline: public reads -> fragment patterns -> scores -> LoD tables and
// figures. One run does every step:
//   1. READS_COHORT: SPOT-MAS healthy-donor runs (samples/phase5_samples.tsv) and the HCT116
//      spike source -> <data>/phase5.
//   2. READS_PLASMA: GSE186007 colorectal-cancer and healthy plasma runs
//      (samples/phase6_gse_runs.tsv), the first spike_pairs read pairs of each -> <data>/phase6_gse.
//   3. ANALYSIS: scores, thresholds, LoD, tumour fraction, bootstrap -> results/ and figures/.
// Processing rules are fixed in PREREGISTERED.md and PREREGISTERED_PHASE6.md.
nextflow.enable.dsl = 2

def plasma_runs = file(params.gse_runs).splitCsv(header: true, sep: '\t').collect { it.run }.join(',')

include { READS as READS_COHORT } from './modules/reads.nf'
include { READS as READS_PLASMA } from './modules/reads.nf' addParams(
    outdir: params.gse_outdir, fastq_dir: "${params.gse_outdir}/fastq", ena_report: params.gse_ena_report,
    stream_runs: plasma_runs, spike_run: false, only_runs: null, bam_publish_mode: 'link')
include { ANALYSIS } from './modules/analysis.nf'

workflow {
    READS_COHORT()
    if (params.gse) { READS_PLASMA() }
    if (params.analysis) {
        if (!params.gse) { error('--analysis needs the patient-plasma reads: drop --gse false or set --analysis false') }
        ANALYSIS(READS_COHORT.out.qc, READS_PLASMA.out.qc)
    }
}
