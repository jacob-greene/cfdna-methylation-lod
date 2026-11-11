// Analysis: fragment patterns -> scores, thresholds, LoD and FPR -> donor bootstrap -> figures.
// Each process runs one script of scripts/ inside the repository. The scripts read the read
// outputs under params.data_dir ($CFDNA_DATA) and write their tables to results/ and figures/.
// The processes are not cached: a rerun recomputes every table from the current code.
// Scoring rules are fixed in pipeline/PREREGISTERED_PHASE6.md; the scoring scripts refuse to
// run if those rules or the scoring code are uncommitted or changed.
nextflow.enable.dsl = 2

// Run a script from the repository root with the project's Python environment.
def py(cmd) {
    """
    cd ${params.repo}
    export CFDNA_DATA=${params.data_dir}
    uv run --offline python ${cmd}
    """
}

// publishDir copies run asynchronously: wait up to 10 min for the read outputs the scripts read.
def wait_for(paths) {
    def tries = 0
    while (paths.any { !file(it).exists() } && tries++ < 60) { sleep(10000) }
    def missing = paths.findAll { !file(it).exists() }
    if (missing) { error("read outputs missing: ${missing.join(', ')}") }
    true
}

// Length-matched HCT116 spike (amendment 1 of the Phase 6 pre-registration).
process LENMATCH {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    """
    cd ${params.repo}
    REPO_DIR=${params.repo} CFDNA_DATA=${params.data_dir} bash scripts/phase6_lenmatch.sbatch
    """
}

// Classifiers trained on the array simulation, applied to real reads in the Phase 5 scores.
process TRANSFER_TRAIN {
    tag "$features"
    label 'analysis'
    cache false
    input:
    val features
    output:
    val true
    script:
    py("scripts/phase5_transfer_train.py --features ${features} --threads ${task.cpus}")
}

// Phase 5: mix the HCT116 spike into donor pools, score, thresholds on the 14 threshold donors.
process PHASE5_SCORE {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    py("scripts/phase5_score.py --workers ${task.cpus} && uv run --offline python scripts/phase5_summary.py")
}

// Phase 6: FRAG and WZ null model on the threshold donors, then the 16 test donors.
process PHASE6_SCORE {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    py("scripts/phase6_score.py --stage threshold && uv run --offline python scripts/phase6_score.py --stage eval"
        + " && uv run --offline python scripts/phase6_summary.py")
}

// Tumour fraction of each patient run from copy number: GC per bin, read counts per bin, TF.
process TF_GC {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    """
    mkdir -p ${params.gse_outdir}/tf
    ${py("scripts/phase6_tf.py gc ${params.ref_dir}/${params.fasta_name} ${params.gse_outdir}/tf/gc_1mb.tsv")}
    """
}

process TF_BINS {
    tag "$run"
    label 'analysis'
    cache false
    input:
    val run
    output:
    val true
    script:
    """
    mkdir -p ${params.gse_outdir}/tf/bins
    ${py("scripts/phase6_tf.py bins ${params.gse_outdir}/bam/${run}.md.bam ${params.gse_outdir}/tf/bins/${run}.bins.tsv.gz")}
    """
}

process TF_CALL {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    def g = params.gse_outdir
    py("scripts/phase6_tf.py tf ${params.gse_runs} ${g}/tf/bins ${g}/tf/gc_1mb.tsv ${g}/tf")
}

// Patient-plasma spike (amendment 2): score each patient run in the test pools, then summarise.
process PLASMA_SCORE {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    """
    cp ${params.gse_runs} ${params.gse_outdir}/runs.tsv
    ${py("scripts/phase6_gse_score.py && uv run --offline python scripts/phase6_gse_summary.py")}
    """
}

// Donor bootstrap: intervals for detection, FPR and LoD.
process BOOTSTRAP {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    py("scripts/bootstrap_lod.py")
}

process FIGURES {
    label 'analysis'
    cache false
    input:
    val ready
    output:
    val true
    script:
    py("scripts/key_findings_figure.py && uv run --offline python scripts/phase6_figures.py")
}

workflow ANALYSIS {
    take:
    cohort_qc   // QC table of the cohort read run: all cohort outputs are written
    plasma_qc   // QC table of the patient-plasma read run
    main:
    def p5 = "${params.data_dir}/phase5"
    def g = params.gse_outdir
    def cohort = file(params.samples).splitCsv(header: true, sep: '\t').collect { it.run } + [params.spike_run]
    def plasma = file(params.gse_runs).splitCsv(header: true, sep: '\t').collect { it.run }
    cohort_ready = cohort_qc.map { wait_for(cohort.collect { r -> "${p5}/pat/${r}.dedup.pat.gz" } + ["${p5}/blocks_cpg.tsv"]) }
    plasma_ready = plasma_qc.map { wait_for(plasma.collect { r -> "${g}/pat/${r}.dedup.pat.gz" } + ["${g}/blocks_cpg.tsv"]) }

    TRANSFER_TRAIN(channel.of('atlas', 'tumour'))
    PHASE5_SCORE(cohort_ready.combine(TRANSFER_TRAIN.out.collect()).first())
    LENMATCH(cohort_ready)
    PHASE6_SCORE(PHASE5_SCORE.out.combine(LENMATCH.out).first())

    TF_GC(plasma_ready)
    TF_BINS(plasma_ready.combine(channel.fromList(plasma)).map { it[1] })
    TF_CALL(TF_GC.out.combine(TF_BINS.out.collect()).first())
    PLASMA_SCORE(PHASE6_SCORE.out.combine(TF_CALL.out).first())

    BOOTSTRAP(PLASMA_SCORE.out)
    FIGURES(BOOTSTRAP.out)
}
