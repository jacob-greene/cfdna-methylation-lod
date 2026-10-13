#!/usr/bin/env nextflow
// Phase 5 read-level pipeline: FASTQ -> trimmed reads -> bwa-meth (hg38) -> duplicate-marked BAM
// -> MethylDackel per-CpG counts + M-bias -> wgbstools fragment patterns -> per-sample QC table.
// A run's raw and trimmed FASTQ are deleted once its BAM passes BAM_CHECK (--keep_fastq keeps them).
// Every source (SPOT-MAS runs and the HCT116 spike source) goes through the same steps.
// Processing rules are fixed in pipeline/PREREGISTERED.md.
nextflow.enable.dsl = 2

// Run table: SPOT-MAS runs from the pre-registered sample sheet plus the spike source.
def ena_rows() {
    def rows = [:]
    file(params.ena_report).splitCsv(header: true, sep: '\t').each { r -> rows[r.run_accession] = r }
    rows
}

// wgbstools finds its reference genome relative to its own source tree, so every wgbstools task
// works on a private copy of the tree with the genome directory linked in.
def wt_setup(ref) {
    """
    mkdir -p wt/references && cp -r ${params.wgbs_home}/src wt/
    ${ref ? "ln -s \$(readlink -f ${ref}) wt/references/${params.wgbs_genome}" : ''}
    # GenomicRegion() reads references/default before it reads --genome (wgbstools 0.3.0).
    ln -sfn ${params.wgbs_genome} wt/references/default
    wt="${params.wgbs_python} wt/src/python/wgbs_tools.py"
    """
}

// Delete a run's raw and trimmed FASTQ after BAM_CHECK passed and the published BAM copy is
// complete (same size as the checked BAM). Anything else keeps the files. Every decision is
// appended to <outdir>/fastq_cleanup.tsv.
def clean_fastq(run, fastqs, status) {
    def log_tsv = file("${params.outdir}/fastq_cleanup.tsv")
    def f = status.trim().tokenize(' ')
    def pub = file("${params.outdir}/bam/${run}.md.bam")
    def ok = f[0] == 'pass'
    if (ok) {
        // publishDir copies run asynchronously; wait up to 10 min for the copy to finish.
        def tries = 0
        while (!(pub.exists() && pub.size() == (f[2] as long)) && tries++ < 60) { sleep(10000) }
        ok = pub.exists() && pub.size() == (f[2] as long)
        if (!ok) { f = ['fail', 'published_bam_incomplete'] }
    }
    fastqs.each { q ->
        if (!q.exists()) { return }   // deleted on an earlier run (-resume replays BAM_CHECK)
        def action = ok ? (q.delete() ? 'deleted' : 'delete_failed') : 'kept'
        if (!ok) { log.warn "FASTQ kept for ${run} (${f[1]}): ${q}" }
        log_tsv.append("${run}\t${q}\t${pub}\t${f.join(' ')}\t${action}\n")
    }
}

process FETCH_ENA {
    tag "$run"
    label 'fetch'
    storeDir params.fastq_dir
    input:
    tuple val(run), val(urls), val(md5s)
    output:
    tuple val(run), path("${run}_1.fastq.gz"), path("${run}_2.fastq.gz")
    script:
    """
    curl -fsSL --retry 5 --retry-delay 10 -o ${run}_1.fastq.gz https://${urls[0]}
    curl -fsSL --retry 5 --retry-delay 10 -o ${run}_2.fastq.gz https://${urls[1]}
    printf '%s  %s\\n%s  %s\\n' ${md5s[0]} ${run}_1.fastq.gz ${md5s[1]} ${run}_2.fastq.gz | md5sum -c -
    """
}

// Spike source and --stream_runs: the first N read pairs of each mate, streamed (a full run
// can be much larger than N pairs).
process FETCH_STREAM {
    tag "$run"
    label 'fetch'
    // A broken stream fails the line-count check below; fetch again.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2
    storeDir params.fastq_dir
    input:
    tuple val(run), val(urls)
    output:
    tuple val(run), path("${run}_head${params.spike_pairs}_1.fastq.gz"), path("${run}_head${params.spike_pairs}_2.fastq.gz")
    script:
    def n = (params.spike_pairs as long) * 4
    """
    set +o pipefail   # head closes the pipe early by design
    for m in 1 2; do
      url=\$([ \$m = 1 ] && echo ${urls[0]} || echo ${urls[1]})
      curl -fsSL --retry 5 https://\$url | gzip -dc | head -n ${n} | gzip -1 > part_\$m.fastq.gz
      got=\$(gzip -dc part_\$m.fastq.gz | wc -l)
      [ "\$got" -eq ${n} ] || { echo "mate \$m: \$got lines, expected ${n}" >&2; exit 1; }
      mv part_\$m.fastq.gz ${run}_head${params.spike_pairs}_\$m.fastq.gz
    done
    """
}

process FASTP {
    tag "$run"
    label 'fastp'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.json'
    input:
    tuple val(run), path(r1), path(r2)
    output:
    tuple val(run), path("${run}_1.trim.fq.gz"), path("${run}_2.trim.fq.gz"), emit: reads
    path "${run}.fastp.json", emit: json
    script:
    """
    fastp -i $r1 -I $r2 -o ${run}_1.trim.fq.gz -O ${run}_2.trim.fq.gz \\
      --detect_adapter_for_pe --disable_quality_filtering --length_required 30 \\
      --trim_front2 ${params.clip_r2_5} --trim_tail1 ${params.clip_r1_3} --trim_tail2 ${params.clip_r2_3} \\
      --max_len1 ${params.max_len} --max_len2 ${params.max_len} \\
      -w ${task.cpus} -j ${run}.fastp.json -h ${run}.fastp.html
    """
}

process BWAMETH_INDEX {
    label 'bwameth'
    storeDir "${params.ref_dir}"
    input:
    path fasta
    output:
    path "${fasta.name}.bwameth.c2t*"
    script:
    """
    bwameth.py index ${fasta}
    """
}

// Align, fix mate tags, sort and mark duplicates (marked, not removed).
process BWAMETH_ALIGN {
    tag "$run"
    label 'bwameth'
    // 'link' (a hard link) avoids a second BAM copy when the work dir is deleted after the run.
    publishDir "${params.outdir}/bam", mode: params.bam_publish_mode, pattern: '*.{bam,bai,txt}'
    input:
    tuple val(run), path(r1), path(r2)
    path ref_dir
    output:
    tuple val(run), path("${run}.md.bam"), path("${run}.md.bam.bai"), emit: bam
    path "${run}.{markdup,stats}.txt", emit: qc
    script:
    def fa = "${ref_dir}/${params.fasta_name}"
    """
    bwameth.py --threads ${task.cpus} --reference ${fa} --read-group '@RG\\tID:${run}\\tSM:${run}' $r1 $r2 \\
      | samtools fixmate -m -@ 2 - - \\
      | samtools sort -@ 4 -m 1500M -T tmp_${run} - \\
      | samtools markdup -@ 4 -s -f ${run}.markdup.txt - ${run}.md.bam
    samtools index -@ 4 ${run}.md.bam
    samtools stats -@ 4 ${run}.md.bam > ${run}.stats.txt
    """
}

// Gate for FASTQ deletion: the BAM must pass samtools quickcheck and have mapped reads.
// Prints "pass <mapped reads> <BAM bytes>" or "fail <reason>"; a failed check keeps the FASTQ.
process BAM_CHECK {
    tag "$run"
    label 'samtools'
    input:
    tuple val(run), path(bam), path(bai)
    output:
    tuple val(run), stdout
    script:
    """
    if ! samtools quickcheck ${bam}; then echo "fail quickcheck"; exit 0; fi
    mapped=\$(samtools idxstats ${bam} | awk '{s+=\$3} END{print s+0}')
    if [ "\$mapped" -gt 0 ]; then echo "pass \$mapped \$(stat -L -c %s ${bam})"; else echo "fail no_mapped_reads"; fi
    """
}

// Per-CpG counts (duplicates ignored, MethylDackel default), M-bias, and CHH calls on one
// chromosome for the bisulfite conversion rate.
process METHYLDACKEL {
    tag "$run"
    label 'methyldackel'
    publishDir "${params.outdir}/methyldackel", mode: 'copy'
    input:
    tuple val(run), path(bam), path(bai)
    path ref_dir
    output:
    tuple val(run), path("${run}_CpG.bedGraph.gz"), path("${run}.mbias.txt"), path("${run}.chh.txt")
    script:
    def fa = "${ref_dir}/${params.fasta_name}"
    """
    MethylDackel extract -@ ${task.cpus} -o ${run} ${fa} ${bam}
    gzip -f ${run}_CpG.bedGraph
    MethylDackel mbias -@ ${task.cpus} --txt ${fa} ${bam} ${run} > ${run}.mbias.txt
    MethylDackel extract -@ ${task.cpus} --CHH --noCpG -r ${params.conversion_chrom} -o conv ${fa} ${bam}
    awk 'NR>1{m+=\$5; u+=\$6} END{printf "chh_meth\\t%d\\nchh_unmeth\\t%d\\n", m, u}' conv_CHH.bedGraph > ${run}.chh.txt
    """
}

// Fragment patterns in marker and off-marker blocks: duplicates removed (wgbstools default
// flags 1796) and duplicates kept (flags 772), for rung (4) of the binomial gap ladder.
process BAM2PAT {
    tag "$run"
    label 'wgbstools'
    // bam2pat's temp-dir cleanup can hit "Directory not empty" on NFS (.nfs files); retry.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2
    publishDir "${params.outdir}/pat", mode: 'copy'
    input:
    tuple val(run), path(bam), path(bai)
    path blocks
    path offmarker
    path wgbs_ref
    output:
    tuple val(run), path("${run}.dedup.pat.gz"), path("${run}.withdup.pat.gz")
    script:
    """
    mkdir -p dd wd   # bam2pat requires an existing output directory
    cat ${blocks} ${offmarker} | cut -f1-3 | sort -k1,1 -k2,2n > regions.bed
    ${wt_setup(wgbs_ref)}
    \$wt bam2pat --genome ${params.wgbs_genome} -L regions.bed -@ ${task.cpus} -F 1796 -o dd ${bam}
    \$wt bam2pat --genome ${params.wgbs_genome} -L regions.bed -@ ${task.cpus} -F 772 -o wd ${bam}
    mv dd/${bam.baseName}.pat.gz ${run}.dedup.pat.gz
    mv wd/${bam.baseName}.pat.gz ${run}.withdup.pat.gz
    """
}

// Technical features of read-1 fragments in off-marker blocks, for the batch control:
// insert size, mapping quality, mismatches (NM) and read length. Duplicates excluded.
process FRAG_FEATURES {
    tag "$run"
    label 'samtools'
    publishDir "${params.outdir}/frag_features", mode: 'copy'
    input:
    tuple val(run), path(bam), path(bai)
    path offmarker
    output:
    tuple val(run), path("${run}.frag_features.tsv.gz")
    script:
    """
    samtools view -f 67 -F 1796 -q 10 -L ${offmarker} ${bam} \\
      | awk -v OFS='\\t' 'BEGIN{print "chrom","pos","tlen","mapq","nm","read_len"}
          {nm="NA"; for(i=12;i<=NF;i++) if(substr(\$i,1,5)=="NM:i:") nm=substr(\$i,6);
           t=\$9<0?-\$9:\$9; print \$3,\$4,t,\$5,nm,length(\$10)}' \\
      | gzip > ${run}.frag_features.tsv.gz
    """
}

// Initialise a wgbstools genome from the FASTA (only when params.wgbs_ref is unset, e.g. tests).
process WGBS_INIT {
    label 'wgbstools'
    input:
    path ref_dir
    output:
    path "wgbs_ref"
    script:
    """
    ${wt_setup(null)}
    \$wt init_genome ${params.wgbs_genome} --fasta_path \$(readlink -f ${ref_dir}/${params.fasta_name}) --no_default -@ ${task.cpus}
    cp -rL wt/references/${params.wgbs_genome} wgbs_ref
    """
}

// Marker and off-marker blocks with their CpG index ranges (wgbstools convert), for QC and
// scoring.
process BLOCK_INDEX {
    label 'wgbstools'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path blocks
    path offmarker
    path wgbs_ref
    output:
    path "blocks_cpg.tsv", emit: blocks
    path "offmarker_cpg.tsv", emit: offmarker
    script:
    """
    ${wt_setup(wgbs_ref)}
    \$wt convert -L ${blocks} --genome ${params.wgbs_genome} --no_anno -p > blocks_cpg.tsv
    \$wt convert -L ${offmarker} --genome ${params.wgbs_genome} --no_anno -p > offmarker_cpg.tsv
    """
}

process QC_TABLE {
    label 'wgbstools'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path files
    path blocks_cpg
    path qc_script  // staged, so the script is visible inside a container
    output:
    path "qc_table.tsv"
    script:
    """
    ${params.wgbs_python} ${qc_script} --blocks ${blocks_cpg} --out qc_table.tsv .
    """
}

workflow {
    ref_dir = file(params.ref_dir, checkIfExists: true)
    bams = null
    reused = null
    reused_qc = null
    if (params.input_bams) {
        // Smoke-test entry: start from aligned, duplicate-marked BAMs.
        bams = channel.fromPath(params.input_bams)
            .map { b -> tuple(b.simpleName, b, file("${b}.bai")) }
    } else {
        if (params.input_fastq) {
            // Local paired FASTQ ('dir/*_{1,2}.fastq.gz'); the test_fastq profile uses this.
            reads = channel.fromFilePairs(params.input_fastq, checkIfExists: true)
                .map { run, fq -> tuple(run, fq[0], fq[1]) }
        } else {
            ena = ena_rows()
            // --stream_runs: runs from the ENA report, the first spike_pairs read pairs each,
            // in place of the sample sheet.
            stream = params.stream_runs ? params.stream_runs.tokenize(',') : []
            runs = stream ?: file(params.samples).splitCsv(header: true, sep: '\t').collect { it.run }
            if (params.only_runs) { runs = runs.findAll { it in params.only_runs.tokenize(',') } }
            // `--spike_run ''` and `--spike_run false` reach here as booleans: skip the spike.
            def spike_run = (params.spike_run instanceof CharSequence && params.spike_run) ? params.spike_run : null
            // A run whose BAM is already published (its FASTQ may be deleted) starts from that BAM.
            def published = { r -> file("${params.outdir}/bam/${r}.md.bam").exists() && file("${params.outdir}/bam/${r}.md.bam.bai").exists() }
            done = params.reuse_bams ? (runs + (spike_run ? [spike_run] : [])).findAll(published) : []
            if (done) { log.info "Starting ${done.size()} runs from their published BAMs: ${done.join(',')}" }
            todo = (runs - done).collect { r -> tuple(r, ena[r] ?: error("run $r missing from ${params.ena_report}")) }
            streamed = stream ? todo.collect { r, e -> tuple(r, e.fastq_ftp.tokenize(';')) } : []
            if (spike_run && !(spike_run in done)) { streamed << tuple(spike_run, params.spike_urls.tokenize(';')) }
            spotmas = channel.fromList(stream ? [] : todo.collect { r, e ->
                tuple(r, e.fastq_ftp.tokenize(';'), e.fastq_md5.tokenize(';'))
            })
            reads = FETCH_ENA(spotmas)
            if (streamed) { reads = reads.mix(FETCH_STREAM(channel.fromList(streamed))) }
            reused = channel.fromList(done.collect { r ->
                tuple(r, file("${params.outdir}/bam/${r}.md.bam"), file("${params.outdir}/bam/${r}.md.bam.bai"))
            })
            reused_qc = channel.fromList(done.collect { r ->
                [file("${params.outdir}/fastp/${r}.fastp.json"), file("${params.outdir}/bam/${r}.markdup.txt"),
                 file("${params.outdir}/bam/${r}.stats.txt")]
            }).flatten().filter { it.exists() }
        }
        c2t = file("${params.ref_dir}/${params.fasta_name}.bwameth.c2t.sa")
        index_ready = c2t.exists() ? channel.value(ref_dir)
            : BWAMETH_INDEX(file("${params.ref_dir}/${params.fasta_name}")).map { ref_dir }
        FASTP(reads)
        BWAMETH_ALIGN(FASTP.out.reads, index_ready)
        bams = BWAMETH_ALIGN.out.bam
        if (!params.keep_fastq) {
            BAM_CHECK(BWAMETH_ALIGN.out.bam)
            // Local --input_fastq files may be the only copy: they are deleted only with
            // --delete_input_fastq. Downloaded FASTQ and trimmed FASTQ are always deleted.
            def raw_too = !params.input_fastq || params.delete_input_fastq
            reads.join(FASTP.out.reads).join(BAM_CHECK.out)
                .subscribe { run, r1, r2, t1, t2, status ->
                    clean_fastq(run, (raw_too ? [r1, r2] : []) + [t1, t2], status)
                }
        }
        if (reused != null) { bams = bams.mix(reused) }
    }
    blocks = file(params.blocks, checkIfExists: true)
    offmarker = file(params.offmarker, checkIfExists: true)
    wgbs_ref = params.wgbs_ref ? channel.value(file(params.wgbs_ref, checkIfExists: true))
        : WGBS_INIT(ref_dir)
    METHYLDACKEL(bams, ref_dir)
    BAM2PAT(bams, blocks, offmarker, wgbs_ref)
    BLOCK_INDEX(blocks, offmarker, wgbs_ref)
    FRAG_FEATURES(bams, offmarker)
    qc_in = METHYLDACKEL.out.map { it[1..3] }.flatten()
        .mix(BAM2PAT.out.map { it[1] })
    if (!params.input_bams) {
        qc_in = qc_in.mix(FASTP.out.json, BWAMETH_ALIGN.out.qc.flatten())
        if (reused_qc != null) { qc_in = qc_in.mix(reused_qc) }
    }
    QC_TABLE(qc_in.collect(), BLOCK_INDEX.out.blocks, file("${projectDir}/bin/phase5_qc.py"))
}
