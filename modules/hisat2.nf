process index_reference {
    input:
        val params
    output:
        val "$params.reference"
    script:
    """if [ ! -e ${params.reference}.1.ht2 ]; then hisat2-build -p 24 ${params.reference} ${params.reference}; fi"""
}

// sample_info.csv strand
// 0: nonstrand
// 1: first
// 2: second

process hisat2 {
    storeDir    "$params.directory/$params.analysis_name/02.mapping"
    maxForks    Math.floor("nproc".execute().text.toInteger()/params.cpus).toInteger()
    fair        true
    tag         "${meta.sample}"

    input:
        val(meta)
        path(reads)
        val(reference) // treat reference from index_reference as a value, actually， it's a filepath
    output:
        val(meta), emit: meta
        tuple path("${meta.sample}.sorted.bam"), path("${meta.sample}.sorted.bam.bai"), emit: bam
        path("${meta.sample}_summary.log"), emit: log
    script:
    def read1_para = meta.paired ? " -1 ${reads[0]} " : " -U ${reads[0]} "
    def read2_para = meta.paired ? " -2 ${reads[1]} " : " "
    def strandNum = meta.strand.toInteger()

    // for single end reads, strandedness is not supported for pipeline
    if (!meta.paired && strandNum != 0) { 
        log.error("Single end standed library is not supported")
        exit(2)
    }
    def strand = ""
    if (strandNum == 1) { strand = "--rna-strandness RF " }
    if (strandNum == 2) { strand = "--rna-strandness FR " }

    """hisat2 \
    -x $reference \
    ${read1_para} \
    ${read2_para} \
    --threads $params.cpus \
    --new-summary \
    ${strand} \
    --summary-file ${meta.sample}_summary.log \
    $params.extra_hisat2_args | \
    samtools view -@ 4 -bS > ${meta.sample}.raw.bam && \
    samtools sort -@ $params.cpus -o ${meta.sample}.sorted.bam ${meta.sample}.raw.bam && \
    samtools index -@ 2 ${meta.sample}.sorted.bam && \
    rm ${meta.sample}.raw.bam"""

    stub:
    """touch ${meta.sample}.sorted.bam && \
    touch ${meta.sample}.sorted.bam.bai && 
    touch ${meta.sample}_summary.log"""
}