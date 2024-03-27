// sample_info.csv strand
// 0: nonstrand
// 1: first
// 2: second

process featureCounts {
    storeDir    "$params.directory/$params.analysis_name/03.counts"
    tag         "$meta.sample"
    fair        true

    input:
        val(meta)
        tuple path(bam), path(bai)
    output:
        val(meta), emit: meta
        path("${meta.sample}.counts.tsv"), emit: count
        path("${meta.sample}.counts.tsv.summary"), emit: summary
        path("${meta.sample}.log"), emit: log
    
    script:
    /*
    Count fragments which have both ends successfully aligned without considering the fragment length constraint
    -p  If specified, libraries are assumed to contain paired-end reads
    --countReadPairs  If specified, fragments (or templates) will be counted instead of reads. For paired-end reads, if not specified, count may be twice
    -B  Only count read pairs that have both ends aligned
    */
    def refParentPath = file(params.reference).getParent()
    def paired_para = meta.paired ? "-p --countReadPairs -B" : ""
    """featureCounts \
    ${paired_para} \
    -a $refParentPath/genomic.gtf \
    $params.feature_attr \
    -o ${meta.sample}.counts.tsv \
    $bam 2> ${meta.sample}.log
    """

    stub:
    """touch ${meta.sample}.counts.tsv && \
    touch ${meta.sample}.counts.tsv.summary && \
    touch ${meta.sample}.log"""

}