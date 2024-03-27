process fastp {
    storeDir "$params.directory/$params.analysis_name/01.clean_data"
    fair true
    tag "${meta.sample}"

    input:
        val(meta)

    output:
        val(meta), emit: meta
        path("${meta.sample}.clean_?.fastq.gz"), emit: clean_reads
        path("${meta.sample}.log")
        path("${meta.sample}.json")
        path("${meta.sample}.html")

    script:
    """fastp \
    -i ${meta.fq1} \
    -I ${meta.fq2} \
    -o ${meta.sample}.clean_1.fastq.gz \
    -O ${meta.sample}.clean_2.fastq.gz \
    --thread $params.cpus \
    --json ${meta.sample}.json \
    --html ${meta.sample}.html \
    $params.extra_fastp_args 2> ${meta.sample}.log """

    stub:
    """touch ${meta.sample}.clean_1.fastq.gz && \
    touch ${meta.sample}.clean_2.fastq.gz && \
    touch ${meta.sample}.json && \
    touch ${meta.sample}.html && \
    touch ${meta.sample}.log """
}

