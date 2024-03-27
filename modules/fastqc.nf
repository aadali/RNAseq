process fastqc {
    storeDir    "$params.directory/$params.analysis_name/00.qc"
    fair        true
    input:
        val(meta)
        path(reads)
    
    output:
        path("fastqc*.zip"), emit: zip
    script:
    """fastqc   --noextract -t $params.cpus ${reads}"""
}