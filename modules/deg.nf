process get_comparisons {
    tag     "get_comparisions"
    input:
        path(comparison_info)
    output:

}

process deg {
    storeDir    "$params.directory/$params.analysis_name/05.deg/$comparison"
    fair        true
    conda       "$params.analysis_name => deg"
    tag         "$params.analysis_name => $comparison"

    input:
        val(comparison)
        path(raw_counts)
}