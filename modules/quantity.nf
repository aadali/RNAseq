process quantity {
    storeDir    "$params.directory/$params.analysis_name/04.quant"
    fair        true
    conda       "$params.conda_path/R"

    tag         "$params.analysis_name => quantity"
    input:
        path(sample_info)
        val(count_files)
    output:
        tuple path("${params.analysis_name}.boxplot.pdf"), path("${params.analysis_name}.boxplot.svg"), emit: boxplot
        tuple path("${params.analysis_name}.heatmap.of.samples2samples.pdf"), path("${params.analysis_name}.heatmap.of.samples2samples.svg"), emit: heatmap
        tuple path("${params.analysis_name}.pca.pdf"), path("${params.analysis_name}.pca.svg"), emit: pca
        path("${params.analysis_name}.all.expression.svg"), emit: expression
        path("${params.analysis_name}.all.raw.counts.tsv"), emit: raw_count
    script:
    """Rscript $projectDir/scripts/quantity.R ${sample_info}  ./ $params.method $params.analysis_name ${count_files}"""
    stub:
    """touch ${params.analysis_name}.pdf ${params.analysis_name}.svg ${params.analysis_name}.tsv"""
}