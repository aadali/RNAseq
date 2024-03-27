include {fastqc as fastqcPre; fastqc  as fastqcPost} from './modules/fastqc'
include {fastp} from './modules/fastp'
include {index_reference} from './modules/hisat2'
include {hisat2} from './modules/hisat2'
include {featureCounts} from './modules/featureCounts'
include {quantity} from './modules/quantity.nf'


workflow fromFastqToCounts {
    // sample_info = Channel.fromPath("$params.sample_info", checkIfExists: true)
    take:
        sample_info
    main:
        sample_info.splitCsv(header: true).map {
            row -> {
                if (!row.sample.startsWith("#")) {
                    if (!row.fq1) { log.error("fq1 must be set"); eixt(2) }
                    if (row.fq2 == '') { row['paired'] =  false } else { row['paired'] = true }
                    return row
                } 
            }
        }.set{fastpIn}
        fastp(fastpIn).set{fastpOut} // fastpOut: meta, [clean_reads], log, json, html
        index_reference(params).set{reference}
        hisat2(fastpOut.meta, fastpOut.clean_reads, reference).set{hisat2Out} // hisat2Out: meta, [bam, bai], log
        featureCounts(hisat2Out.meta, hisat2Out.bam).set{featureCountsOut} // featureCountsOut: meta, count.tsv, summary, log
    emit:
        bam_bai = hisat2Out.bam
        count_tsv = featureCountsOut.count
}

workflow {
    sample_info = Channel.fromPath("$params.sample_info", checkIfExists: true)
    a = fromFastqToCounts(sample_info)
    a.count_tsv.collect().map{it -> it.join(",")}.set{count_files}

    // count_files.view()
    quantity(sample_info, count_files).set{quantityOut}
    quantityOut[0].view()
    quantityOut[1].view()
    quantityOut[2].view()

}
