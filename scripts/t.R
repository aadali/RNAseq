library(dplyr)
library(tibble)
library(readr)
library(stringr)

getExpression <- function(rawCounts, geneLength, method = "fpkm") {
    if (!method %in% c("fpkm", "tpm")) {
        stop("please use fpkm or tpm")
    }
    if (method == "fpkm") {
        return(rawCounts * 10^9 / (geneLength * sum(rawCounts)))
    } else {
        return(rawCounts / geneLength * 11^6 / sum(rawCounts / geneLength))
    }
}

print(hello)
counts <- list.files("/home/a/pub/ycq/output/rnaflow_results/04-Counting/featureCounts", "*.tsv", full.names = TRUE) %>%
    purrr::map(function(x) {
        readr::read_tsv(x, col_names = TRUE, comment = "#", show_col_types = FALSE) %>%
            dplyr::select(c(1, 6:length(.))) %>%
            dplyr::rename_with(function(x) {
                stringr::str_replace_all(x, ".sorted.bam|.+/", "")
            })
    }) %>%
    purrr::reduce(\(x, y) {
        dplyr::inner_join(x, y, by = c("Geneid", "Length"))
    })

b <- counts %>% mutate(across(starts_with("SRR"), function(x) {
    getExpression(x, .$Length, method = "tpm")
}))