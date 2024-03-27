library(optparse)
# library(stringr)
library(dplyr)
# library(DESeq2)

get.args <- function() {
    options.list <- list(
        make_option(
            c('-o', '--output'),
            type = "character",
            default = ".",
            help = "the output directory"
        ),
        make_option(opt_str = c('-s', '--sample'),
                    type = "character",
                    help = "the sample info csv"),
        make_option(opt_str = c('-c', '--comparison'),
                    type = "character",
                    help = "the comparison groups info csv"),
        make_option('--count',
                    type = "character",
                    help = "the count file from featureCounts, multi files must be separated by comma"),
        make_option('--expression_threshold',
                    type = "double",
                    help = "the min expression threshold of fpkm or tpm"),
        make_option('--method',
                    type = "character",
                    default = "fpkm",
                    help = "which method to quantify the expression, must be one of (fpkm, tpm, normalization)"
        )
    )
    parser <- OptionParser(
        option_list = options.list,
        add_help_option = TRUE,
        usage = "Rscript %prog -s sample_info.tsv -c counts1.tsv[,counts2.tsv,counts2.tsv...]",
    )
    args <- parse_args(parser)

    if.file.exists <- function(x) {
        paths <- c(file.path(getwd(), x), x)
        a <- any(file.exists(paths))
        return(a)
    }

    if (is.null(args$sample)) {
        stop("--sample must be set")
    }
    if (!if.file.exists(args$sample)) {
        stop(stringr::str_glue("No such file for --sample {args$sample}"))
    }
    if (!is.null(args$comparison) &&
        !if.file.exists(args$comparison)) {
        stop(stringr::str_glue("No such file for --comparison {args$comparison}"))
    }
    if (!args$method %in% c("fpkm", "tpm", "normalization")) {
        stop(stringr::str_glue("method must be one of (fpkm, tpm, normalization), not {args$method}"))
    }
    return(args)
}

get.comparisons <- function(sample, comparison = NULL) {
    if (!is.null(comparison)) {
        comp <- readr::read_csv(comparison, comment = "#", show_col_types = FALSE, col_names = TRUE) %>%
            select(1, 2)
        colnames(comp) <- c("control", "treated") # the first column shoud be control and the second is treated
        a <- purrr::map2(comp$control, comp$treated, function(x, y) {c("condition", x, y)})
        return(a)
    }
    sample.group <- readr::read_csv(sample, comment = "#", show_col_types = FALSE)
    #todo: the colnames of sample shoud be: c("sample", "condition", "source", "fq1", "fq2", "strand")
    conditions <- unique(sample.group$condition)
    a <- list()
    while (1) {
        purrr::map2(conditions[1], conditions[2:length(conditions)],
                    function(x, y) {a[[length(a) + 1]] <<- c("condition", x, y)})
        conditions <- conditions[2:length(conditions)]
        if (length(conditions) == 1) {break}
    }
    return(a)
}

get.counts <- function(count.files) {
                        #' @param cont.files: a  character vector of outfiles from featureCounts
                        #' @return a list, 1st element is data.frame counts whose rownames and colnames are sample and gene names,
                        #' and 2nd element is the length of genes sorted by the counts rownames
    counts <- purrr::map(count.files, \(x) {
        readr::read_tsv(
            x,
            col_names = TRUE,
            comment = "#",
            show_col_types = FALSE
        ) %>%
            dplyr::select(c(1, 6:length(.))) %>%
            dplyr::rename_with(\(x) {
                # the format of bam's name shoud be ".+/{sampleName}.sorted.bam"
                # sampleName means the name of sample in the sample_info.csv and used to distinguish different samples"
                stringr::str_replace_all(x, ".sorted.bam|.+/", "")
            })
    }) %>%
        purrr::reduce(\(x, y) {dplyr::inner_join(x, y, by = c("Geneid", "Length"))})
    readr::write_tsv(counts, "all.raw.counts.tsv") # 保存原始的counts数到本地
    keep <- rowSums((counts > 5)) / length(counts) > 0.25 # 对于某个基因，如果feature counts数大于5的样本数量小于总样本数的1/4，则舍弃该基因
    # print(keep)
    counts <- counts[keep,] %>%
        tibble::column_to_rownames("Geneid")
    return(list(counts = counts[, names(counts) != 'Length'],
                length = counts$Length))
}

get.expression <- function(raw.counts, genes.length, method = "fpkm") {
    if (!method %in% c("fpkm", "tpm")) {
        stop("method must be one of ('fpkm', 'tpm')")
    }
    if (method == "fpkm") {
        return(raw.counts * 10^9 / (genes.length * sum(raw.counts)))
    } else {
        return(raw.counts / genes.length * 11^6 / sum(raw.counts / genes.length))
    }
}

save.expression <- function(counts, lengths, outdir, name, method = "fpkm") {
        #' @param counts: the counts data.frame from get.counts
        #' @param lengths: the gene lengths from from get.counts
        #' @param outdir: the output directory
        #' @param name: the output file's name
        #' @return the path of output file
        #' save fpkm or tpm results to disk
    expression <- dplyr::mutate_all(counts, function(x) {get.expression(x, lengths, method = method)}) %>%
        dplyr::rename_with(function(x) {paste(x, method, sep = "__")})
    outfile <- file.path(outdir, name)
    readr::write_tsv(expression, outfile)
    return(outfile)
}

expression.boxplot <- function(sample, expression.file, outdir, name) {
    groups <- readr::read_tsv(sample, col_names = TRUE, show_col_types = FALSE, col_select = c("sample", "condition"))
    readr::read_tsv(expression.file, col_names =)
}

deg <- function(counts, sample.info, comparison = NULL) {
    # col_names = c("sample", "condition", "fq1", "fq2", "strand")
    col.data <- readr::read_csv(sample.info, col_names = TRUE, comment = "#", show_col_types = FALSE) %>%
        dplyr::select(c("sample", "condition")) %>%
        dplyr::mutate(condition = factor(.$condition))
    dds <- DESeq2::DESeqDataSetFromMatrix(counts[, col.data$sample],
                                          colData = col.data,
                                          design = ~condition)
    return(dds)
}


main <- function() {
    message("\n\n\n=====================================START====================================")
    # args <- get.args()
    files <- c(
        "./test_data/SRR23185523.counts.tsv",
        "./test_data/SRR23185524.counts.tsv",
        "./test_data/SRR23185525.counts.tsv",
        "./test_data/SRR23185544.counts.tsv",
        "./test_data/SRR23185545.counts.tsv",
        "./test_data/SRR23185546.counts.tsv",
        "./test_data/SRR23185547.counts.tsv",
        "./test_data/SRR23185548.counts.tsv",
        "./test_data/SRR23185549.counts.tsv"
    )
    args <- list(
        count = paste(files, collapse = ","),
        sample = "sample_info.csv"
    )
    count.files <- stringr::str_split(args$count, ",")[[1]]
    raw.counts <- get.counts(count.files)
    counts <- raw.counts$counts
    length <- raw.counts$length
    dds <- deg(counts = counts,
               sample.info = args$sample,
               comparison = NULL)
    return(dds)
}

files <- c(
    "./test_data/SRR23185523.counts.tsv",
    "./test_data/SRR23185524.counts.tsv",
    "./test_data/SRR23185525.counts.tsv",
    "./test_data/SRR23185544.counts.tsv",
    "./test_data/SRR23185545.counts.tsv",
    "./test_data/SRR23185546.counts.tsv",
    "./test_data/SRR23185547.counts.tsv",
    "./test_data/SRR23185548.counts.tsv",
    "./test_data/SRR23185549.counts.tsv"
)
counts <- get.counts(files)
# get.counts(counts)
save.expression(counts$count, counts$length, outdir = "./test_data", name = "all.fpkm.tsv")

# dds <- main()
#
# dds <- DESeq(dds)
# # dds <- results(dds)
# dds.res <- results(dds, contrast = c("condition", "WT", "CrxE80A_A"), alpha = 0.05, lfcThreshold = 1)
# resultsNames(dds)
# summary(dds.res)
# dds.res[c('Tiam2', 'Pde6h'),]
# df <- read_tsv("./test_data/merged.counts.tsv", col_names = TRUE) %>%
#     column_to_rownames("geneid")
# col.data <- read_csv('sample_info.csv', comment = "#") %>%
#     select(sample, condition) %>%
#     mutate(condition = factor(.$condition))
#
# df <- df[, col.data$sample]
# keep <- (df > 5) %>% apply(1, sum) / ncol(df) > 0.25
# dds <- DESeqDataSetFromMatrix(df[keep,], colData = col.data, design = ~condition)
# dds <- DESeq(dds)
# dds.res <- results(dds, contrast = c("condition", 'WT', 'CrxE80A_A'))
# dds.res[c('Tiam2', 'Pde6h'),]
a <- readr::read_csv("sample_info.csv", col_names = TRUE, col_select = c("sample", "condition"), show_col_types = FALSE)
head(a)
expressions <- readr::read_tsv("./test_data/all.fpkm.tsv", show_col_types = FALSE, col_names = TRUE) %>%
    rename_with(function(x) {stringr::str_replace_all(x, "__fpkm", "")}) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sample", values_to = "fpkm") %>%
    dplyr::mutate(fpkm=log2(1+fpkm))
# head(a)
b <- inner_join(a, expressions, by = "sample")
# head(b)
library(ggplot2)
ggplot(data = b) +
    geom_boxplot(mapping = aes(x = sample, y = fpkm, color = condition))+
    ylab("log2(fpkm+1)") +
    xlab("The box plot of gene expression distribution") +
    theme_classic()+
    theme(
        panel.grid.major = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
        axis.title = element_text(size=20)
    )

expre <- readr::read_tsv("./test_data/all.fpkm.tsv", show_col_types = FALSE, col_names= TRUE)
expre <- dplyr::rename_with(expre, function (x){stringr::str_replace_all(x, "__fpkm", "")})
expre

cor(expre)

