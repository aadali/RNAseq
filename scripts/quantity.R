library("tibble")
library("ggplot2")



#' Title
#'
#' @param sample_info get the groups of samples
#'
#' @return a tibble of tow columns named as c(sample, condition)
get_samples <- function(sample_info) {
    samples_df <- readr::read_csv(sample_info, comment = "#", show_col_types = FALSE, col_names = TRUE) %>%
        dplyr::select(c("sample", "condition"))
    return(samples_df) # colnames of samples_df: [sample, condition]
}

#' Title
#'
#' @param raw_counts raw counts vector
#' @param genes_length lengths of genes vector, its length must be equal raw_counts
#' @param method tpm or fpkm
#'
#' @return expression of tpm or fpkm vector
calculate_expression <- function(raw_counts, genes_length, method = "fpkm") {
    if (!method %in% c("fpkm", "tpm")) {
        stop("method must be one of ('fpkm', 'tpm')")
    }
    if (method == "fpkm") {
        return(raw_counts * 10^9 / (genes_length * sum(raw_counts)))
    } else {
        return(raw_counts / genes_length * 11^6 / sum(raw_counts / genes_length))
    }
}

#' Title
#'
#' @param count_files count file path from featureCount, character vector
#' @param outdir the output directory
#' @param analysis_name the prefix of output file name
#' merge count files of samples from featureCount into one tsv file whose columns will be
#' c("Geneid", "Length", sample1, sample2, smaple3...), file will be named as {analysis_name}.all.raw.counts.tsv
#' And return this data.frame
#' @return a tibble of all samples' raw count
#'
#' @examples
get_counts <- function(count_files, outdir, analysis_name) {
    count_df_list <- purrr::map(count_files, function(x) {
        readr::read_tsv(x, show_col_types = FALSE, comment = "#") %>%
            dplyr::select(c(1, 6:length(.))) %>%
            dplyr::rename_with(.cols = dplyr::ends_with(".sorted.bam"),
                               .fn = function(y) {stringr::str_replace_all(y, ".sorted.bam|.+/", "")})}
    )
    counts_df <- purrr::reduce(count_df_list, function(x, y) {
        dplyr::inner_join(x, y, by = c("Geneid", "Length"))
    })
    readr::write_tsv(counts_df, file.path(outdir, paste0(analysis_name, ".all.raw.counts.tsv")))
    return(counts_df)
}

#' Title
#'
#' @param counts_df the raw counts of all  samples' from get_counts
#' @param samples_df the group infomation of samples from get_samples
#' @param method fpkm or tpm that will be used in get_expression
#' @param outdir the output directory
#' @param analysis_name prefix name of output file
#' read counts_df from get_counts, and calculate expression in place, 
#' but rename the sample name with adding the __{method} suffix
#' @return a tibble object of expresion whose colnames are c("Geneid", "Length", sample1__{method}, sample2__{method}, ...)
get_and_save_expression <- function(counts_df, samples_df, method, outdir, analysis_name) {
    # colnames of counts_df: Geneid, Length, sample1, sample2, sample3...
    keep <- rowSums(counts_df[, 3:length(counts_df)] > 5) > 3
    expression <- counts_df[keep,] %>%
        dplyr::mutate(dplyr::across(samples_df$sample, function(x) {
            calculate_expression(.$x, .$Length, method = method)
        })) %>%
        dplyr::rename_with(function(x) {paste0(x, "__", method)}, 3:length(.)) # add suffix: __{method} for each sample name
    readr::write_tsv(expression, file.path(outdir, paste0(analysis_name, ".all.expression.tsv")))
    return(expression)
}

#' Title
#'
#' @param expression expression tibble from get_and_save_expression
#' @param method fpkm or tpm
#' @param samples_df group of samples
#' @param outdir output directory
#' @param plot_name the prefix of plot name
#' @param boxplot whether make boxplot plot?
#' @param heatmap whether make heatmap plot?
#' @param pca whether make pca plot?
#' make some figures, such as boxplot of expression for all samples, sampels2samples correlation heatmap 
#' and pca plot of different groups. Save them in disk
#' @return NULL
#' @export
#'
#' @examples
plot <- function(expression, 
                 method, 
                 samples_df, 
                 outdir, 
                 plot_name, 
                 boxplot = TRUE, 
                 heatmap = TRUE, 
                 pca = TRUE) {
    expression_df <- dplyr::rename_with(expression,
                                        .fn = function(x) {stringr::str_replace_all(x, paste0("__", method), "")},
                                        .cols = 3:length(expression))
    # colnames of expression_df: sample1, sample2, sample3...
    expression_df <- dplyr::select(expression_df, 3:length(expression_df)) %>%
        dplyr::rename_with(.fn = function(x) {stringr::str_replace_all(x, paste0("__", method), "")}) # restore sample names
    
    if (boxplot) {
        expression_long <- tidyr::pivot_longer(expression_df,
                                               cols = dplyr::everything(),
                                               names_to = "sample",
                                               values_to = "expression")
        expression_samples_long <- dplyr::inner_join(expression_long, samples_df, by = 'sample')
        expression_boxplot <- ggplot(expression_samples_long) +
            geom_boxplot(mapping = aes(x = sample, y = log2(expression + 1), color = condition)) +
            labs(color = "Group", x = NULL, y = paste0("log2(", method, "+1)"), title = "The boxplot of gene expression") +
            scale_color_brewer(type = "qual", palette = "Dark2") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5),
                  axis.text.y = element_text(size = 14),
                  plot.title = element_text(hjust = 0.5, size = 18, margin = margin(t = 20, b = 10)),
            )
        for (fmt in c('pdf', 'svg')){
            ggsave(file.path(outdir, paste0(plot_name, ".boxplot.", fmt)), plot = expression_boxplot, height = 10, width = 10, units = "in")
        }
    }

    if (heatmap) {
        corr_df <- cor(log2(expression_df + 0.01))^2
        coor_long_df <- data.frame(corr_df) %>%
            tibble::rownames_to_column(var = "row") %>%
            tidyr::pivot_longer(cols = colnames(.)[2:length(.)], 
                                names_to = "col", 
                                values_to = "R2")
        samples_heatmap <- ggplot(coor_long_df) +
            geom_tile(mapping = aes(row, col, fill = R2)) +
            scale_fill_gradient(low = "white", high = "#bd0026") +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            geom_text(mapping = aes(row, col, label = round(R2, 3)), size = 3) +
            labs(fill = quote(R^2), x = NULL, y = NULL, title = "Correlation of samples to samples") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
                  axis.ticks = element_blank(),
                  plot.title = element_text(hjust = 0.5, size = 14))
        for (fmt in c('pdf', 'svg')) {
            ggsave(file.path(outdir, paste0(plot_name, ".heatmap.of.samples2samples.", fmt)), plot = samples_heatmap, width = 10, height = 10, units = "in")
        }
    }

    if (pca) {
        # use SDV to make PCA
        pca_res <- prcomp(t(expression_df), scale. = TRUE, center = TRUE)
        # calculate the contribution of each component to total variance,
        # but only get the large 2 sdev that means we get the first and second element
        sdev_ratio <- (pca_res$sdev^2/sum(pca_res$sdev^2))[1:2]
        pc_dim2 <- as.data.frame(pca_res$x[,c('PC1', 'PC2')]) %>% 
            tibble::rownames_to_column("sample") %>% 
            dplyr::inner_join(samples_df, by="sample")
        pca_plot <- ggplot(data=pc_dim2) + 
            geom_point(mapping = aes(PC1, PC2, color=condition)) + 
            theme_classic() + 
            labs(
                x=paste0("PC1 (", round(sdev_ratio[1]*100, 2), "%)"),
                y=paste0("PC2 (", round(sdev_ratio[2]*100,2), "%)"),
                color="group"
            )
        for (fmt in c('pdf', 'svg')){
            ggsave(file.path(outdir, paste0(plot_name, ".pca", fmt)), plot = pca_plot)
        }
    }
}

main <- function() {
    args <- commandArgs(TRUE)
    samples <- args[1]
    # comparisons <- args[2]
    outdir <- args[2]
    method <- args[3]
    analysis_name <- args[4]
    count_files <- args[5]
    if (length(args) != 5) {
        stop(
            "Useage:
            Rscript quantity.R <samples> <comparisons> <outdir> <method> <analysis_name> <count_files>
            --------------------------------------------------------------------------
            samples: the sample_info.csv
            comparisons: the comparisons.csv
            outdir: the output directory
            method: fpkm or tpm
            analysis_name: prefix of files or plots
            count_files: one or more count file from featureCounts, seprated by comma
            ")
    }

    
    # print(getwd())
    # samples <- "gaopeng_samples.csv"
    # comparisons <- "comparisons.csv"
    # outdir <- "."
    # method <- "fpkm"
    # count_files <- list.files("./test_data/gaopeng/", full.names = TRUE)
    
    count_files <- stringr::str_split(count_files, ",")[[1]]
    if (!method %in% c("fpkm", "tpm")) {stop("method must be one of ('fpkm', 'tpm')")}
    samples_df <- get_samples(samples)
    counts_df <- get_counts(
        count_files = count_files, 
        analysis_name = analysis_name, 
        outdir = outdir)
    
    expression <- get_and_save_expression(
        counts_df = counts_df, 
        samples_df = samples_df, 
        method = method, 
        analysis_name = analysis_name, 
        outdir = outdir)
    plot(
        expression = expression, 
        method = method, 
        samples_df = samples_df, 
        outdir = outdir, 
        plot_name = analysis_name,
        boxplot = TRUE,
        heatmap = TRUE,
        pca = TRUE
        )
}

main()
