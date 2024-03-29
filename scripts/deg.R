library(optparse)
library(readr)
library(dplyr)
library(purrr)
library(DESeq2)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)

count_file <- "/home/a/pub/ycq/output/gaopeng2/04.quant/gaopeng2.all.raw.counts.tsv"
samples_info <- "/home/a/big/ycq/projects/RNAseq/gaopeng_samples.csv"
comparisons_info <- "/home/a/big/ycq/projects/RNAseq/comparisons.csv"
outdir <- "./test_data"
min_lfc <- 1
alpha <- 0.05
shrink <- FALSE
outdir <- "./test_data"


setwd("/home/a/big/ycq/projects/RNAseq")
# which combination will be analysised different genes, this is contrast param of DESeq2::results()
comparisons_df <- read_csv(comparisons_info, show_col_types = FALSE, comment = "#")

# contrast: c('condition', 'treated', 'control')
contrasts <- map2(comparisons_df$treated, comparisons_df$control, function(x, y) {c("condition", x, y)})

# get colData param of DESeqDataSetFromMatrix
col_data <- read_csv(samples_info, show_col_types = FALSE, comment = "#") %>%
    select("condition", "sample") %>%
    arrange(condition, sample) %>%
    column_to_rownames("sample")

# get raw countData param of DESeqDataSetFromMatrix
count_df <- read_tsv(count_file, show_col_types = FALSE, comment = "#") %>%
    select(!"Length") %>%
    column_to_rownames("Geneid")

# pre filt countData
count_df <- count_df[keep, rownames(col_data)]

col_data$condition <- factor(col_data$condition)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_df, colData = col_data, design = ~condition)
dds <- DESeq2::DESeq(dds)

# normalized count from DESeq2
normalized_count_df <- DESeq2::counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("geneid")
deg_info_cols <- c("log2FoldChange", "pvalue", "padj")

# to save number of up or down genes of each comparison
all_comparisons_deg_info <- tibble(comparison=character(), up=integer(), down=integer(), no=integer())
for (contrast in contrasts) {
    # this comparison deg results
    result <- DESeq2::results(dds, contrast = contrast, lfcThreshold = min_lfc, alpha = alpha)
    result_df <- data.frame(result) %>% rownames_to_column("geneid")
    this_group_samples <- data.frame(colData(dds)) %>% filter(condition %in% contrast[2:3])

    # the output tsv columns: c("geneid", sample1, sample2..., log2FoldChange, pvalue, padj)
    # TODO add additional columns about genes annotation
    result_columns <- c("geneid", rownames(this_group_samples), deg_info_cols)

    comparison_result <- left_join(normalized_count_df, result_df, by = "geneid") %>%
        select(result_columns)
    comparison_result <- na.omit(comparison_result)

    # export result into specified directory with specified file name. {TreatedName}__VS__{ControlName}
    result_name_prefix <- paste0(contrast[2], "__VS__", contrast[3])
    comparison_name <- str_glue("{treated}__VS__{control}", treated = contrast[2], control = contrast[3])
    print(str_glue("{outdir}/{comparison_name}"))
    comparison_outdir <- str_glue("{outdir}/{comparison_name}")
    prefix <- str_glue("{comparison_outdir}/{comparison_name}")
    dir.create(comparison_outdir, showWarnings = FALSE)

    # save the result of deseq2
    write_tsv(comparison_result, file = str_glue("{prefix}.deg2_results.tsv")) # TODO
    comparison_deg <- filter(comparison_result, abs(log2FoldChange) > min_lfc & padj < alpha) %>%
        select(result_columns)
    # save all deg genes
    write_tsv(comparison_deg, paste0(prefix, ".all_deg.tsv"))
    # save up deg genes
    comparison_up_deg <- filter(comparison_deg, log2FoldChange > 1)
    write_tsv(comparison_up_deg, paste0(prefix, ".up_deg.tsv"))
    # save down deg genes
    comparison_down_deg <- filter(comparison_deg, log2FoldChange < -1)
    write_tsv(comparison_down_deg, paste0(prefix, ".down_deg.tsv"))

    comparison_result <- mutate(comparison_result,
                                change = case_when(log2FoldChange > min_lfc & padj < alpha ~ "UP",
                                                   log2FoldChange < alpha & padj < alpha ~ "DOWN",
                                                   .default = "NO"),
                                log10padj = log10(padj))
    up_number <- sum(comparison_result$change == "UP")
    down_number <- sum(comparison_result$change == "DOWN")
    # some low counts will be dropped. this genes will be counted for  no changing genes
    no_number <- dim(comparison_result)[1] - up_number - down_number
    all_comparisons_deg_info <- add_row(
        all_comparisons_deg_info,
        comparison = str_glue("{treated}/{control}", treated=contrast[2], control=contrast[3]),
        up = up_number,
        down = down_number,
        no = no_number
    )
    next

    # make volcano plot
    volcano <- ggplot(comparison_result) +
        geom_point(mapping = aes(x = log2FoldChange, y = -log10padj, color = change)) +
        geom_vline(xintercept = c(-min_lfc, min_lfc), linetype = "dashed", color = "grey") +
        geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey") +
        scale_y_continuous(breaks = c(-log10(alpha), 20, 50, 80, 100, 150),
                           labels = c(round(-log10(alpha), 2), 20, 50, 80, 100, 150)) +
        theme_classic() +
        scale_color_manual(values = c("UP" = "red", "DOWN" = "green", "NO" = "grey"),
                           breaks = c("UP", "DOWN", "NO"),
                           labels = c("UP" = paste0("UP (", up_number, ")"),
                                      "DOWN" = paste0("DOWN (", down_number, ")"),
                                      "NO" = paste0("NO (", no_number, ")"))) +
        labs(x = "log2FoldChange", y = "-log10(padj)", color = "group")
    ggsave(str_glue("{prefix}.volcona.pdf"), plot = volcano, height = 10, width = 6, units = 'in')
    ggsave(str_glue("{prefix}.volcona.png"), plot = volcano, height = 10, width = 6, units = 'in')
    ggsave(str_glue("{prefix}.volcona.svg"), plot = volcano, height = 10, width = 6, units = 'in')
}


if (shrink) {
    deg_res <- lfcShrink(dds = dds, contrast = contrasts[[1]], type = "normal")
}
if (1) {
    a <- pivot_longer(all_comparisons_deg_info, cols = c('up', 'down', 'no'), names_to = 'change', values_to = 'number')
}
head(a)
ggplot(a) +
    geom_bar(mapping = aes(comparison, number, fill=change), stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90))

apply(all_comparisons_deg_info[,2:4], 1, sum)
b <- mutate(all_comparisons_deg_info, total=down + no + up)
little_com <- DESeq2::results(dds, contrast = c('condition', 'TTL', 'TCL'), lfcThreshold = min_lfc, alpha = alpha)
little_com
summary(little_com)
dim(little_com)

log2(4.33/0.49)
print("hello world")