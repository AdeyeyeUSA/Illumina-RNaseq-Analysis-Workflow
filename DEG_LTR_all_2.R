# ~/Renemma_2/Analyzed_R/DEG_LTR_all.R

pacman::p_load(
  "tximeta",
  "AnnotationDbi",
  "GenomicFeatures",
  "txdbmaker",
  "DESeq2",
  "ggplot2"
)

gtf_file <- "Rattus_norvegicus.mRatBN7.2.113.gtf"
sample_metadata <- "sample_original1.tsv"
output_dir <- "DESeq2_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

if (!file.exists(sample_metadata)) {
  stop("Sample metadata file not found: ", sample_metadata)
}

coldata <- read.delim(sample_metadata, stringsAsFactors = TRUE)

required_cols <- c("files", "condition")
missing_cols <- setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in sample metadata: ", paste(missing_cols, collapse = ", "))
}

coldata$files <- as.character(coldata$files)

if (!all(file.exists(coldata$files))) {
  missing_files <- coldata$files[!file.exists(coldata$files)]
  stop("Some Salmon quantification files are missing:\n", paste(missing_files, collapse = "\n"))
}

txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)

tx_ids <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = tx_ids,
  columns = "GENEID",
  keytype = "TXNAME"
)

tx2gene <- unique(tx2gene)
tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), , drop = FALSE]
tx2gene$TXNAME <- sub("\\.\\d+$", "", tx2gene$TXNAME)

gse <- tximeta(
  coldata,
  skipMeta = TRUE,
  txOut = FALSE,
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE,
  ignoreTxVersion = TRUE
)

dds <- DESeqDataSet(gse, design = ~ condition)
dds$condition <- droplevels(dds$condition)

dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds)
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df$gene_id <- rownames(res_lrt_df)

lrt_file <- file.path(output_dir, "DESeq2_results_LRT.csv")
write.csv(res_lrt_df, lrt_file, row.names = FALSE)

deg_threshold <- 0.05
num_deg <- sum(!is.na(res_lrt_df$padj) & res_lrt_df$padj < deg_threshold)

message("Number of DEGs with padj < ", deg_threshold, ": ", num_deg)
message("LRT results saved to: ", lrt_file)

vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_line(color = "grey80")
  )

print(p)

pca_file <- file.path(output_dir, "PCA_plot.png")
ggsave(pca_file, plot = p, width = 8, height = 6, dpi = 300)

message("PCA plot saved to: ", pca_file)