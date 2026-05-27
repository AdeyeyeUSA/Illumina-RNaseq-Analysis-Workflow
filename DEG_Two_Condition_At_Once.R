# /home/ade/Renemma_2/Analyzed_R/Pairwise_compare/DESeq2_pairwise_ar_reg.R

pacman::p_load(
  "tximeta",
  "AnnotationDbi",
  "GenomicFeatures",
  "txdbmaker",
  "DESeq2",
  "ggplot2",
  "tools",
  "org.Rn.eg.db"
)

gtf_file <- "/home/ade/Renemma_2/Analyzed_R/Pairwise_compare/Rattus_norvegicus.mRatBN7.2.113.gtf"
sample_metadata <- "/home/ade/Renemma_2/Analyzed_R/Pairwise_compare/sample_ar_reg.tsv"

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

if (!file.exists(sample_metadata)) {
  stop("Sample metadata file not found: ", sample_metadata)
}

meta_name <- file_path_sans_ext(basename(sample_metadata))
output_dir <- file.path(dirname(sample_metadata), paste0("DESeq2_results_", meta_name))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

coldata <- read.delim(sample_metadata, stringsAsFactors = TRUE)

required_cols <- c("files", "condition")
missing_cols <- setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in metadata: ", paste(missing_cols, collapse = ", "))
}

coldata$files <- path.expand(as.character(coldata$files))

if (!all(file.exists(coldata$files))) {
  missing_files <- coldata$files[!file.exists(coldata$files)]
  stop("Some Salmon quantification files are missing:\n", paste(missing_files, collapse = "\n"))
}

txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)

k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = k,
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

dds <- DESeq(dds)
res <- results(dds)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

gene_ids_clean <- sub("\\.\\d+$", "", res_df$gene_id)

gene_symbols <- AnnotationDbi::mapIds(
  org.Rn.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df$gene_name <- unname(gene_symbols[gene_ids_clean])
res_df$gene_name <- ifelse(
  is.na(res_df$gene_name) | res_df$gene_name == "",
  res_df$gene_id,
  res_df$gene_name
)

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$gene_id <- rownames(normalized_counts_df)

all_genes_data <- merge(
  res_df,
  normalized_counts_df,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

all_genes_file <- file.path(
  output_dir,
  "Normalized_Counts_All_Genes_with_Statistics_and_GeneNames.csv"
)
write.csv(all_genes_data, all_genes_file, row.names = FALSE)

group_means <- as.data.frame(t(apply(normalized_counts, 1, function(x) {
  tapply(x, coldata$condition, mean, na.rm = TRUE)
})))
group_means$gene_id <- rownames(group_means)

mean_colnames <- colnames(group_means)
mean_colnames[mean_colnames != "gene_id"] <- paste0("mean_", mean_colnames[mean_colnames != "gene_id"])
colnames(group_means) <- mean_colnames

deg_threshold <- 0.05
deg_filtered <- res_df[!is.na(res_df$padj) & res_df$padj < deg_threshold, , drop = FALSE]

if (nrow(deg_filtered) > 0) {
  deg_genes <- deg_filtered$gene_id
  
  deg_counts <- as.data.frame(normalized_counts[deg_genes, , drop = FALSE])
  deg_counts$gene_id <- rownames(deg_counts)
  
  deg_table <- merge(deg_filtered, deg_counts, by = "gene_id", all.x = TRUE, sort = FALSE)
  deg_table <- merge(deg_table, group_means, by = "gene_id", all.x = TRUE, sort = FALSE)
} else {
  deg_table <- deg_filtered
}

deg_file <- file.path(output_dir, "DESeq2_DEGs_Normalized_Counts_with_Means_and_Symbols.csv")
write.csv(deg_table, deg_file, row.names = FALSE)

message("Saved normalized counts for all genes with statistics and gene names to: ", all_genes_file)
message("Saved DEG results with normalized counts, group means, and gene symbols to: ", deg_file)

deg_count <- nrow(deg_filtered)
upregulated <- sum(deg_filtered$log2FoldChange > 0, na.rm = TRUE)
downregulated <- sum(deg_filtered$log2FoldChange < 0, na.rm = TRUE)

message("Summary of Differentially Expressed Genes:")
message("Total DEGs (padj < 0.05): ", deg_count)
message("Upregulated DEGs: ", upregulated)
message("Downregulated DEGs: ", downregulated)
