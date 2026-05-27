# /home/ade/Renemma_2/Analyzed_R/Pairwise_compare/run_all_pairwise_deseq2.R

pacman::p_load(
  "tximeta",
  "AnnotationDbi",
  "GenomicFeatures",
  "txdbmaker",
  "DESeq2",
  "tools",
  "org.Rn.eg.db"
)

gtf_file <- "/home/ade/Renemma_2/Analyzed_R/Pairwise_compare/Rattus_norvegicus.mRatBN7.2.113.gtf"
sample_metadata <- "/home/ade/Renemma_2/Analyzed_R/Pairwise_compare/sample.tsv"
base_output_dir <- "/home/ade/Renemma_2/Analyzed_R/Pairwise_compare_New"

deg_threshold <- 0.1

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

if (!file.exists(sample_metadata)) {
  stop("Sample metadata file not found: ", sample_metadata)
}

if (!dir.exists(base_output_dir)) {
  dir.create(base_output_dir, recursive = TRUE)
}

coldata_all <- read.delim(sample_metadata, stringsAsFactors = FALSE)

required_cols <- c("files", "condition")
missing_cols <- setdiff(required_cols, colnames(coldata_all))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in metadata: ", paste(missing_cols, collapse = ", "))
}

coldata_all$files <- path.expand(as.character(coldata_all$files))
coldata_all$condition <- as.character(coldata_all$condition)

if (!all(file.exists(coldata_all$files))) {
  missing_files <- coldata_all$files[!file.exists(coldata_all$files)]
  stop("Some Salmon quantification files are missing:\n", paste(missing_files, collapse = "\n"))
}

all_conditions <- sort(unique(coldata_all$condition))

if (length(all_conditions) < 2) {
  stop("Need at least two unique conditions for pairwise comparisons.")
}

message("Building transcript-to-gene map once from GTF...")
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

map_gene_symbols <- function(gene_ids) {
  gene_ids_clean <- sub("\\.\\d+$", "", gene_ids)
  
  gene_symbols <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = gene_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  gene_names <- unname(gene_symbols[gene_ids_clean])
  gene_names[is.na(gene_names) | gene_names == ""] <- gene_ids[is.na(gene_names) | gene_names == ""]
  gene_names
}

run_pairwise_comparison <- function(condition_1, condition_2, coldata_all, tx2gene, base_output_dir, deg_threshold) {
  message("\n==============================")
  message("Running pairwise comparison: ", condition_1, " vs ", condition_2)
  message("==============================")
  
  pair_coldata <- coldata_all[coldata_all$condition %in% c(condition_1, condition_2), , drop = FALSE]
  
  if (nrow(pair_coldata) == 0) {
    stop("No samples found for comparison: ", condition_1, " vs ", condition_2)
  }
  
  pair_coldata$condition <- factor(pair_coldata$condition, levels = c(condition_1, condition_2))
  
  pair_counts <- table(pair_coldata$condition)
  if (any(pair_counts == 0)) {
    stop("One condition has zero samples in comparison: ", condition_1, " vs ", condition_2)
  }
  
  pair_dir <- file.path(base_output_dir, paste0("Pairwise_", condition_1, "_", condition_2))
  if (!dir.exists(pair_dir)) {
    dir.create(pair_dir, recursive = TRUE)
  }
  
  pair_metadata_file <- file.path(pair_dir, paste0("metadata_", condition_1, "_vs_", condition_2, ".tsv"))
  write.table(pair_coldata, pair_metadata_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  gse <- tximeta(
    pair_coldata,
    skipMeta = TRUE,
    txOut = FALSE,
    tx2gene = tx2gene,
    ignoreAfterBar = TRUE,
    ignoreTxVersion = TRUE
  )
  
  dds <- DESeqDataSet(gse, design = ~ condition)
  dds$condition <- droplevels(dds$condition)
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("condition", condition_2, condition_1))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- map_gene_symbols(res_df$gene_id)
  
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
  
  group_means <- as.data.frame(t(apply(normalized_counts, 1, function(x) {
    tapply(x, pair_coldata$condition, mean, na.rm = TRUE)
  })))
  group_means$gene_id <- rownames(group_means)
  
  mean_colnames <- colnames(group_means)
  mean_colnames[mean_colnames != "gene_id"] <- paste0("mean_", mean_colnames[mean_colnames != "gene_id"])
  colnames(group_means) <- mean_colnames
  
  all_genes_data <- merge(
    all_genes_data,
    group_means,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  all_genes_file <- file.path(
    pair_dir,
    paste0("All_Genes_", condition_1, "_vs_", condition_2, "_with_statistics.csv")
  )
  write.csv(all_genes_data, all_genes_file, row.names = FALSE)
  
  deg_filtered <- res_df[!is.na(res_df$padj) & res_df$padj < deg_threshold, , drop = FALSE]
  
  if (nrow(deg_filtered) > 0) {
    deg_counts <- as.data.frame(normalized_counts[deg_filtered$gene_id, , drop = FALSE])
    deg_counts$gene_id <- rownames(deg_counts)
    
    deg_table <- merge(deg_filtered, deg_counts, by = "gene_id", all.x = TRUE, sort = FALSE)
    deg_table <- merge(deg_table, group_means, by = "gene_id", all.x = TRUE, sort = FALSE)
  } else {
    deg_table <- deg_filtered
  }
  
  deg_file <- file.path(
    pair_dir,
    paste0("DEGs_", condition_1, "_vs_", condition_2, ".csv")
  )
  write.csv(deg_table, deg_file, row.names = FALSE)
  
  summary_df <- data.frame(
    comparison = paste0(condition_1, "_vs_", condition_2),
    reference_condition = condition_1,
    test_condition = condition_2,
    total_samples = nrow(pair_coldata),
    samples_reference = sum(pair_coldata$condition == condition_1),
    samples_test = sum(pair_coldata$condition == condition_2),
    total_genes_tested = nrow(res_df),
    deg_threshold = deg_threshold,
    total_degs = nrow(deg_filtered),
    upregulated_in_test_condition = sum(deg_filtered$log2FoldChange > 0, na.rm = TRUE),
    downregulated_in_test_condition = sum(deg_filtered$log2FoldChange < 0, na.rm = TRUE)
  )
  
  summary_file <- file.path(
    pair_dir,
    paste0("Summary_", condition_1, "_vs_", condition_2, ".csv")
  )
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  message("Completed: ", condition_1, " vs ", condition_2)
  message("Results folder: ", pair_dir)
  
  summary_df
}

condition_pairs <- combn(all_conditions, 2, simplify = FALSE)

all_summaries <- lapply(condition_pairs, function(pair) {
  run_pairwise_comparison(
    condition_1 = pair[1],
    condition_2 = pair[2],
    coldata_all = coldata_all,
    tx2gene = tx2gene,
    base_output_dir = base_output_dir,
    deg_threshold = deg_threshold
  )
})

combined_summary <- do.call(rbind, all_summaries)

combined_summary_file <- file.path(base_output_dir, "Pairwise_Comparison_Summary.csv")
write.csv(combined_summary, combined_summary_file, row.names = FALSE)

message("\nAll pairwise comparisons completed.")
message("Combined summary saved to: ", combined_summary_file)

