# ~/Renemma_2/Analyzed_R/heatmap_all_significant_genes.R

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(grid)

deg_threshold <- 0.05
significant_genes <- rownames(res_lrt)[!is.na(res_lrt$padj) & res_lrt$padj < deg_threshold]

if (length(significant_genes) == 0) {
  stop("No significant genes found with the given criteria.")
}

message("Performing variance stabilizing transformation...")
vsd <- vst(dds, blind = FALSE)
sig_vsd_counts <- assay(vsd)[significant_genes, , drop = FALSE]

cat("Dimensions of VST counts for significant genes:", dim(sig_vsd_counts), "\n")

gene_sd <- apply(sig_vsd_counts, 1, sd, na.rm = TRUE)
sig_vsd_counts <- sig_vsd_counts[gene_sd > 0, , drop = FALSE]

if (nrow(sig_vsd_counts) == 0) {
  stop("All significant genes have zero variance after VST; cannot create heatmap.")
}

scaled_sig_vsd_counts <- t(scale(t(sig_vsd_counts), center = TRUE, scale = TRUE))
scaled_sig_vsd_counts[is.na(scaled_sig_vsd_counts)] <- 0

cat("Range of scaled VST counts:\n")
print(range(scaled_sig_vsd_counts, na.rm = TRUE))

desired_order <- c("reg", "ar", "ssar", "esar")
sample_condition <- factor(as.character(colData(dds)$condition), levels = desired_order)

if (any(is.na(sample_condition))) {
  stop(
    "Some samples have condition values not found in desired_order.\nFound: ",
    paste(sort(unique(as.character(colData(dds)$condition))), collapse = ", ")
  )
}

annotation_col <- data.frame(
  Condition = sample_condition,
  row.names = colnames(scaled_sig_vsd_counts)
)

stopifnot(identical(rownames(annotation_col), colnames(scaled_sig_vsd_counts)))

sample_order <- order(annotation_col$Condition)
scaled_sig_vsd_counts <- scaled_sig_vsd_counts[, sample_order, drop = FALSE]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

breakpoints <- seq(-2, 2, length.out = 9)
heat_colors <- rev(c(
  "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF",
  "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4"
))
col_fun <- colorRamp2(breakpoints, heat_colors)

condition_colors <- c(
  reg = "#4DAF4A",
  ar = "#377EB8",
  ssar = "#984EA3",
  esar = "#E41A1C"
)

ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(Condition = condition_colors),
  annotation_name_gp = gpar(fontsize = 10),
  show_legend = FALSE
)

ht <- Heatmap(
  scaled_sig_vsd_counts,
  name = "Scaled Expression",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = annotation_col$Condition,
  cluster_column_slices = FALSE,
  top_annotation = ha,
  show_column_dend = FALSE,
  heatmap_legend_param = list(
    title = "Scaled Expression Levels",
    color_bar = "continuous"
  ),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 6)
)

png_file <- file.path(output_dir, "heatmap_all_significant_genes.png")
png(filename = png_file, width = 3000, height = 2500, res = 600, bg = "white")
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

pdf_file <- file.path(output_dir, "heatmap_all_significant_genes.pdf")
pdf(file = pdf_file, width = 12, height = 10, useDingbats = FALSE, bg = "white")
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")

message("Heatmap PNG saved to: ", png_file)
message("Heatmap PDF saved to: ", pdf_file)