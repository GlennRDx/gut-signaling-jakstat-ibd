# Load required libraries
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Load Seurat object with single-cell RNA-seq data
sobj <- readRDS('/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/04_li/GSE266546.rds')

# Pseudobulk the counts based on condition-donor-celltype
pseudo_sobj <- AggregateExpression(sobj, assays = "RNA", return.seurat = T, 
                                   group.by = c("inferred.state", "Sample.name", "annotation_V2"))

# Get count data from pseudobulk object
counts <- GetAssayData(pseudo_sobj, layer = "counts", assay = "RNA")

# Extract metadata from column names
meta_info <- data.frame(colnames = colnames(counts))
meta_info <- meta_info %>% separate(colnames, into = c("inferred.state", "Sample.name", "annotation_V2"), sep = "_")

# Create metadata for DESeq2
meta_data <- data.frame(
  row.names = colnames(counts),
  condition = meta_info$inferred.state,
  sample = meta_info$Sample.name,
  celltype = meta_info$annotation_V2
)

# Get all unique cell types
cell_types <- unique(meta_data$celltype)

# Define results directory path
results_base_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DEG_Results"

# Create results directories
dir.create(results_base_dir, showWarnings = FALSE, recursive = TRUE)
comparisons <- c("CD_vs_Control", "nonInflamed_vs_Control")
for (comparison in comparisons) {
  dir.create(paste0(results_base_dir, "/", comparison), showWarnings = FALSE)
  for (cell_type in cell_types) {
    # Sanitize cell type names to be safe for file paths (replace /, \, and other problematic characters)
    safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
    dir.create(paste0(results_base_dir, "/", comparison, "/", safe_cell_type), showWarnings = FALSE, recursive = TRUE)
  }
}

# Create PCA plot
create_pca_plot <- function(dds, cell_type, comparison, output_dir) {
  vst <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vst, intgroup = c("condition"), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(paste0("PCA Plot: ", cell_type, " - ", comparison)) +
    theme_minimal()
  
  pca_file <- paste0(output_dir, "/PCA_plot.pdf")
  ggsave(pca_file, pca_plot, width = 10, height = 8)
  return(pca_plot)
}

# Create Dispersion plot
create_dispersion_plot <- function(dds, cell_type, comparison, output_dir) {
  disp_file <- paste0(output_dir, "/Dispersion_plot.pdf")
  pdf(disp_file, width = 10, height = 8)
  plotDispEsts(dds, main = paste0("Dispersion Plot: ", cell_type, " - ", comparison))
  dev.off()
}

# Create MA plot
create_ma_plot <- function(res, cell_type, comparison, output_dir) {
  ma_file <- paste0(output_dir, "/MA_plot.pdf")
  pdf(ma_file, width = 10, height = 8)
  plotMA(res, main = paste0("MA Plot: ", cell_type, " - ", comparison), ylim = c(-5, 5), cex = 0.8)
  abline(h = 0, col = "gray", lwd = 2)
  dev.off()
}

# Perform DESeq2 analysis
perform_DESeq2 <- function(counts_matrix, metadata, condition1, condition2, cell_type) {
  # Define comparison name
  if (condition2 == "CD") {
    comparison <- "CD_vs_Control"
  } else if (condition2 == "non-inflamed") {
    comparison <- "nonInflamed_vs_Control"
  } else {
    comparison <- paste0(condition2, "_vs_", condition1)
  }
  
  # Sanitize cell type name for directory creation
  safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
  
  # Create output directory
  output_dir <- paste0(results_base_dir, "/", comparison, "/", safe_cell_type)
  
  # Subset data for the two conditions to compare
  keep <- metadata$condition %in% c(condition1, condition2)
  counts_subset <- counts_matrix[, keep]
  meta_subset <- metadata[keep, ]
  
  # Check if we have enough samples
  if(nrow(meta_subset) < 3) {
    return(NULL)
  }
  
  # Set reference level
  meta_subset$condition <- factor(meta_subset$condition, levels = c(condition1, condition2))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta_subset,
    design = ~ condition
  )
  
  # Filter low count genes
  keep_genes <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_genes, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", condition2, condition1))
  
  # Create diagnostic plots
  create_pca_plot(dds, cell_type, comparison, output_dir)
  create_dispersion_plot(dds, cell_type, comparison, output_dir)
  create_ma_plot(res, cell_type, comparison, output_dir)
  
  # Create volcano plot
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  volcano_file <- paste0(output_dir, "/Volcano_plot.pdf")
  pdf(volcano_file, width = 10, height = 8)
  print(EnhancedVolcano(res_df,
                        lab = res_df$gene,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste0(cell_type, ': ', condition2, ' vs ', condition1),
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        pointSize = 3.0,
                        labSize = 4.0))
  dev.off()
  
  # Save results table
  output_file <- paste0(output_dir, "/DEG_results.csv")
  write.csv(res_df, file = output_file)
  
  return(res_df)
}

# Run DESeq2 for each cell type for both comparisons
results_CD_vs_Control_by_celltype <- list()
results_nonInflamed_vs_Control_by_celltype <- list()
failed_analyses <- list()

for (cell_type in cell_types) {
  # Subset data for this cell type
  cell_type_indices <- meta_data$celltype == cell_type
  counts_cell_type <- counts[, cell_type_indices]
  meta_cell_type <- meta_data[cell_type_indices, ]
  
  # Compare CD vs Control
  cd_result <- tryCatch({
    perform_DESeq2(counts_cell_type, meta_cell_type, "Control", "CD", cell_type)
  }, error = function(e) {
    failed_analyses[[paste0(cell_type, "_CD_vs_Control")]] <- e$message
    return(NULL)
  })
  results_CD_vs_Control_by_celltype[[cell_type]] <- cd_result
  
  # Compare non-inflamed vs Control
  ni_result <- tryCatch({
    perform_DESeq2(counts_cell_type, meta_cell_type, "Control", "non-inflamed", cell_type)
  }, error = function(e) {
    failed_analyses[[paste0(cell_type, "_nonInflamed_vs_Control")]] <- e$message
    return(NULL)
  })
  results_nonInflamed_vs_Control_by_celltype[[cell_type]] <- ni_result
}

# Create summary table
summary_table <- data.frame(
  cell_type = cell_types,
  CD_vs_Control_up = sapply(cell_types, function(ct) {
    res <- results_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res)) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  }),
  CD_vs_Control_down = sapply(cell_types, function(ct) {
    res <- results_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res)) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  }),
  nonInflamed_vs_Control_up = sapply(cell_types, function(ct) {
    res <- results_nonInflamed_vs_Control_by_celltype[[ct]]
    if(is.null(res)) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  }),
  nonInflamed_vs_Control_down = sapply(cell_types, function(ct) {
    res <- results_nonInflamed_vs_Control_by_celltype[[ct]]
    if(is.null(res)) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  })
)

# Save summary table
write.csv(summary_table, paste0(results_base_dir, "/DEG_summary_by_celltype.csv"), row.names = FALSE)

# Print failed analyses
if (length(failed_analyses) > 0) {
  cat("\nCell types that failed to produce DGE results:\n")
  for (name in names(failed_analyses)) {
    cat(paste0(name, ": ", failed_analyses[[name]], "\n"))
  }
}