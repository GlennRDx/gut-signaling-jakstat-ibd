# Load required libraries
library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#' Perform DESeq2 differential expression analysis on single-cell RNA-seq data
#'
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param condition_col Column name in metadata for condition/treatment groups
#' @param sample_col Column name in metadata for sample identifiers
#' @param celltype_col Column name in metadata for cell type annotations
#' @param condition1 Reference condition (typically control group)
#' @param condition2 Test condition (typically treatment group)
#' @param output_dir Base directory path to save analysis results
#' @param return_results Logical, whether to return results in memory (default: FALSE)
#' @param folder_prefix Optional prefix to add to the comparison folder name (default: NULL)
#'
#' @return If return_results = TRUE, returns list containing:
#'   - results_by_celltype: DESeq2 results for each cell type
#'   - summary_table: Summary of up/down regulated genes per cell type
#'   - failed_analyses: List of cell types that failed analysis with error messages
#'   If return_results = FALSE, returns NULL (results only saved to files)
#'
#' @examples
#' # Example usage (files only):
#' run_deseq2_comparison(
#'   seurat_obj = sobj,
#'   condition_col = "inferred.state",
#'   sample_col = "Sample.name", 
#'   celltype_col = "annotation_V2",
#'   condition1 = "Control",
#'   condition2 = "CD",
#'   output_dir = "/path/to/results"
#' )
#'
#' # Example usage (with in-memory results):
#' results <- run_deseq2_comparison(
#'   seurat_obj = sobj,
#'   condition_col = "inferred.state",
#'   sample_col = "Sample.name", 
#'   celltype_col = "annotation_V2",
#'   condition1 = "Control",
#'   condition2 = "CD",
#'   output_dir = "/path/to/results",
#'   return_results = TRUE
# Example usage with folder prefix:
# run_deseq2_comparison(
#   seurat_obj = sobj2,
#   condition_col = "CD_Status",
#   sample_col = "PatientID",
#   celltype_col = "Celltypes",
#   condition1 = "Control",
#   condition2 = "Active_CD",
#   output_dir = "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/",
#   folder_prefix = "li"
# )
run_deseq2_comparison <- function(seurat_obj, condition_col, sample_col, celltype_col, 
                                  condition1, condition2, output_dir, return_results = FALSE,
                                  folder_prefix = NULL) {
  
  # Create comparison name for directory structure
  comparison_name <- paste0(condition2, "_vs_", condition1)
  
  # Add prefix if provided
  if (!is.null(folder_prefix) && folder_prefix != "") {
    comparison_name <- paste0(folder_prefix, "_", comparison_name)
  }
  
  # Pseudobulk the counts based on condition-donor-celltype
  cat("Performing pseudobulk aggregation...\n")
  pseudo_sobj <- AggregateExpression(seurat_obj, assays = "RNA", return.seurat = TRUE, 
                                     group.by = c(condition_col, sample_col, celltype_col))
  
  # Get count data from pseudobulk object
  counts <- GetAssayData(pseudo_sobj, layer = "counts", assay = "RNA")
  
  # Extract metadata from column names
  meta_info <- data.frame(colnames = colnames(counts))
  meta_info <- meta_info %>% 
    separate(colnames, into = c("condition", "sample", "celltype"), sep = "_")
  
  # Create metadata for DESeq2
  meta_data <- data.frame(
    row.names = colnames(counts),
    condition = meta_info$condition,
    sample = meta_info$sample,
    celltype = meta_info$celltype
  )
  
  # Get all unique cell types
  cell_types <- unique(meta_data$celltype)
  cat("Found", length(cell_types), "cell types:\n")
  cat(paste(cell_types, collapse = ", "), "\n\n")
  
  # Create results directories
  results_base_dir <- file.path(output_dir, comparison_name)
  dir.create(results_base_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (cell_type in cell_types) {
    # Sanitize cell type names to be safe for file paths
    safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
    dir.create(file.path(results_base_dir, safe_cell_type), showWarnings = FALSE, recursive = TRUE)
  }
  
  # Helper function to create PCA plot
  create_pca_plot <- function(dds, cell_type, comparison, output_dir) {
    vst_data <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vst_data, intgroup = c("condition"), returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))
    
    pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
      geom_point(size = 4) +
      xlab(paste0("PC1: ", percent_var[1], "% variance")) +
      ylab(paste0("PC2: ", percent_var[2], "% variance")) +
      ggtitle(paste0("PCA Plot: ", cell_type, " - ", comparison)) +
      theme_minimal()
    
    pca_file <- file.path(output_dir, "PCA_plot.pdf")
    ggsave(pca_file, pca_plot, width = 10, height = 8)
    return(pca_plot)
  }
  
  # Helper function to create Dispersion plot
  create_dispersion_plot <- function(dds, cell_type, comparison, output_dir) {
    disp_file <- file.path(output_dir, "Dispersion_plot.pdf")
    pdf(disp_file, width = 10, height = 8)
    plotDispEsts(dds, main = paste0("Dispersion Plot: ", cell_type, " - ", comparison))
    dev.off()
  }
  
  # Helper function to create MA plot
  create_ma_plot <- function(res, cell_type, comparison, output_dir) {
    ma_file <- file.path(output_dir, "MA_plot.pdf")
    pdf(ma_file, width = 10, height = 8)
    plotMA(res, main = paste0("MA Plot: ", cell_type, " - ", comparison), ylim = c(-5, 5), cex = 0.8)
    abline(h = 0, col = "gray", lwd = 2)
    dev.off()
  }
  
  # Helper function to perform DESeq2 analysis for a single cell type
  perform_DESeq2_celltype <- function(counts_matrix, metadata, condition1, condition2, cell_type) {
    # Sanitize cell type name for directory creation
    safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
    
    # Create output directory
    output_dir <- file.path(results_base_dir, safe_cell_type)
    
    # Subset data for the two conditions to compare
    keep <- metadata$condition %in% c(condition1, condition2)
    counts_subset <- counts_matrix[, keep]
    meta_subset <- metadata[keep, ]
    
    # Check if we have enough samples
    if(nrow(meta_subset) < 3) {
      stop("Not enough samples (minimum 3 required)")
    }
    
    # Set reference level (condition1 as reference)
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
    create_pca_plot(dds, cell_type, comparison_name, output_dir)
    create_dispersion_plot(dds, cell_type, comparison_name, output_dir)
    create_ma_plot(res, cell_type, comparison_name, output_dir)
    
    # Create volcano plot
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    volcano_file <- file.path(output_dir, "Volcano_plot.pdf")
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
    output_file <- file.path(output_dir, "DEG_results.csv")
    write.csv(res_df, file = output_file, row.names = FALSE)
    
    return(res_df)
  }
  
  # Run DESeq2 for each cell type
  cat("Running DESeq2 analysis for each cell type...\n")
  results_by_celltype <- list()
  failed_analyses <- list()
  
  for (cell_type in cell_types) {
    cat("Processing cell type:", cell_type, "\n")
    
    # Subset data for this cell type
    cell_type_indices <- meta_data$celltype == cell_type
    counts_cell_type <- counts[, cell_type_indices]
    meta_cell_type <- meta_data[cell_type_indices, ]
    
    # Perform DESeq2 analysis
    result <- tryCatch({
      perform_DESeq2_celltype(counts_cell_type, meta_cell_type, condition1, condition2, cell_type)
    }, error = function(e) {
      failed_analyses[[cell_type]] <<- e$message
      cat("  Failed:", e$message, "\n")
      return(NULL)
    })
    
    results_by_celltype[[cell_type]] <- result
  }
  
  # Create summary table
  cat("Creating summary table...\n")
  summary_table <- data.frame(
    cell_type = cell_types,
    upregulated = sapply(cell_types, function(ct) {
      res <- results_by_celltype[[ct]]
      if(is.null(res)) return(NA)
      sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
    }),
    downregulated = sapply(cell_types, function(ct) {
      res <- results_by_celltype[[ct]]
      if(is.null(res)) return(NA)
      sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
    }),
    total_significant = sapply(cell_types, function(ct) {
      res <- results_by_celltype[[ct]]
      if(is.null(res)) return(NA)
      sum(res$padj < 0.05, na.rm = TRUE)
    })
  )
  
  # Add comparison information to summary
  summary_table$comparison <- comparison_name
  summary_table$condition1 <- condition1
  
  # Add sample counts for condition1
  summary_table$condition1_sample_count <- sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    meta_cell_type <- meta_data[cell_type_indices, ]
    sum(meta_cell_type$condition == condition1)
  })
  
  summary_table$condition2 <- condition2
  
  # Add sample counts for condition2
  summary_table$condition2_sample_count <- sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    meta_cell_type <- meta_data[cell_type_indices, ]
    sum(meta_cell_type$condition == condition2)
  })
  
  # Add gene-specific columns for JAK1, STAT1, STAT2, SOCS1
  target_genes <- c("JAK1", "JAK2", "STAT1", "SOCS1")
  
  for (gene in target_genes) {
    # Add padj column
    summary_table[[paste0(gene, "_padj")]] <- sapply(cell_types, function(ct) {
      res <- results_by_celltype[[ct]]
      if(is.null(res)) return(NA)
      gene_idx <- which(res$gene == gene)
      if(length(gene_idx) == 0) return(NA)
      return(res$padj[gene_idx])
    })
    
    # Add log2FC column
    summary_table[[paste0(gene, "_log2FC")]] <- sapply(cell_types, function(ct) {
      res <- results_by_celltype[[ct]]
      if(is.null(res)) return(NA)
      gene_idx <- which(res$gene == gene)
      if(length(gene_idx) == 0) return(NA)
      
      log2fc <- res$log2FoldChange[gene_idx]
      padj <- res$padj[gene_idx]
      
      # If significant (padj < 0.05), format as bold
      if(!is.na(padj) && padj < 0.05) {
        return(paste0("**", round(log2fc, 3), "**"))
      } else {
        return(round(log2fc, 3))
      }
    })
  }
  
  # Sort by total_significant in descending order
  summary_table <- summary_table[order(-summary_table$total_significant, na.last = TRUE), ]
  
  # Save summary table
  summary_file <- file.path(results_base_dir, "DEG_summary_by_celltype.csv")
  write.csv(summary_table, summary_file, row.names = FALSE)
  
  # Print summary
  cat("\nAnalysis Summary:\n")
  cat("Comparison:", comparison_name, "\n")
  cat("Successfully analyzed cell types:", sum(!is.na(summary_table$total_significant)), "\n")
  cat("Failed analyses:", length(failed_analyses), "\n")
  
  if (length(failed_analyses) > 0) {
    cat("\nCell types that failed analysis:\n")
    for (name in names(failed_analyses)) {
      cat(paste0("  ", name, ": ", failed_analyses[[name]], "\n"))
    }
  }
  
  # Print results summary
  successful_results <- summary_table[!is.na(summary_table$total_significant), ]
  if(nrow(successful_results) > 0) {
    cat("\nTop 5 cell types by number of significant DEGs:\n")
    top_results <- successful_results[order(-successful_results$total_significant), ][1:min(5, nrow(successful_results)), ]
    print(top_results[, c("cell_type", "upregulated", "downregulated", "total_significant")])
  }
  
  # Return results based on user preference
  if (return_results) {
    return(list(
      results_by_celltype = results_by_celltype,
      summary_table = summary_table,
      failed_analyses = failed_analyses,
      comparison_name = comparison_name
    ))
  } else {
    cat("\nAnalysis complete. All results saved to:", results_base_dir, "\n")
    return(invisible(NULL))
  }
}