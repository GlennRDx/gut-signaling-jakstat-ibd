# SCRIPT 2 - IMPROVED FOR COMPREHENSIVE DGE ANALYSIS

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
sobj <- readRDS('/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/04_li/compatible.rds')
sobj <- subset(sobj, subset = Tissue == "Terminal Ileum")

# Create a new column in your Seurat object metadata called "CD_Status"
sobj$CD_Status <- "Unknown"  # Initialize with a default value

# Map the histology classifications to the new groups
sobj$CD_Status[sobj$Histology == "Control"] <- "Control"
# Active CD: mild, moderate, and severe
sobj$CD_Status[sobj$Histology %in% c("Mild", "Moderate", "Severe")] <- "Active_CD"
# Inactive CD: normal and quiescent
sobj$CD_Status[sobj$Histology %in% c("Normal_CD", "Quiescent")] <- "Inactive_CD"

# Verify the new annotation
table(sobj$CD_Status)

# Pseudobulk the counts based on condition-donor-celltype
sobj <- AggregateExpression(sobj, assays = "RNA", return.seurat = T, 
                            group.by = c("CD_Status", "PatientID", "Celltypes"))

# Extract pseudobulk count matrix and metadata
bulk_counts <- sobj@assays$RNA@layers$counts
rownames(bulk_counts) <- rownames(sobj@assays$RNA)
bulk_metadata <- sobj@meta.data

colnames(bulk_counts) <- colnames(sobj)

# Add identifiers to metadata for easier filtering
bulk_metadata$sample_id <- rownames(bulk_metadata)
bulk_metadata$CD_Status <- sapply(strsplit(bulk_metadata$sample_id, "_"), `[`, 1)
bulk_metadata$PatientID <- sapply(strsplit(bulk_metadata$sample_id, "_"), `[`, 2)
bulk_metadata$Celltypes <- sapply(strsplit(bulk_metadata$sample_id, "_"), `[`, 3)

# Use the bulk_counts and bulk_metadata for further analysis
counts <- bulk_counts
meta_data <- data.frame(
  row.names = rownames(bulk_metadata),
  condition = bulk_metadata$CD_Status,
  sample = bulk_metadata$PatientID,
  celltype = bulk_metadata$Celltypes
)

# Get all unique cell types
cell_types <- unique(meta_data$celltype)

# Define results directory path
results_base_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DEG_Results_Li_Dataset"

# Create results directories
dir.create(results_base_dir, showWarnings = FALSE, recursive = TRUE)
comparisons <- c("Active_CD_vs_Control", "Inactive_CD_vs_Control")
for (comparison in comparisons) {
  dir.create(paste0(results_base_dir, "/", comparison), showWarnings = FALSE)
  for (cell_type in cell_types) {
    safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
    dir.create(paste0(results_base_dir, "/", comparison, "/", safe_cell_type), showWarnings = FALSE, recursive = TRUE)
  }
}

# IMPROVED FUNCTIONS WITH BETTER ERROR HANDLING

# Create PCA plot with error handling for low samples/genes
create_pca_plot <- function(dds, cell_type, comparison, output_dir) {
  tryCatch({
    # Check if we have enough samples for meaningful PCA
    if (ncol(dds) < 3) {
      cat("      Skipping PCA: too few samples (", ncol(dds), ")\n")
      return(NULL)
    }
    
    vst_data <- tryCatch({
      vst(dds, blind = FALSE)
    }, error = function(e) {
      # If vst fails, try varianceStabilizingTransformation directly
      varianceStabilizingTransformation(dds, blind = FALSE)
    })
    
    pca_data <- plotPCA(vst_data, intgroup = c("condition"), returnData = TRUE)
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
  }, error = function(e) {
    cat("      Error creating PCA plot:", e$message, "\n")
    return(NULL)
  })
}

# Create Dispersion plot with error handling
create_dispersion_plot <- function(dds, cell_type, comparison, output_dir) {
  tryCatch({
    disp_file <- paste0(output_dir, "/Dispersion_plot.pdf")
    pdf(disp_file, width = 10, height = 8)
    plotDispEsts(dds, main = paste0("Dispersion Plot: ", cell_type, " - ", comparison))
    dev.off()
  }, error = function(e) {
    cat("      Error creating dispersion plot:", e$message, "\n")
    if (dev.cur() > 1) dev.off()  # Close PDF device if open
  })
}

# Create MA plot with error handling
create_ma_plot <- function(res, cell_type, comparison, output_dir) {
  tryCatch({
    ma_file <- paste0(output_dir, "/MA_plot.pdf")
    pdf(ma_file, width = 10, height = 8)
    plotMA(res, main = paste0("MA Plot: ", cell_type, " - ", comparison), ylim = c(-5, 5), cex = 0.8)
    abline(h = 0, col = "gray", lwd = 2)
    dev.off()
  }, error = function(e) {
    cat("      Error creating MA plot:", e$message, "\n")
    if (dev.cur() > 1) dev.off()  # Close PDF device if open
  })
}

# IMPROVED DESeq2 function with better filtering and error handling
perform_DESeq2 <- function(counts_matrix, metadata, condition1, condition2, cell_type, comparison_type = "") {
  # Define comparison name
  if (condition2 == "Active_CD" || condition2 == "Active-CD") {
    comparison <- "Active_CD_vs_Control"
  } else if (condition2 == "Inactive_CD" || condition2 == "Inactive-CD") {
    comparison <- "Inactive_CD_vs_Control"
  } else {
    comparison <- paste0(gsub("-", "_", condition2), "_vs_", condition1)
  }
  
  # Sanitize cell type name for directory creation
  safe_cell_type <- gsub("[/\\\\?%*:|\"<>]", "_", gsub(" ", "_", cell_type))
  
  # Create output directory
  output_dir <- paste0(results_base_dir, "/", comparison, "/", safe_cell_type)
  
  # Subset data for the two conditions to compare
  keep <- metadata$condition %in% c(condition1, condition2)
  counts_subset <- counts_matrix[, keep]
  meta_subset <- metadata[keep, ]
  
  cat("    Found", sum(keep), "samples for comparison\n")
  cat("    Condition counts:", table(meta_subset$condition), "\n")
  
  # Enhanced sample count checking
  condition_counts <- table(meta_subset$condition)
  min_samples_per_group <- 3  # Minimum samples per group for DESeq2
  
  if(any(condition_counts < min_samples_per_group)) {
    error_msg <- paste0("Insufficient samples per group (need >=", min_samples_per_group, " per group)")
    cat("   ", error_msg, "\n")
    return(list(error = error_msg))
  }
  
  # Set reference level
  meta_subset$condition <- factor(meta_subset$condition, levels = c(condition1, condition2))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = meta_subset,
    design = ~ condition
  )
  
  # IMPROVED gene filtering strategy
  # Keep genes with at least 5 counts in at least 2 samples
  keep_genes <- rowSums(counts(dds) >= 5) >= 2
  dds_filtered <- dds[keep_genes, ]
  
  cat("    Genes before filtering:", nrow(dds), "\n")
  cat("    Genes after filtering:", nrow(dds_filtered), "\n")
  
  # Check if we still have enough genes after filtering
  if(nrow(dds_filtered) < 100) {  # Need reasonable number of genes
    # Try more permissive filtering
    cat("    Trying more permissive filtering...\n")
    keep_genes_permissive <- rowSums(counts(dds) >= 1) >= 2
    dds_filtered <- dds[keep_genes_permissive, ]
    cat("    Genes after permissive filtering:", nrow(dds_filtered), "\n")
    
    if(nrow(dds_filtered) < 50) {
      error_msg <- paste0("Too few genes after filtering (", nrow(dds_filtered), " genes, need >=50)")
      cat("   ", error_msg, "\n")
      return(list(error = error_msg))
    }
  }
  
  dds <- dds_filtered
  
  # Run DESeq2 with error handling
  cat("    Running DESeq2...\n")
  tryCatch({
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", condition2, condition1))
    
    # Create diagnostic plots (with error handling built into functions)
    cat("    Creating plots...\n")
    create_pca_plot(dds, cell_type, comparison, output_dir)
    create_dispersion_plot(dds, cell_type, comparison, output_dir)
    create_ma_plot(res, cell_type, comparison, output_dir)
    
    # Create volcano plot with error handling
    tryCatch({
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      volcano_file <- paste0(output_dir, "/Volcano_plot.pdf")
      pdf(volcano_file, width = 10, height = 8)
      print(EnhancedVolcano(res_df,
                            lab = res_df$gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            title = paste0(cell_type, ': ', condition2, ' vs ', condition1),
                            subtitle = 'Ascending Colon - Li Dataset',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            pointSize = 3.0,
                            labSize = 4.0))
      dev.off()
    }, error = function(e) {
      cat("      Error creating volcano plot:", e$message, "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    # Save results table
    res_df <- as.data.frame(res)
    output_file <- paste0(output_dir, "/DEG_results.csv")
    write.csv(res_df, file = output_file)
    
    # Calculate summary statistics
    sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
    sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
    
    cat("    Analysis completed successfully!\n")
    cat("    Total genes analyzed:", nrow(res_df), "\n")
    cat("    Significant upregulated:", sig_up, "\n")
    cat("    Significant downregulated:", sig_down, "\n")
    
    return(res_df)
    
  }, error = function(e) {
    error_msg <- paste0("DESeq2 analysis failed: ", e$message)
    cat("   ", error_msg, "\n")
    return(list(error = error_msg))
  })
}

# Initialize results storage and logging
results_Active_CD_vs_Control_by_celltype <- list()
results_Inactive_CD_vs_Control_by_celltype <- list()
failed_analyses <- list()

# Initialize log vectors
active_cd_log <- c()
inactive_cd_log <- c()

# Print initial data summary
cat("=== INITIAL DATA SUMMARY ===\n")
cat("Cell types found:", length(cell_types), "\n")
cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Total samples:", nrow(meta_data), "\n")
cat("Condition distribution across all samples:\n")
print(table(meta_data$condition))
cat("\n")

# Analyze sample distribution by cell type
cat("=== SAMPLE DISTRIBUTION BY CELL TYPE ===\n")
for (cell_type in cell_types) {
  cell_type_indices <- meta_data$celltype == cell_type
  meta_cell_type <- meta_data[cell_type_indices, ]
  condition_table <- table(meta_cell_type$condition)
  cat(sprintf("%-20s: %s\n", cell_type, paste(names(condition_table), "=", condition_table, collapse = ", ")))
}
cat("\n")

cat("=== STARTING DGE ANALYSIS ===\n")

# FIXED: Function to handle single comparison analysis
run_comparison <- function(cell_type, counts_cell_type, meta_cell_type, control_condition, test_condition, comparison_name) {
  result <- tryCatch({
    available_conditions <- unique(meta_cell_type$condition)
    cat("    Available conditions:", paste(available_conditions, collapse = ", "), "\n")
    
    if(control_condition %in% available_conditions && test_condition %in% available_conditions) {
      # Check sample counts
      control_samples <- sum(meta_cell_type$condition == control_condition)
      test_samples <- sum(meta_cell_type$condition == test_condition)
      
      cat("    Sample counts:", control_condition, "=", control_samples, ",", test_condition, "=", test_samples, "\n")
      
      result <- perform_DESeq2(counts_cell_type, meta_cell_type, control_condition, test_condition, cell_type, comparison_name)
      
      # Check if result is an error
      if(is.list(result) && !is.null(result$error)) {
        # Return error information instead of using return()
        list(success = FALSE, error = result$error)
      } else {
        list(success = TRUE, result = result)
      }
    } else {
      missing_conditions <- c()
      if(!control_condition %in% available_conditions) missing_conditions <- c(missing_conditions, control_condition)
      if(!test_condition %in% available_conditions) missing_conditions <- c(missing_conditions, test_condition)
      
      error_msg <- paste0("Missing conditions - ", paste(missing_conditions, collapse = ", "))
      cat("   ", error_msg, "\n")
      list(success = FALSE, error = error_msg)
    }
  }, error = function(e) {
    error_msg <- paste0("Unexpected error: ", e$message)
    cat("   ", error_msg, "\n")
    list(success = FALSE, error = error_msg)
  })
  
  return(result)
}

for (cell_type in cell_types) {
  cat("Processing cell type:", cell_type, "\n")
  
  # Subset data for this cell type
  cell_type_indices <- meta_data$celltype == cell_type
  counts_cell_type <- counts[, cell_type_indices]
  meta_cell_type <- meta_data[cell_type_indices, ]
  
  cat("  Total samples for", cell_type, ":", nrow(meta_cell_type), "\n")
  
  # Compare Active_CD vs Control
  cat("  Running Active_CD vs Control comparison...\n")
  # Determine the correct condition name
  available_conditions <- unique(meta_cell_type$condition)
  active_condition <- ifelse("Active-CD" %in% available_conditions, "Active-CD", "Active_CD")
  
  active_result <- run_comparison(cell_type, counts_cell_type, meta_cell_type, "Control", active_condition, "Active_CD_vs_Control")
  
  if(active_result$success) {
    results_Active_CD_vs_Control_by_celltype[[cell_type]] <- active_result$result
  } else {
    results_Active_CD_vs_Control_by_celltype[[cell_type]] <- NULL
    active_cd_log <- c(active_cd_log, paste0(cell_type, ": ", active_result$error))
  }
  
  # Compare Inactive_CD vs Control
  cat("  Running Inactive_CD vs Control comparison...\n")
  # Determine the correct condition name
  inactive_condition <- ifelse("Inactive-CD" %in% available_conditions, "Inactive-CD", "Inactive_CD")
  
  inactive_result <- run_comparison(cell_type, counts_cell_type, meta_cell_type, "Control", inactive_condition, "Inactive_CD_vs_Control")
  
  if(inactive_result$success) {
    results_Inactive_CD_vs_Control_by_celltype[[cell_type]] <- inactive_result$result
  } else {
    results_Inactive_CD_vs_Control_by_celltype[[cell_type]] <- NULL
    inactive_cd_log <- c(inactive_cd_log, paste0(cell_type, ": ", inactive_result$error))
  }
  
  cat("  Completed", cell_type, "\n\n")
}

# Write detailed log files
active_cd_log_file <- paste0(results_base_dir, "/Active_CD_vs_Control/analysis_log.txt")
if(length(active_cd_log) > 0) {
  log_content <- c(
    "=== ACTIVE CD vs CONTROL ANALYSIS LOG ===",
    paste("Generated on:", Sys.time()),
    "",
    "Failed analyses:",
    active_cd_log,
    "",
    paste("Total failed:", length(active_cd_log)),
    paste("Total attempted:", length(cell_types))
  )
  writeLines(log_content, active_cd_log_file)
  cat("Active_CD vs Control log saved to:", active_cd_log_file, "\n")
} else {
  writeLines("All cell types successfully analyzed for Active_CD vs Control comparison", active_cd_log_file)
}

inactive_cd_log_file <- paste0(results_base_dir, "/Inactive_CD_vs_Control/analysis_log.txt")
if(length(inactive_cd_log) > 0) {
  log_content <- c(
    "=== INACTIVE CD vs CONTROL ANALYSIS LOG ===",
    paste("Generated on:", Sys.time()),
    "",
    "Failed analyses:",
    inactive_cd_log,
    "",
    paste("Total failed:", length(inactive_cd_log)),
    paste("Total attempted:", length(cell_types))
  )
  writeLines(log_content, inactive_cd_log_file)
  cat("Inactive_CD vs Control log saved to:", inactive_cd_log_file, "\n")
} else {
  writeLines("All cell types successfully analyzed for Inactive_CD vs Control comparison", inactive_cd_log_file)
}

# Create enhanced summary table WITHOUT status columns
cat("Creating summary table...\n")
summary_table <- data.frame(
  cell_type = cell_types,
  Active_CD_vs_Control_up = sapply(cell_types, function(ct) {
    res <- results_Active_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res) || (is.list(res) && !is.null(res$error))) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  }),
  Active_CD_vs_Control_down = sapply(cell_types, function(ct) {
    res <- results_Active_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res) || (is.list(res) && !is.null(res$error))) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  }),
  Inactive_CD_vs_Control_up = sapply(cell_types, function(ct) {
    res <- results_Inactive_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res) || (is.list(res) && !is.null(res$error))) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  }),
  Inactive_CD_vs_Control_down = sapply(cell_types, function(ct) {
    res <- results_Inactive_CD_vs_Control_by_celltype[[ct]]
    if(is.null(res) || (is.list(res) && !is.null(res$error))) return(NA)
    sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  }),
  Total_samples_available = sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    sum(cell_type_indices)
  }),
  Control_samples = sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    meta_cell_type <- meta_data[cell_type_indices, ]
    sum(meta_cell_type$condition == "Control")
  }),
  Active_CD_samples = sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    meta_cell_type <- meta_data[cell_type_indices, ]
    sum(meta_cell_type$condition %in% c("Active_CD", "Active-CD"))
  }),
  Inactive_CD_samples = sapply(cell_types, function(ct) {
    cell_type_indices <- meta_data$celltype == ct
    meta_cell_type <- meta_data[cell_type_indices, ]
    sum(meta_cell_type$condition %in% c("Inactive_CD", "Inactive-CD"))
  })
)

# Save summary table
write.csv(summary_table, paste0(results_base_dir, "/DEG_summary_by_celltype.csv"), row.names = FALSE)

# Print final summary
cat("\n=== FINAL ANALYSIS SUMMARY ===\n")
cat("Results saved to:", results_base_dir, "\n")
cat("Total cell types analyzed:", length(cell_types), "\n")
cat("Comparisons performed: Active_CD vs Control, Inactive_CD vs Control\n\n")

# Summary statistics (calculated from non-NA values)
successful_active <- sum(!is.na(summary_table$Active_CD_vs_Control_up))
successful_inactive <- sum(!is.na(summary_table$Inactive_CD_vs_Control_up))

cat("Analysis Success Rate:\n")
cat("  Active_CD vs Control:", successful_active, "/", length(cell_types), 
    "(", round(100*successful_active/length(cell_types), 1), "%)\n")
cat("  Inactive_CD vs Control:", successful_inactive, "/", length(cell_types), 
    "(", round(100*successful_inactive/length(cell_types), 1), "%)\n\n")

print("Detailed summary by cell type:")
print(summary_table)

# Clean up
rm(sobj, counts)
gc()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Check the results directory for all outputs and detailed logs.\n")
cat("Summary table saved as: DEG_summary_by_celltype.csv\n")