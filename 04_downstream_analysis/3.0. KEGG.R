# Automated pathway heatmap generation for all datasets - Multiple KEGG pathways
library(pheatmap)
library(KEGGREST)

# Function to get KEGG pathway name
get_kegg_pathway_name <- function(pid) {
  pathway_info <- keggGet(pid)
  pathway_name <- pathway_info[[1]]$NAME
  pathway_name <- gsub(" - Mus musculus \\(house mouse\\)", "", pathway_name)
  return(pathway_name)
}

# Find all genes associated with a KEGG pathway
get_gene_list <- function(identifier, type) {
  library(org.Hs.eg.db)
  library(KEGGREST)
  library(AnnotationDbi)
  if (type == "GO") {
    genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keytype = "GOALL", 
                                   keys = identifier, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else if (type == "KEGG") {
    kegg_genes <- gsub("hsa:", "", keggLink("hsa", identifier))
    genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keytype = "ENTREZID", 
                                   keys = kegg_genes, 
                                   columns = c("ENSEMBL", "SYMBOL"))
  } else {
    stop("Invalid type. Use 'GO' or 'KEGG'.")
  }
  return(genes)
}

# Simplified pathway heatmap function
pathway_heatmap <- function(df_list, 
                            pid = NULL, 
                            scale_to_one = FALSE, 
                            remove_na_rows = FALSE, 
                            order_by_sum = TRUE, 
                            custom_gene_list = NULL,
                            output_dir = ".",
                            file_name = "heatmap.png",
                            fontsize = 10,
                            w = 1500,
                            h = 2200,
                            heatmap_title = "Heatmap") {
  
  # Get gene list
  if (!is.null(custom_gene_list)) {
    gene_list <- unique(custom_gene_list)
  } else if (startsWith(pid, "hsa")) {
    gene_data <- get_gene_list(pid, 'KEGG')
    gene_list <- unique(gene_data$SYMBOL[!is.na(gene_data$SYMBOL)])
  } else if (startsWith(pid, "GO")) {
    gene_data <- get_gene_list(pid, 'GO')
    gene_list <- unique(gene_data$SYMBOL[!is.na(gene_data$SYMBOL)])
  }
  
  # Function to convert cell type names to proper capitalization
  format_cell_type_name <- function(name) {
    # Convert to proper case (first letter capitalized, rest lowercase)
    paste0(toupper(substr(name, 1, 1)), tolower(substr(name, 2, nchar(name))))
  }
  
  # Initialize matrices
  p_val <- matrix(NA, nrow = length(gene_list), ncol = length(df_list))
  logFC <- matrix(NA, nrow = length(gene_list), ncol = length(df_list))
  rownames(p_val) <- rownames(logFC) <- gene_list
  
  # Use full cell type names with proper capitalization instead of last 3 letters
  colnames(p_val) <- colnames(logFC) <- sapply(names(df_list), format_cell_type_name)
  
  # Populate matrices
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    for (j in seq_along(gene_list)) {
      gene <- gene_list[j]
      if (gene %in% df$gene) {
        gene_idx <- which(df$gene == gene)[1]
        p_val[j, i] <- df$padj[gene_idx]
        if (!is.na(df$padj[gene_idx]) && df$padj[gene_idx] < 0.05) {
          logFC[j, i] <- df$log2FoldChange[gene_idx]
        }
      }
    }
  }
  
  # Convert to data frames
  p_val <- as.data.frame(p_val)
  logFC <- as.data.frame(logFC)
  
  # Handle NA rows
  if (remove_na_rows) {
    non_na_rows <- apply(logFC, 1, function(row) !all(is.na(row)))
    p_val <- p_val[non_na_rows, , drop = FALSE]
    logFC <- logFC[non_na_rows, , drop = FALSE]
  } else {
    na_rows <- apply(p_val, 1, function(row) all(is.na(row)))
    p_val <- rbind(p_val[!na_rows, , drop = FALSE], p_val[na_rows, , drop = FALSE])
    logFC <- logFC[rownames(p_val), , drop = FALSE]
  }
  
  # Clean data
  p_val_clean <- p_val
  p_val_clean[is.na(p_val_clean)] <- 1
  logFC_clean <- logFC
  logFC_clean[is.na(logFC_clean)] <- 0
  
  # Order by sum if requested
  if (order_by_sum) {
    significant_mask <- as.matrix(p_val_clean) < 0.05
    significant_mask[is.na(significant_mask)] <- FALSE
    significant_logFC <- logFC_clean
    significant_logFC[!significant_mask] <- 0
    logFC_sums <- rowSums(as.matrix(significant_logFC), na.rm = TRUE)
    order_idx <- order(logFC_sums, decreasing = TRUE)
    logFC_clean <- logFC_clean[order_idx, , drop = FALSE]
    p_val_clean <- p_val_clean[order_idx, , drop = FALSE]
  }
  
  # Prepare for plotting
  logFC_matrix <- as.matrix(logFC_clean)
  p_val_matrix <- as.matrix(p_val_clean)
  
  # Handle edge cases for matrix values
  logFC_matrix[is.na(logFC_matrix)] <- 0
  logFC_matrix[is.infinite(logFC_matrix)] <- 0
  
  # Check if we have any valid data to plot
  if (all(logFC_matrix == 0) || nrow(logFC_matrix) == 0) {
    warning("No valid data to plot for this pathway")
    return(NULL)
  }
  
  max_abs_logFC <- max(abs(logFC_matrix), na.rm = TRUE)
  
  # Handle case where max_abs_logFC is 0 or very small
  if (max_abs_logFC == 0 || is.na(max_abs_logFC)) {
    max_abs_logFC <- 1
  }
  
  color_limits <- if (scale_to_one && max_abs_logFC <= 1) c(-1, 1) else c(-max_abs_logFC, max_abs_logFC)
  breaks <- seq(color_limits[1], color_limits[2], length.out = 101)
  color_palette <- colorRampPalette(c("red", "white", "green"))(100)
  
  # Create significance symbols
  pval_to_significance <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.0001) return("****")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("")
  }
  significance_symbols <- apply(p_val_matrix, c(1, 2), pval_to_significance)
  
  # Save heatmap
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  file_path <- file.path(output_dir, file_name)
  
  # Additional safety check before plotting
  if (any(is.na(breaks)) || any(is.infinite(breaks))) {
    warning("Invalid breaks detected, skipping heatmap generation")
    return(NULL)
  }
  
  png(filename = file_path, width = w, height = h, res = 200)
  tryCatch({
    pheatmap(logFC_matrix, 
             main = heatmap_title, 
             cluster_rows = !order_by_sum,
             cluster_cols = T, 
             display_numbers = significance_symbols, 
             color = color_palette, 
             breaks = breaks,
             border_color = NA,
             fontsize = fontsize,
             labels_col = colnames(logFC_matrix),
             angle_col = 45)  # Add angled column labels
  }, error = function(e) {
    dev.off()  # Make sure to close the PNG device even if pheatmap fails
    stop(paste("Error in pheatmap:", e$message))
  })
  dev.off()
  
  message("Heatmap saved to ", file_path)
}

# Load DEG results
load_deg_results <- function(cell_type_dir) {
  deg_file <- file.path(cell_type_dir, "DEG_results.csv")
  if (file.exists(deg_file)) {
    return(read.csv(deg_file, stringsAsFactors = FALSE))
  }
  return(NULL)
}

# Modified function to process all datasets for multiple KEGG pathways
process_all_datasets <- function(base_dir, pathway_ids = c('hsa04630')) {
  
  comparison_dirs <- c("ActiveCD_vs_Control", "InactiveCD_vs_Control")
  
  # Loop through each pathway ID
  for (pid in pathway_ids) {
    message(paste("Processing pathway:", pid))
    
    # Get the actual KEGG pathway name
    tryCatch({
      # Suppress the select() warning messages
      suppressMessages({
        pathway_name <- get_kegg_pathway_name(pid)
      })
      # Clean the pathway name for use in filename (remove special characters)
      pathway_name_clean <- gsub("[^A-Za-z0-9_-]", "_", pathway_name)
      pathway_name_clean <- gsub("_+", "_", pathway_name_clean)  # Replace multiple underscores with single
      pathway_name_clean <- gsub("^_|_$", "", pathway_name_clean)  # Remove leading/trailing underscores
    }, error = function(e) {
      message(paste("Error retrieving pathway name for", pid, "- using pathway ID as name"))
      pathway_name <- pid
      pathway_name_clean <- pid
    })
    
    for (comparison in comparison_dirs) {
      comparison_path <- file.path(base_dir, comparison)
      if (!dir.exists(comparison_path)) next
      
      dataset_dirs <- list.dirs(comparison_path, recursive = FALSE, full.names = FALSE)
      dataset_dirs <- dataset_dirs[grepl("^[0-9]+_", dataset_dirs)]
      
      for (dataset in dataset_dirs) {
        dataset_path <- file.path(comparison_path, dataset)
        cell_type_dirs <- list.dirs(dataset_path, recursive = FALSE, full.names = FALSE)
        
        df_list <- list()
        for (cell_type in cell_type_dirs) {
          cell_type_path <- file.path(dataset_path, cell_type)
          deg_data <- load_deg_results(cell_type_path)
          if (!is.null(deg_data)) df_list[[cell_type]] <- deg_data
        }
        
        if (length(df_list) > 0) {
          tryCatch({
            # Suppress select() warning messages from get_gene_list
            suppressMessages({
              result <- pathway_heatmap(df_list, 
                                        pid = pid, 
                                        scale_to_one = TRUE, 
                                        remove_na_rows = TRUE, 
                                        order_by_sum = TRUE,
                                        output_dir = dataset_path,
                                        file_name = paste0(pathway_name_clean, ".png"),
                                        heatmap_title = paste(pathway_name, "-", dataset, "-", comparison),
                                        fontsize = 12,
                                        w = 3000,
                                        h = 2600)
            })
          }, error = function(e) {
            message(paste("Error processing pathway", pid, "for dataset", dataset, "comparison", comparison, ":", e$message))
          })
        }
      }
    }
    message(paste("Finished processing pathway:", pid))
  }
  message("All pathways and datasets processed!")
}

# Usage examples:

# Option 1: Define your kegg_list and use it directly
kegg_list <- c('hsa04630', 'hsa04136', 'hsa04010', 'hsa04064', 'hsa04530')
base_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results"
process_all_datasets(base_dir, pathway_ids = kegg_list)

# Option 2: Call with individual pathways
# process_all_datasets(base_dir, pathway_ids = c('hsa04630', 'hsa04060'))

# Option 3: Process single pathway (backwards compatible)
# process_all_datasets(base_dir, pathway_ids = 'hsa04630')