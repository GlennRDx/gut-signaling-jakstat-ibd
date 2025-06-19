# Pathway-based hexbin plot script for RDS datasets
# Output: Hexbin plots colored by average pathway expression
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(org.Hs.eg.db)
library(KEGGREST)
library(AnnotationDbi)
library(viridis)  # Added for magma color palette

# Paths
data_path <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/'
output_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/plots/pathway_analysis/'

# Function to get gene list from pathway
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

# Function to parse dataset name and get location mapping
parse_dataset_info <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_")[[1]]
  id <- parts[1]
  location_code <- parts[2]
  cellgroup <- parts[3]
  
  # Map location codes to full names
  location_full <- ifelse(location_code == "CO", "Colon", "Ileum")
  
  return(list(id = id, location_code = location_code, location_full = location_full, cellgroup = cellgroup))
}

# Function to calculate average pathway expression
calculate_pathway_expression <- function(seurat_obj, pathway_genes) {
  # Get genes that are present in the Seurat object
  available_genes <- intersect(pathway_genes$SYMBOL, rownames(seurat_obj))
  
  if(length(available_genes) == 0) {
    warning("No pathway genes found in the dataset")
    return(rep(0, ncol(seurat_obj)))
  }
  
  cat("Found", length(available_genes), "out of", nrow(pathway_genes), "pathway genes in dataset\n")
  cat("Pathway genes found:", paste(available_genes, collapse = ", "), "\n")
  
  # Extract expression data for pathway genes
  tryCatch({
    expr_data <- FetchData(seurat_obj, vars = available_genes, layer = "data")
  }, error = function(e1) {
    tryCatch({
      if ("data" %in% Layers(seurat_obj[["RNA"]])) {
        expr_data <- t(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[available_genes, , drop = FALSE])
      } else {
        expr_data <- t(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[available_genes, , drop = FALSE])
      }
    }, error = function(e2) {
      expr_data <- t(GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE])
    })
  })
  
  # Calculate average expression across pathway genes for each cell
  if(length(available_genes) == 1) {
    avg_expression <- as.numeric(expr_data)
  } else {
    avg_expression <- rowMeans(expr_data, na.rm = TRUE)
  }
  
  return(avg_expression)
}

# Global variable to store user's overwrite choice
overwrite_choice <- NULL

# Function to check if plot files exist and handle overwrite decision
check_and_handle_existing_plots <- function(plot_output_dir, dataset_name, pathway_id) {
  # Define expected output files
  hex_filename <- file.path(plot_output_dir, paste0(dataset_name, "_", pathway_id, "_pathway_hexbin.png"))
  
  # Check if plot file already exists
  if (file.exists(hex_filename)) {
    # If overwrite choice hasn't been made yet, prompt user
    if (is.null(overwrite_choice)) {
      cat("Plot results already exist. Overwrite existing plots? (Y/N): ")
      user_input <- toupper(trimws(readline()))
      if (user_input %in% c("Y", "YES")) {
        overwrite_choice <<- TRUE
      } else if (user_input %in% c("N", "NO")) {
        overwrite_choice <<- FALSE
      } else {
        cat("Invalid input. Defaulting to No (N).\n")
        overwrite_choice <<- FALSE
      }
    }
    
    # If user chose not to overwrite, skip this dataset
    if (!overwrite_choice) {
      cat("Skipping (plot exists):", dataset_name, "\n")
      return(FALSE)  # Skip processing
    }
  }
  
  return(TRUE)  # Proceed with processing
}

# Main processing function
process_pathway_hexbin <- function(pathway_id, pathway_type = "KEGG") {
  # Get pathway genes
  cat("Retrieving genes for pathway:", pathway_id, "\n")
  pathway_genes <- get_gene_list(pathway_id, pathway_type)
  
  if(nrow(pathway_genes) == 0) {
    stop("No genes found for pathway: ", pathway_id)
  }
  
  cat("Found", nrow(pathway_genes), "genes in pathway", pathway_id, "\n")
  
  # Get RDS files
  rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
  dataset_names <- tools::file_path_sans_ext(basename(rds_files))
  
  cat("Starting pathway hexbin plot generation for", length(dataset_names), "datasets...\n")
  success <- 0
  total <- 0
  
  seurat_obj <- NULL
  
  for(i in seq_along(rds_files)) {
    total <- total + 1
    dataset_name <- dataset_names[i]
    
    # Parse dataset information
    dataset_info <- parse_dataset_info(dataset_name)
    
    # Create output directory structure for plots
    plot_output_dir <- file.path(output_base, pathway_type, pathway_id, dataset_info$location_full, dataset_info$cellgroup)
    dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Check if plots exist and handle overwrite decision
    if (!check_and_handle_existing_plots(plot_output_dir, dataset_name, pathway_id)) {
      next  # Skip this dataset
    }
    
    cat("Processing:", dataset_name, "\n")
    
    # Clean up previous object
    if (!is.null(seurat_obj)) { rm(seurat_obj); gc(verbose = FALSE) }
    
    tryCatch({
      seurat_obj <- readRDS(rds_files[i])
      
      # Set working directory
      setwd(plot_output_dir)
      
      # Calculate pathway expression
      pathway_expression <- calculate_pathway_expression(seurat_obj, pathway_genes)
      
      # Prepare plot data
      seurat_obj$CD_Status <- factor(seurat_obj$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
      umap_coords <- Embeddings(seurat_obj, reduction = "umap")
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
      metadata <- seurat_obj@meta.data
      plot_data <- cbind(umap_coords, metadata)
      plot_data$pathway_expression <- pathway_expression
      
      # Create UMAP plot showing cell types
      umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Celltypes)) +
        geom_point(size = 0.1, alpha = 0.8) +
        facet_wrap(~ CD_Status, nrow = 1) +
        theme_void() +
        theme(strip.text = element_text(size = 10), 
              legend.position = "right",
              legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6),
              plot.margin = margin(2,2,2,2, "pt"), 
              panel.spacing = unit(0.1, "lines")) +
        coord_fixed() +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 4, shape = 15))) +
        ggtitle("Cell Types")
      
      # Create hexbin plot for pathway expression with magma color gradient
      max_pathway_expr <- max(plot_data$pathway_expression, na.rm = TRUE)
      
      pathway_hex_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, z = pathway_expression)) +
        stat_summary_hex(bins = 50, fun = mean, alpha = 0.8) +
        scale_fill_viridis_c(option = "rocket", 
                             name = "Avg Pathway\nExpression", 
                             limits = c(0, max_pathway_expr)) +
        facet_wrap(~ CD_Status, nrow = 1) +
        theme_void() +
        theme(strip.text = element_text(size = 10), 
              legend.position = "right",
              legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6),
              plot.margin = margin(2,2,2,2, "pt"), 
              panel.spacing = unit(0.1, "lines")) +
        coord_fixed() +
        ggtitle(paste("Pathway", pathway_id, "Expression"))
      
      # Combine plots
      combined_plot <- wrap_plots(list(umap_plot, pathway_hex_plot), ncol = 1, heights = c(1, 1))
      
      # Add overall title
      final_plot <- combined_plot + 
        plot_annotation(
          title = paste("Pathway Analysis:", pathway_id, "-", dataset_name),
          subtitle = paste("Dataset:", dataset_info$location_full, dataset_info$cellgroup),
          theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                        plot.subtitle = element_text(size = 12, hjust = 0.5))
        )
      
      # Save plot
      hex_filename <- paste0(dataset_name, "_", pathway_id, "_pathway_hexbin.png")
      ggsave(hex_filename, plot = final_plot, width = 14, height = 8, units = "in", dpi = 300, bg = "white")
      
      # Clean up
      rm(umap_plot, pathway_hex_plot, combined_plot, final_plot, plot_data)
      
      success <- success + 1
      cat("Completed:", dataset_name, "\n")
      
    }, error = function(e) {
      cat("Error processing", dataset_name, ":", e$message, "\n")
    })
    
    gc(verbose = FALSE)
  }
  
  # Clean up
  if (!is.null(seurat_obj)) { rm(seurat_obj) }
  gc(verbose = FALSE)
  
  cat("Pathway hexbin plot generation completed:", success, "of", total, "datasets successful\n")
}

# Execute the analysis for the specified pathway
# You can change the pathway_id and pathway_type as needed
pathway_id <- "hsa04136"
pathway_type <- "KEGG"

cat("Starting pathway analysis for:", pathway_id, "\n")
process_pathway_hexbin(pathway_id, pathway_type)