#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Set the directory paths
data_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/"
output_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/QA"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get all .rds files in the directory
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

# Check if any files were found
if (length(rds_files) == 0) {
  stop("No .rds files found in the specified directory")
}

# Initialize results dataframe
results_summary <- data.frame(
  Dataset = character(),
  Total_Cells = integer(),
  Mean_Depth = numeric(),
  Median_Depth = numeric(),
  Min_Depth = numeric(),
  Max_Depth = numeric(),
  SD_Depth = numeric(),
  stringsAsFactors = FALSE
)

# Initialize list to store individual cell data for plotting
all_cell_data <- list()

# First pass: collect all sequencing depth data to determine global limits
cat("=== FIRST PASS: Collecting data for global axis limits ===\n")
global_max_depth <- 0
global_max_count <- 0

for (i in seq_along(rds_files)) {
  file_path <- rds_files[i]
  dataset_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat("Scanning dataset:", dataset_name, "for axis limits...\n")
  
  tryCatch({
    seurat_obj <- readRDS(file_path)
    
    if (!inherits(seurat_obj, "Seurat")) {
      rm(seurat_obj)
      gc()
      next
    }
    
    # Get raw counts for axis limit calculation
    if ("RNA" %in% names(seurat_obj@assays)) {
      assay_name <- "RNA"
    } else {
      assay_name <- DefaultAssay(seurat_obj)
    }
    
    if (inherits(seurat_obj@assays[[assay_name]], "Assay5")) {
      raw_counts <- LayerData(seurat_obj, assay = assay_name, layer = "counts")
    } else {
      raw_counts <- seurat_obj@assays[[assay_name]]@counts
    }
    
    sequencing_depth <- colSums(raw_counts)
    
    # Update global maximums
    global_max_depth <- max(global_max_depth, max(sequencing_depth))
    
    # Calculate histogram to get max count for y-axis
    hist_data <- hist(sequencing_depth, breaks = 50, plot = FALSE)
    global_max_count <- max(global_max_count, max(hist_data$counts))
    
    # Clean up
    rm(seurat_obj, raw_counts, sequencing_depth, hist_data)
    gc()
    
  }, error = function(e) {
    cat("Error scanning", dataset_name, ":", conditionMessage(e), "\n")
    if (exists("seurat_obj")) rm(seurat_obj)
    if (exists("raw_counts")) rm(raw_counts)
    if (exists("sequencing_depth")) rm(sequencing_depth)
    gc()
  })
}

cat("Global max sequencing depth:", global_max_depth, "\n")
cat("Global max histogram count:", global_max_count, "\n\n")

cat("=== SECOND PASS: Processing datasets with fixed axes ===\n")

# Process each dataset with memory management
for (i in seq_along(rds_files)) {
  file_path <- rds_files[i]
  dataset_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat("Processing dataset:", dataset_name, "\n")
  cat("Memory usage before loading:", round(object.size(ls())/1024/1024, 2), "MB\n")
  
  tryCatch({
    # Load the Seurat object
    seurat_obj <- readRDS(file_path)
    
    # Ensure it's a Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
      warning(paste("File", dataset_name, "is not a Seurat object. Skipping."))
      rm(seurat_obj)  # Clean up
      gc()  # Force garbage collection
      next
    }
    
    # Calculate sequencing depth using RAW counts
    # Handle both Seurat v4 and v5 assay structures
    if ("RNA" %in% names(seurat_obj@assays)) {
      assay_name <- "RNA"
    } else {
      assay_name <- DefaultAssay(seurat_obj)
    }
    
    # Get raw counts - different methods for v4 vs v5
    if (inherits(seurat_obj@assays[[assay_name]], "Assay5")) {
      # Seurat v5 - use LayerData or GetAssayData
      raw_counts <- LayerData(seurat_obj, assay = assay_name, layer = "counts")
    } else {
      # Seurat v4 and earlier - use @counts slot
      raw_counts <- seurat_obj@assays[[assay_name]]@counts
    }
    
    # Calculate sequencing depth per cell
    sequencing_depth <- colSums(raw_counts)
    
    # Calculate summary statistics immediately
    total_cells <- length(sequencing_depth)
    mean_depth <- mean(sequencing_depth)
    median_depth <- median(sequencing_depth)
    min_depth <- min(sequencing_depth)
    max_depth <- max(sequencing_depth)
    sd_depth <- sd(sequencing_depth)
    
    # Add to results summary
    results_summary <- rbind(results_summary, data.frame(
      Dataset = dataset_name,
      Total_Cells = total_cells,
      Mean_Depth = round(mean_depth, 2),
      Median_Depth = round(median_depth, 2),
      Min_Depth = min_depth,
      Max_Depth = max_depth,
      SD_Depth = round(sd_depth, 2)
    ))
    
    # Store only essential data for plotting (not full cell IDs to save memory)
    all_cell_data[[dataset_name]] <- data.frame(
      Dataset = dataset_name,
      Sequencing_Depth = sequencing_depth
    )
    
    # Create individual histogram for this dataset with fixed axes
    p <- ggplot(data.frame(depth = sequencing_depth), aes(x = depth)) +
      geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "black") +
      labs(title = paste("Sequencing Depth Distribution -", dataset_name),
           x = "Total UMIs per cell",
           y = "Number of cells") +
      xlim(0, global_max_depth) +
      ylim(0, global_max_count) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Save individual plot to QA directory
    ggsave(filename = file.path(output_dir, paste0("sequencing_depth_", dataset_name, ".png")),
           plot = p, width = 8, height = 6, dpi = 300)
    
    cat("Completed:", dataset_name, "- Cells:", total_cells, "- Mean depth:", round(mean_depth, 2), "\n")
    
    # MEMORY CLEANUP - Critical for large datasets
    rm(seurat_obj, raw_counts, sequencing_depth, p)
    gc()  # Force garbage collection
    
    cat("Memory usage after cleanup:", round(object.size(ls())/1024/1024, 2), "MB\n\n")
    
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", conditionMessage(e), "\n")
    # Clean up on error
    if (exists("seurat_obj")) rm(seurat_obj)
    if (exists("raw_counts")) rm(raw_counts)
    if (exists("sequencing_depth")) rm(sequencing_depth)
    if (exists("p")) rm(p)
    gc()
    cat("\n")
  })
}

# Print summary table
cat("=== SEQUENCING DEPTH SUMMARY ===\n")
print(results_summary)

# Save summary table to QA directory
write.csv(results_summary, file.path(output_dir, "sequencing_depth_summary.csv"), row.names = FALSE)

# Create combined plot if multiple datasets
if (length(all_cell_data) > 1) {
  combined_data <- do.call(rbind, all_cell_data)
  
  # Combined histogram with fixed axes
  p_combined <- ggplot(combined_data, aes(x = Sequencing_Depth, fill = Dataset)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(title = "Sequencing Depth Distribution - All Datasets",
         x = "Total UMIs per cell",
         y = "Number of cells") +
    xlim(0, global_max_depth) +
    ylim(0, global_max_count) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, "sequencing_depth_combined.png"), plot = p_combined, 
         width = 10, height = 8, dpi = 300)
  
  # Faceted plot with fixed axes
  p_faceted <- ggplot(combined_data, aes(x = Sequencing_Depth)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "black") +
    facet_wrap(~Dataset, scales = "fixed") +  # Changed from "free" to "fixed"
    labs(title = "Sequencing Depth Distribution by Dataset",
         x = "Total UMIs per cell",
         y = "Number of cells") +
    xlim(0, global_max_depth) +
    ylim(0, global_max_count) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir, "sequencing_depth_faceted.png"), plot = p_faceted, 
         width = 12, height = 8, dpi = 300)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_dir, "\n")
cat("- sequencing_depth_summary.csv\n")
cat("- Individual histograms for each dataset\n")
if (length(all_cell_data) > 1) {
  cat("- sequencing_depth_combined.png\n")
  cat("- sequencing_depth_faceted.png\n")
}
cat("- Original .rds files left unchanged\n")