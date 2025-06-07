# Batch processing script for Seurat objects
# Applies normalization, feature selection, scaling, PCA, and UMAP to each RDS file

library(Seurat)

# Define input and output directories
input_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/"
output_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Get list of all RDS files in input directory
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)

if (length(rds_files) == 0) {
  stop("No RDS files found in the input directory!")
}

cat("Found", length(rds_files), "RDS files to process:\n")
cat(paste(rds_files, collapse = "\n"), "\n\n")

# Process each RDS file
for (i in seq_along(rds_files)) {
  file_name <- rds_files[i]
  input_path <- file.path(input_dir, file_name)
  output_path <- file.path(output_dir, file_name)
  
  cat("Processing file", i, "of", length(rds_files), ":", file_name, "\n")
  
  tryCatch({
    # Clear environment (except essential variables)
    rm(list = setdiff(ls(), c("input_dir", "output_dir", "rds_files", "i", "file_name", "input_path", "output_path")))
    gc()
    
    # Load the Seurat object
    cat("  Loading RDS file...\n")
    seurat_obj <- readRDS(input_path)
    
    # Apply normalization and processing steps
    cat("  Normalizing data...\n")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    
    cat("  Finding variable features...\n")
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    
    cat("  Scaling data...\n")
    seurat_obj <- ScaleData(seurat_obj)
    
    cat("  Running PCA...\n")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    
    cat("  Running UMAP...\n")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
    
    # Save the processed object
    cat("  Saving processed object...\n")
    saveRDS(seurat_obj, file = output_path)
    
    # Clear memory
    rm(seurat_obj)
    gc()
    
    cat("  Successfully processed:", file_name, "\n\n")
    
  }, error = function(e) {
    cat("  ERROR processing", file_name, ":", e$message, "\n\n")
    
    # Clean up in case of error
    if (exists("seurat_obj")) {
      rm(seurat_obj)
    }
    gc()
  })
}

cat("Batch processing completed!\n")
cat("Processed files saved to:", output_dir, "\n")