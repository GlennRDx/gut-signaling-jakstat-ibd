# Define directories
data_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/"
output_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/"

# Get RDS files and process each one
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

# Active vs Control
for(file_path in rds_files) {
  # Load Seurat object
  sobj <- readRDS(file_path)
  
  # Get file name without extension for folder prefix
  folder_prefix <- gsub("\\.rds$", "", basename(file_path))
  
  # Run DGE analysis
  run_deseq2_comparison(
    seurat_obj = sobj,
    condition_col = "CD_Status",
    sample_col = "PatientID",
    celltype_col = "Celltypes",
    condition1 = "Control",
    condition2 = "ActiveCD",
    output_dir = output_dir,
    folder_prefix = folder_prefix
  )
}

# Inactive vs Control
for(file_path in rds_files) {
  # Load Seurat object
  sobj <- readRDS(file_path)
  
  # Get file name without extension for folder prefix
  folder_prefix <- gsub("\\.rds$", "", basename(file_path))
  
  # Run DGE analysis
  run_deseq2_comparison(
    seurat_obj = sobj,
    condition_col = "CD_Status",
    sample_col = "PatientID",
    celltype_col = "Celltypes",
    condition1 = "Control",
    condition2 = "InactiveCD",
    output_dir = output_dir,
    folder_prefix = folder_prefix
  )
}