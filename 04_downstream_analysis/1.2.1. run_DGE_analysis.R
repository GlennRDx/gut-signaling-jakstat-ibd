# Define directories
data_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/"
output_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results"

# Get RDS files and process each one
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

cat("Found", length(rds_files), "RDS files to process:\n")
cat(paste(basename(rds_files), collapse = ", "), "\n\n")

# ActiveCD vs Control
cat("=== RUNNING ActiveCD vs Control COMPARISONS ===\n\n")
for(file_path in rds_files) {
  # Load Seurat object
  cat("Loading:", basename(file_path), "\n")
  sobj <- readRDS(file_path)
  
  # Get dataset name without extension
  dataset_name <- gsub("\\.rds$", "", basename(file_path))
  
  cat("Processing dataset:", dataset_name, "for ActiveCD vs Control\n")
  
  # Run DGE analysis
  tryCatch({
    run_deseq2_comparison(
      seurat_obj = sobj,
      condition_col = "CD_Status",
      sample_col = "PatientID",
      celltype_col = "Celltypes",
      condition1 = "Control",
      condition2 = "ActiveCD",
      output_dir = output_dir,
      dataset_name = dataset_name
    )
    cat("✓ Successfully completed ActiveCD vs Control analysis for", dataset_name, "\n\n")
  }, error = function(e) {
    cat("✗ Failed ActiveCD vs Control analysis for", dataset_name, ":", e$message, "\n\n")
  })
  
  # Clean up memory
  rm(sobj)
  gc()
}

cat("\n=== RUNNING InactiveCD vs Control COMPARISONS ===\n\n")

# InactiveCD vs Control
for(file_path in rds_files) {
  # Load Seurat object
  cat("Loading:", basename(file_path), "\n")
  sobj <- readRDS(file_path)
  
  # Get dataset name without extension
  dataset_name <- gsub("\\.rds$", "", basename(file_path))
  
  cat("Processing dataset:", dataset_name, "for InactiveCD vs Control\n")
  
  # Run DGE analysis
  tryCatch({
    run_deseq2_comparison(
      seurat_obj = sobj,
      condition_col = "CD_Status",
      sample_col = "PatientID",
      celltype_col = "Celltypes",
      condition1 = "Control",
      condition2 = "InactiveCD",
      output_dir = output_dir,
      dataset_name = dataset_name
    )
    cat("✓ Successfully completed InactiveCD vs Control analysis for", dataset_name, "\n\n")
  }, error = function(e) {
    cat("✗ Failed InactiveCD vs Control analysis for", dataset_name, ":", e$message, "\n\n")
  })
  
  # Clean up memory
  rm(sobj)
  gc()
}

cat("=== ALL ANALYSES COMPLETE ===\n")
cat("Results are organized in the following structure:\n")
cat("└── 05_results_repository/\n")
cat("    ├── ActiveCD_vs_Control/\n")
cat("    │   ├── Colon/\n")
cat("    │   │   ├── EPI/\n")
cat("    │   │   ├── IMM/\n")
cat("    │   │   └── STR/\n")
cat("    │   └── Ileum/\n")
cat("    │       ├── EPI/\n")
cat("    │       ├── IMM/\n")
cat("    │       └── STR/\n")
cat("    └── InactiveCD_vs_Control/\n")
cat("        ├── Colon/\n")
cat("        │   ├── EPI/\n")
cat("        │   ├── IMM/\n")
cat("        │   └── STR/\n")
cat("        └── Ileum/\n")
cat("            ├── EPI/\n")
cat("            ├── IMM/\n")
cat("            └── STR/\n")