library(Seurat)

# Base paths
base_path <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/expression"
norm_path <- file.path(base_path, "Normalised data")
raw_path  <- file.path(base_path, "Raw data")
output_path <- file.path(base_path, "Preprocessed Raw Data")

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Subdirectory base names
samples <- c("CO_EPI", "CO_IMM", "CO_STR", "TI_EPI", "TI_IMM", "TI_STR")

# Process each pair
for (sample in samples) {
  norm_file <- list.files(file.path(norm_path, paste0(sample, "_NRM")), pattern = "\\.rds$", full.names = TRUE)
  raw_file  <- list.files(file.path(raw_path,  paste0(sample, "_RAW")), pattern = "\\.rds$", full.names = TRUE)
  
  if (length(norm_file) == 1 && length(raw_file) == 1) {
    message("Processing: ", sample)
    
    # Load Seurat objects
    seurat_norm <- readRDS(norm_file)
    seurat_raw  <- readRDS(raw_file)
    
    # Match cells AND features (genes)
    common_cells <- intersect(colnames(seurat_norm), colnames(seurat_raw))
    common_genes <- intersect(rownames(seurat_norm), rownames(seurat_raw))
    
    # Subset both to common features/cells
    seurat_norm <- subset(seurat_norm, features = common_genes, cells = common_cells)
    seurat_raw  <- subset(seurat_raw,  features = common_genes, cells = common_cells)
    
    # Replace counts layer
    raw_counts <- GetAssayData(seurat_raw, assay = "RNA", layer = "counts")
    seurat_norm <- SetAssayData(seurat_norm, assay = "RNA", layer = "counts", new.data = raw_counts)
    
    # Clear data and scale.data layers
    seurat_norm <- SetAssayData(seurat_norm, assay = "RNA", layer = "data", new.data = matrix(nrow = 0, ncol = 0))
    seurat_norm <- SetAssayData(seurat_norm, assay = "RNA", layer = "scale.data", new.data = matrix(nrow = 0, ncol = 0))
    
    # Save updated object
    output_file <- file.path(output_path, paste0(sample, "_preprocessed_raw.rds"))
    saveRDS(seurat_norm, output_file)
    
  } else {
    warning("Missing or ambiguous RDS files for sample: ", sample)
  }
}
