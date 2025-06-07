library(Seurat)

# Define paths
rds_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/"
umap_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/cluster/"
output_dir <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/test'

# List .rds and umap.coord.txt files
rds_files_raw <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
umap_files_raw <- list.files(umap_dir, pattern = "_umap\\.coords\\.txt$", full.names = TRUE)

# Create maps using base names (without extensions)
rds_files <- setNames(rds_files_raw, tools::file_path_sans_ext(basename(rds_files_raw)))
umap_files <- setNames(umap_files_raw, sub("_umap\\.coords\\.txt$", "", basename(umap_files_raw)))

# Find matching files by base name
common_bases <- intersect(names(rds_files), names(umap_files))

# Process each matched pair
for (base_name in common_bases) {
  message("Processing: ", base_name)
  
  # Load Seurat object
  seurat_obj <- readRDS(rds_files[[base_name]])
  
  # Load UMAP coordinates
  umap_df <- read.table(umap_files[[base_name]], header = TRUE, stringsAsFactors = FALSE)
  
  # Set NAME column as rownames, and keep only X and Y
  rownames(umap_df) <- umap_df$NAME
  umap_df <- umap_df[, c("X", "Y")]
  
  # Rename and convert to numeric
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$UMAP_1 <- as.numeric(umap_df$UMAP_1)
  umap_df$UMAP_2 <- as.numeric(umap_df$UMAP_2)
  
  # Filter to matching cells
  common_cells <- intersect(rownames(umap_df), colnames(seurat_obj))
  umap_matrix <- as.matrix(umap_df[common_cells, , drop = FALSE])
  
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  
  # Create and add DimReduc object
  seurat_obj[["umap"]] <- CreateDimReducObject(
    embeddings = umap_matrix,
    key = "UMAP_",
    assay = DefaultAssay(seurat_obj)
  )
  
  # Save updated Seurat object
  saveRDS(seurat_obj, file = file.path(output_dir, paste0(base_name, "_with_umap.rds")))
  
  message("✅ Saved: ", base_name, "_with_umap.rds")
}
