library(Seurat)
library(stringr)
library(dplyr)
library(tibble)

# Define paths
base_dir <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/expression/Raw data/'
metadata_path <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/metadata/scp_metadata_combined.v2.txt'

# Load metadata
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)

# Get subdirectories
subdirs <- list.dirs(base_dir, recursive = FALSE)

for (drctry in subdirs) {
  # Assign the three files to variables
  mtx_file <- list.files(drctry, pattern = "\\.mtx$", full.names = TRUE)
  feature_file <- list.files(drctry, pattern = "features.*\\.tsv$", full.names = TRUE)
  barcode_file <- list.files(drctry, pattern = "barcodes.*\\.tsv$", full.names = TRUE)
  
  # Read matrix
  data <- ReadMtx(
    mtx = mtx_file,
    features = feature_file,
    cells = barcode_file,
    feature.column = 2
  )
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  
  # Match metadata to Seurat cells
  cell_barcodes <- colnames(seurat_obj)
  
  matched_metadata <- metadata %>%
    filter(NAME %in% cell_barcodes) %>%
    column_to_rownames("NAME") %>%
    .[cell_barcodes, , drop = FALSE]  # keep order consistent with Seurat object
  
  # Sanity check
  stopifnot(all(rownames(matched_metadata) == colnames(seurat_obj)))  # optional but good to include
  
  # Add metadata to Seurat object
  seurat_obj <- AddMetaData(seurat_obj, metadata = matched_metadata)
  
  # Save .rds file
  folder_name <- basename(drctry)
  saveRDS(seurat_obj, file = file.path(drctry, paste0(folder_name, ".rds")))
}
