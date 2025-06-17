library(Seurat)

# 1. Set your input directory
# rds_files <- list.files(path = "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/", pattern = "\\.rds$", full.names = TRUE)
rds_files <- list.files(path = "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/", pattern = "\\.rds$", full.names = TRUE)


# 2. Create an in-memory mapping (not saved to disk)
celltype_map <- list()

# 3. Function to interactively assign groups
get_group <- function(celltype) {
  repeat {
    cat(paste0("Assign group for cell type '", celltype, "' (E = EPI, I = IMM, S = STR): "))
    ans <- toupper(trimws(readline()))
    if (ans %in% c("E", "I", "S")) {
      return(switch(ans, E = "EPI", I = "IMM", S = "STR"))
    } else {
      cat("❌ Invalid input. Please enter E, I, or S.\n")
    }
  }
}

# 4. Main loop to process RDS files
for (file in rds_files) {
  cat("\n📂 Processing file:", file, "\n")
  
  # Load the Seurat object
  obj <- readRDS(file)
  
  # Identify unique cell types
  celltypes <- unique(obj$Celltypes)
  
  # Prompt user for any unmapped cell types
  for (ct in celltypes) {
    if (is.null(celltype_map[[ct]])) {
      celltype_map[[ct]] <- get_group(ct)
    }
  }
  
  # Assign Cellgroups
  obj$Cellgroups <- sapply(obj$Celltypes, function(ct) celltype_map[[ct]])
  
  # Save 3 separate objects based on Cellgroups
  original_name <- tools::file_path_sans_ext(basename(file))
  for (grp in c("EPI", "IMM", "STR")) {
    subset_obj <- subset(obj, subset = Cellgroups == grp)
    out_path <- file.path(dirname(file), paste0(original_name, "_", grp, ".rds"))
    saveRDS(subset_obj, out_path)
    rm(subset_obj)  # Free memory
    gc()
  }
  
  # Clear the big object to save memory
  rm(obj)
  gc()
}

cat("\n✅ All files processed. Mapping stored in memory only.\n")
