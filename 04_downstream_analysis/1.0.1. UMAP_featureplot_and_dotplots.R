# Batch processing script for all RDS datasets - FINAL VERSION with New Directory Structure
# Output: Three PNG files per dataset organized by plot type
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)  # Added for magma color scale

# Paths
data_path <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/'
output_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/plots/dot_feature_hexbin/'
dge_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results/'
features <- c('SOCS1', 'STAT1', 'JAK1', 'JAK2')

rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
dataset_names <- tools::file_path_sans_ext(basename(rds_files))

# Global variable to store user's overwrite choice
overwrite_choice <- NULL

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

# Function to create all necessary output directories
create_output_directories <- function(dataset_info) {
  # Create directories for both comparisons and all plot types
  comparisons <- c("ActiveCD_vs_Control", "InactiveCD_vs_Control")
  plot_types <- c("dotplot", "featureplot", "hexbin")
  
  dirs <- list()
  
  for(comparison in comparisons) {
    for(plot_type in plot_types) {
      dir_path <- file.path(output_base, comparison, dataset_info$location_full, dataset_info$cellgroup, plot_type)
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      dirs[[paste(comparison, plot_type, sep = "_")]] <- dir_path
    }
  }
  
  return(dirs)
}

# Function to check if plot files exist and handle overwrite decision
check_and_handle_existing_plots <- function(output_dirs, dataset_name) {
  # Check all possible output files across all directories
  files_to_check <- c()
  
  # Add all possible filenames
  for(dir_name in names(output_dirs)) {
    dir_path <- output_dirs[[dir_name]]
    
    if(grepl("featureplot", dir_name)) {
      files_to_check <- c(files_to_check, file.path(dir_path, paste0(dataset_name, "_umap_features.png")))
    } else if(grepl("hexbin", dir_name)) {
      files_to_check <- c(files_to_check, file.path(dir_path, paste0(dataset_name, "_umap_hexbin.png")))
    } else if(grepl("dotplot", dir_name)) {
      files_to_check <- c(files_to_check, file.path(dir_path, paste0(dataset_name, "_dotplot.png")))
    }
  }
  
  # Check if any plot files already exist
  existing_files <- files_to_check[file.exists(files_to_check)]
  
  if (length(existing_files) > 0) {
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
      cat("Skipping (plots exist):", dataset_name, "\n")
      return(FALSE)  # Skip processing
    }
  }
  
  return(TRUE)  # Proceed with processing
}

read_dge_results <- function(dataset_name) {
  dataset_info <- parse_dataset_info(dataset_name)
  
  # Build paths according to new structure
  # DGE_Results/ActiveCD_vs_Control/Colon/EPI/dataset_02/Enterocytes_BEST4/
  active_base_path <- file.path(dge_base, "ActiveCD_vs_Control", dataset_info$location_full, dataset_info$cellgroup, dataset_name)
  inactive_base_path <- file.path(dge_base, "InactiveCD_vs_Control", dataset_info$location_full, dataset_info$cellgroup, dataset_name)
  
  dge_data <- data.frame()
  
  # Process ActiveCD results
  if(dir.exists(active_base_path)) {
    # Get all cell type folders in this dataset directory
    celltype_folders <- list.dirs(active_base_path, full.names = TRUE, recursive = FALSE)
    
    for(celltype_folder in celltype_folders) {
      celltype_name <- basename(celltype_folder)
      
      # Look for CSV files in this cell type folder
      csv_files <- list.files(celltype_folder, pattern = "\\.csv$", full.names = TRUE)
      
      for(csv_file in csv_files) {
        if(file.exists(csv_file)) {
          tryCatch({
            active_data <- read.csv(csv_file, stringsAsFactors = FALSE)
            active_data$Comparison <- "ActiveCD"
            active_data$cell_type <- celltype_name  # Add cell type info if not present
            dge_data <- rbind(dge_data, active_data)
          }, error = function(e) {
            cat("Warning: Could not read", csv_file, "- Error:", e$message, "\n")
          })
        }
      }
    }
  }
  
  # Process InactiveCD results
  if(dir.exists(inactive_base_path)) {
    # Get all cell type folders in this dataset directory
    celltype_folders <- list.dirs(inactive_base_path, full.names = TRUE, recursive = FALSE)
    
    for(celltype_folder in celltype_folders) {
      celltype_name <- basename(celltype_folder)
      
      # Look for CSV files in this cell type folder
      csv_files <- list.files(celltype_folder, pattern = "\\.csv$", full.names = TRUE)
      
      for(csv_file in csv_files) {
        if(file.exists(csv_file)) {
          tryCatch({
            inactive_data <- read.csv(csv_file, stringsAsFactors = FALSE)
            inactive_data$Comparison <- "InactiveCD"
            inactive_data$cell_type <- celltype_name  # Add cell type info if not present
            dge_data <- rbind(dge_data, inactive_data)
          }, error = function(e) {
            cat("Warning: Could not read", csv_file, "- Error:", e$message, "\n")
          })
        }
      }
    }
  }
  
  return(dge_data)
}

extract_significant_log2fc <- function(dge_data, gene_name, celltype, comparison) {
  log2fc_col <- paste0(gene_name, "_log2FC")
  if(!log2fc_col %in% colnames(dge_data)) return(list(value = NA, significant = FALSE, original = NA))
  matching_rows <- which(dge_data$cell_type == celltype & dge_data$Comparison == comparison)
  if(length(matching_rows) == 0) return(list(value = NA, significant = FALSE, original = NA))
  row_data <- as.character(dge_data[matching_rows[1], log2fc_col])
  if(grepl("\\*", row_data)) {
    numeric_val <- as.numeric(gsub("\\*+", "", row_data))
    return(list(value = numeric_val, significant = TRUE, original = row_data))
  } else {
    return(list(value = as.numeric(row_data), significant = FALSE, original = row_data))
  }
}

seurat_obj <- NULL

cat("Starting plot generation for", length(dataset_names), "datasets...\n")
success <- 0
total <- 0

for(i in seq_along(rds_files)) {
  total <- total + 1
  dataset_name <- dataset_names[i]
  
  # Parse dataset information
  dataset_info <- parse_dataset_info(dataset_name)
  
  # Create all output directories
  output_dirs <- create_output_directories(dataset_info)
  
  # Check if plots exist and handle overwrite decision
  if (!check_and_handle_existing_plots(output_dirs, dataset_name)) {
    next  # Skip this dataset
  }
  
  cat("Processing:", dataset_name, "\n")
  
  # Clean up previous object
  if (!is.null(seurat_obj)) { rm(seurat_obj); gc(verbose = FALSE) }
  
  tryCatch({
    seurat_obj <- readRDS(rds_files[i])
    
    # UMAP and Feature Plots (Combined)
    seurat_obj$CD_Status <- factor(seurat_obj$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
    umap_coords <- Embeddings(seurat_obj, reduction = "umap")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    metadata <- seurat_obj@meta.data
    plot_data <- cbind(umap_coords, metadata)
    
    # UMAP plot - FIXED: Now with faceting to match feature plots
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
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 4, shape = 15)))
    
    # Feature plots with magma color gradient
    feature_plots <- list()
    for(feature in features) {
      tryCatch({
        expr_data <- FetchData(seurat_obj, vars = feature, layer = "data")
        plot_data[[feature]] <- expr_data[, 1]
      }, error = function(e1) {
        tryCatch({
          if ("data" %in% Layers(seurat_obj[["RNA"]])) {
            plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[feature, ]
          } else {
            plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[feature, ]
          }
        }, error = function(e2) {
          plot_data[[feature]] <- GetAssayData(seurat_obj, slot = "data")[feature, ]
        })
      })
      
      p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
        geom_point(size = 0.1, alpha = 0.4) +
        scale_color_viridis_c(option = "magma", name = feature, limits = c(0, 2)) +
        facet_wrap(~ CD_Status, nrow = 1) +
        theme_void() +
        theme(strip.text = element_text(size = 10), legend.position = "right",
              legend.title = element_text(size = 8), legend.text = element_text(size = 6),
              plot.margin = margin(2,2,2,2, "pt"), panel.spacing = unit(0.1, "lines")) +
        coord_fixed()
      feature_plots[[feature]] <- p
    }
    
    full_plot <- wrap_plots(c(list(umap_plot), feature_plots), ncol = 1, heights = c(1.2, rep(1, length(features))))
    
    # Save feature plots to both ActiveCD_vs_Control and InactiveCD_vs_Control directories
    for(comparison in c("ActiveCD_vs_Control", "InactiveCD_vs_Control")) {
      feature_dir <- output_dirs[[paste(comparison, "featureplot", sep = "_")]]
      setwd(feature_dir)
      umap_filename <- paste0(dataset_name, "_umap_features.png")
      ggsave(umap_filename, plot = full_plot, width = 14, height = 3 + 3 * length(features), units = "in")
    }
    
    rm(umap_plot, feature_plots, full_plot)
    
    # Hexagonal binning feature plots with magma color gradient
    hex_plots <- list()
    for(feature in features) {
      # Ensure feature data is available in plot_data
      if(!feature %in% colnames(plot_data)) {
        tryCatch({
          expr_data <- FetchData(seurat_obj, vars = feature, layer = "data")
          plot_data[[feature]] <- expr_data[, 1]
        }, error = function(e1) {
          tryCatch({
            if ("data" %in% Layers(seurat_obj[["RNA"]])) {
              plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[feature, ]
            } else {
              plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[feature, ]
            }
          }, error = function(e2) {
            plot_data[[feature]] <- GetAssayData(seurat_obj, slot = "data")[feature, ]
          })
        })
      }
      
      p_hex <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, z = .data[[feature]])) +
        stat_summary_hex(bins = 50, fun = mean, alpha = 0.8) +
        scale_fill_viridis_c(option = "magma", name = feature, limits = c(0, 2)) +
        facet_wrap(~ CD_Status, nrow = 1) +
        theme_void() +
        theme(strip.text = element_text(size = 10), 
              legend.position = "right",
              legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6),
              plot.margin = margin(2,2,2,2, "pt"), 
              panel.spacing = unit(0.1, "lines")) +
        coord_fixed()
      hex_plots[[feature]] <- p_hex
    }
    
    # UMAP plot for hexbin version (same as before)
    umap_plot_hex <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Celltypes)) +
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
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 4, shape = 15)))
    
    full_hex_plot <- wrap_plots(c(list(umap_plot_hex), hex_plots), ncol = 1, heights = c(1.2, rep(1, length(features))))
    
    # Save hexbin plots to both ActiveCD_vs_Control and InactiveCD_vs_Control directories
    for(comparison in c("ActiveCD_vs_Control", "InactiveCD_vs_Control")) {
      hexbin_dir <- output_dirs[[paste(comparison, "hexbin", sep = "_")]]
      setwd(hexbin_dir)
      hex_filename <- paste0(dataset_name, "_umap_hexbin.png")
      ggsave(hex_filename, plot = full_hex_plot, width = 14, height = 3 + 3 * length(features), units = "in")
    }
    
    rm(umap_plot_hex, hex_plots, full_hex_plot)
    
    # Generate dotplots for each comparison separately
    dge_results <- read_dge_results(dataset_name)
    all_vars <- c(features, "Celltypes", "CD_Status")
    gene_data <- FetchData(seurat_obj, vars = all_vars)
    
    for(comparison in c("ActiveCD_vs_Control", "InactiveCD_vs_Control")) {
      # Prepare plot data for this comparison
      plot_data_comp <- gene_data %>%
        pivot_longer(cols = all_of(features), names_to = "Gene", values_to = "Expression") %>%
        group_by(Gene, Celltypes, CD_Status) %>%
        summarise(avg_exp = mean(Expression), 
                  pct_exp = sum(Expression > 0) / n() * 100, 
                  .groups = 'drop')
      
      plot_data_comp$CD_Status <- factor(plot_data_comp$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
      plot_data_comp$Gene <- factor(plot_data_comp$Gene, levels = features)
      plot_data_comp$log2fc_text <- ""
      plot_data_comp$is_significant <- FALSE
      
      # Add DGE information for the current comparison
      comparison_type <- ifelse(comparison == "ActiveCD_vs_Control", "ActiveCD", "InactiveCD")
      
      if(nrow(dge_results) > 0) {
        for(idx in 1:nrow(plot_data_comp)) {
          gene <- as.character(plot_data_comp$Gene[idx])
          celltype <- as.character(plot_data_comp$Celltypes[idx])
          cd_status <- as.character(plot_data_comp$CD_Status[idx])
          
          if(cd_status == comparison_type) {
            dge_info <- extract_significant_log2fc(dge_results, gene, celltype, comparison_type)
            if(!is.null(dge_info) && is.list(dge_info) && !is.na(dge_info$value) && dge_info$significant) {
              plot_data_comp$log2fc_text[idx] <- paste0(round(dge_info$value, 2), "")
              plot_data_comp$is_significant[idx] <- TRUE
            }
          }
        }
      }
      
      max_avg_exp <- max(plot_data_comp$avg_exp, na.rm = TRUE)
      max_pct_exp <- max(plot_data_comp$pct_exp, na.rm = TRUE)
      n_cd_status <- length(levels(plot_data_comp$CD_Status))
      n_genes <- length(features)
      
      p <- ggplot(plot_data_comp, aes(x = interaction(CD_Status, Gene), y = Celltypes)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        geom_text(aes(label = ifelse(is_significant, log2fc_text, "")), size = 2.5, color = "white", fontface = "bold") +
        scale_color_viridis_c(option = "magma", name = "Average\nExpression", limits = c(0, max_avg_exp)) +
        scale_size_continuous(range = c(0.5, 8), name = "% Cells\nExpressed", limits = c(0, max_pct_exp)) +
        scale_x_discrete(labels = rep(levels(plot_data_comp$CD_Status), n_genes)) +
        theme_minimal() +
        labs(title = paste("Gene Expression Profile -", comparison, "-", dataset_name), 
             x = "", y = "Cell Types",
             caption = "White text shows significant log2FC values") +
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right",
              legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
        coord_cartesian(clip = "off")
      
      for(j in seq_along(features)) {
        if(j %% 2 == 1) {
          x_start <- (j - 1) * n_cd_status + 0.5
          x_end <- j * n_cd_status + 0.5
          p <- p + annotate("rect", xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = "grey95", alpha = 0.3)
        }
      }
      
      p <- p +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        geom_text(aes(label = ifelse(is_significant, log2fc_text, "")), size = 3.5, color = "black", vjust = 0.5, hjust = 0.5)
      
      gene_label_positions <- sapply(seq_along(features), function(j) mean(((j-1) * n_cd_status + 1):(j * n_cd_status)))
      p <- p + annotate("text", x = gene_label_positions, y = Inf, label = features, vjust = -0.5, size = 4, fontface = "bold")
      
      # Save dotplot to the appropriate directory
      dotplot_dir <- output_dirs[[paste(comparison, "dotplot", sep = "_")]]
      setwd(dotplot_dir)
      dotplot_filename <- paste0(dataset_name, "_dotplot.png")
      ggsave(dotplot_filename, plot = p, width = 16, height = 12, units = "in", dpi = 300, bg = "white")
    }
    
    rm(gene_data, dge_results)
    
    success <- success + 1
    cat("Completed:", dataset_name, "\n")
    
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", e$message, "\n")
  })
  
  gc(verbose = FALSE)
}