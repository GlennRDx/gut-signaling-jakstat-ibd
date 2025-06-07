# Batch processing script for all RDS datasets
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)  # Needed for wrap_plots

# Paths
data_path <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/'
output_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/plots/'
features <- c('SOCS1', 'STAT1', 'JAK1', 'JAK2')

# Get all RDS files
rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
dataset_names <- tools::file_path_sans_ext(basename(rds_files))

# Process each dataset
for(i in seq_along(rds_files)) {
  cat("Processing:", dataset_names[i], "\n")
  
  # Load dataset
  seurat_obj <- readRDS(rds_files[i])
  
  # Create output directory
  output_dir <- file.path(output_base, dataset_names[i])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set working directory for plots
  setwd(output_dir)
  
  # UMAP plot
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "Celltypes", raster = F, label = TRUE, label.size = 4) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4, shape = 15)))
  ggsave("umap_plot.pdf", plot = p, width = 13, height = 12, units = "in")
  
  # Feature expression plots using ggplot and patchwork
  seurat_obj$CD_Status <- factor(seurat_obj$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
  
  # Extract and rename UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")  # Ensure consistent naming
  metadata <- seurat_obj@meta.data
  
  # Combine coordinates with metadata
  plot_data <- cbind(umap_coords, metadata)
  
  # Create a list to store individual plots
  plot_list <- list()
  
  # Loop through each feature
  for(feature in features) {
    # Get expression data - Seurat v5 compatible method
    # Try different approaches for getting expression data
    tryCatch({
      # Method 1: Try using FetchData (most compatible)
      expr_data <- FetchData(seurat_obj, vars = feature, layer = "data")
      plot_data[[feature]] <- expr_data[, 1]
    }, error = function(e1) {
      tryCatch({
        # Method 2: Try GetAssayData with specific layer
        if ("data" %in% Layers(seurat_obj[["RNA"]])) {
          plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[feature, ]
        } else {
          # Method 3: Try with counts layer if data doesn't exist
          plot_data[[feature]] <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[feature, ]
        }
      }, error = function(e2) {
        # Method 4: Fallback to older Seurat syntax
        plot_data[[feature]] <- GetAssayData(seurat_obj, slot = "data")[feature, ]
      })
    })
    
    # Create plot
    p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
      geom_point(size = 0.1, alpha = 0.8) +
      scale_color_gradientn(colors = c("grey", "red"), 
                            name = feature,
                            limits = c(0, 2)) +
      facet_wrap(~ CD_Status, nrow = 1) +
      theme_grey() +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.4),
        panel.grid.minor = element_line(color = "white", linewidth = 0.4),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.margin = margin(2,2,2,2, "pt"),
        panel.spacing = unit(0.1, "lines")
      ) +
      coord_fixed()
    
    plot_list[[feature]] <- p
  }
  
  # Combine and save plots
  combined_plot <- wrap_plots(plot_list, ncol = 1)
  ggsave("feature_plot_custom.pdf", plot = combined_plot, width = 15, height = 20, units = "in")
  
  
  # Individual dot plots - also using FetchData for consistency
  for(j in seq_along(features)) {
    # Use FetchData which is more compatible across Seurat versions
    gene_data <- FetchData(seurat_obj, vars = c(features[j], "Celltypes", "CD_Status"))
    colnames(gene_data)[1] <- "expression"
    
    summary_data <- gene_data %>%
      group_by(Celltypes, CD_Status) %>%
      summarise(avg_exp = mean(expression), pct_exp = sum(expression > 0) / n() * 100, .groups = 'drop')
    
    summary_data$CD_Status <- factor(summary_data$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
    
    p <- ggplot(summary_data, aes(x = CD_Status, y = Celltypes)) +
      geom_point(aes(size = pct_exp, color = avg_exp)) +
      scale_color_gradient(low = "lightgrey", high = "blue") +
      scale_size_continuous(range = c(0, 10)) +
      theme_minimal() +
      labs(title = features[j], x = "CD Status", y = "Cell Types", color = "Avg Expression", size = "% Expressed") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(features[j], "_dotplot.pdf"), plot = p, width = 7, height = 10, units = "in")
  }
  
  cat("Completed:", dataset_names[i], "\n")
}