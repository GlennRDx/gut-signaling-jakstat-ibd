plot_seurat_umap_numbered <- function(
    seurat_obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    palette = "Paired",  # Default palette
    figsize = c(9, 9),
    pt.size = 1,
    bg.color = "#f0f0f0",
    label.size = 4,
    title.size = 14,
    axis.text.size = 10,
    legend.text.size = 8,
    title = NULL,
    legend.title = NULL,
    number.clusters = TRUE,
    cluster.text.size = 5,
    cluster.text.color = "black",
    cluster.text.weight = "bold",
    cluster.bg.color = "white",
    cluster.bg.alpha = 0.7,
    custom.order = NULL
) {
  # Load required packages
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(RColorBrewer)
  
  # Extract UMAP coordinates
  umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  
  # Add cell identities
  umap_coords$identity <- seurat_obj@meta.data[[group.by]]
  
  # Get unique identities
  if (!is.null(custom.order)) {
    # Verify all categories are in custom.order
    all_ids <- unique(umap_coords$identity)
    if (!all(all_ids %in% custom.order)) {
      missing <- setdiff(all_ids, custom.order)
      stop(paste("Custom order missing categories:", paste(missing, collapse = ", ")))
    }
    identities <- custom.order
  } else {
    identities <- sort(unique(umap_coords$identity))
  }
  
  # Create color mapping
  num_ids <- length(identities)
  
  # Get palette colors
  if (is.character(palette) && length(palette) == 1) {
    # If the palette is a name, create a color list
    if (palette %in% rownames(brewer.pal.info)) {
      # RColorBrewer palette
      max_colors <- brewer.pal.info[palette, "maxcolors"]
      if (num_ids <= max_colors) {
        colors <- brewer.pal(num_ids, palette)
      } else {
        # Need more colors than the palette provides
        colors <- colorRampPalette(brewer.pal(max_colors, palette))(num_ids)
      }
    } else if (palette == "tab20") {
      # Recreate matplotlib's tab20 palette
      colors <- c(
        "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", 
        "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5", 
        "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", 
        "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
      )
      # Ensure we have enough colors
      if (num_ids > length(colors)) {
        colors <- colorRampPalette(colors)(num_ids)
      } else {
        colors <- colors[1:num_ids]
      }
    } else {
      # Default to a colorRamp if not recognized
      colors <- colorRampPalette(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                                   "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                                   "#bcbd22", "#17becf"))(num_ids)
    }
  } else {
    # Use provided colors
    colors <- palette
    # Ensure we have enough colors
    if (length(colors) < num_ids) {
      colors <- rep(colors, length.out = num_ids)
    } else {
      colors <- colors[1:num_ids]
    }
  }
  
  # Create named vector of colors
  color_mapping <- setNames(colors, identities)
  
  # Calculate cluster centroids for labels
  cluster_centroids <- umap_coords %>%
    group_by(identity) %>%
    summarize(
      UMAP1 = median(UMAP1),
      UMAP2 = median(UMAP2)
    )
  
  # Create numbering for legend
  id_with_numbers <- paste0(seq_along(identities), ". ", identities)
  names(id_with_numbers) <- identities
  
  # Create plot
  p <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = identity)) +
    # Set background color
    theme_bw() +
    theme(
      panel.background = element_rect(fill = bg.color),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.25),
      plot.title = element_text(size = title.size, hjust = 0.5),
      axis.title = element_text(size = axis.text.size),
      axis.text = element_text(size = axis.text.size),
      legend.title = element_text(size = title.size * 0.8),
      legend.text = element_text(size = legend.text.size),
      legend.key.size = unit(0.8, "lines"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.background = element_rect(fill = "white", alpha = 0.8)
    ) +
    # Plot points
    geom_point(size = pt.size, alpha = 0.8) +
    # Set colors
    scale_color_manual(
      values = color_mapping,
      labels = id_with_numbers,
      name = ifelse(is.null(legend.title), group.by, legend.title)
    ) +
    # Add title
    ggtitle(ifelse(is.null(title), paste("UMAP of", group.by), title)) +
    # Set coordinate labels
    xlab("UMAP1") +
    ylab("UMAP2")
  
  # Add numbered cluster labels if requested
  if (number.clusters) {
    # Add cluster numbers
    p <- p + geom_text(
      data = cluster_centroids,
      aes(x = UMAP1, y = UMAP2, label = match(identity, identities)),
      color = cluster.text.color,
      size = cluster.text.size,
      fontface = cluster.text.weight,
      show.legend = FALSE
    )
    
    # Add white circle background for text (to mimic the Python version)
    circle_size <- cluster.text.size * 1.8
    p <- p + geom_point(
      data = cluster_centroids,
      aes(x = UMAP1, y = UMAP2),
      color = "black",
      fill = cluster.bg.color,
      size = circle_size,
      shape = 21,
      alpha = cluster.bg.alpha,
      show.legend = FALSE
    )
    
    # Replot text on top of circles
    p <- p + geom_text(
      data = cluster_centroids,
      aes(x = UMAP1, y = UMAP2, label = match(identity, identities)),
      color = cluster.text.color,
      size = cluster.text.size,
      fontface = cluster.text.weight,
      show.legend = FALSE
    )
  }
  
  # Set plot dimensions
  if (!is.null(figsize)) {
    p <- p + theme(
      aspect.ratio = figsize[2] / figsize[1]
    )
  }
  
  return(p)
}

# Example usage:

rds_path1 <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/expression/Normalised data/CO_IMM_NRM/CO_IMM_NRM.rds"
rds_path2 <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/03_SCP1884/expression/Raw data/CO_IMM_RAW/CO_IMM_RAW.rds"
# Read the Seurat object
seurat_obj1 <- readRDS(rds_path1)
seurat_obj2 <- readRDS(rds_path2)

# Make sure your Seurat object has UMAP and clustering
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# You can rename clusters if needed
new_cluster_ids <- c("BEST4 enterocyte", "Goblet cell", "IL2RG+ enterocyte (M cell)",
                     "Paneth cell", "TA", "Tuft", "crypt", "early enterocyte",
                     "enterocyte", "enteroendocrine")
names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

# Plot with tab20 colors (similar to your Python plot)
plot_seurat_umap_numbered(
  seurat_obj,
  group.by = "ident",  # Use the renamed identities
  palette = "tab20",
  figsize = c(9, 9),
  pt.size = 0.75,
  title = "UMAP of Epithelium Cells with Numbered Clusters"
)

# Save the plot
ggsave("umap_numbered_clusters.png", width = 9, height = 9, dpi = 300)