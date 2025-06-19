# Batch processing script for RDS datasets
# Output: Three PNG files per dataset (featureplot, hexbin, dotplot)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)

# Paths and parameters
data_path <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets_norm/'
output_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/plots/dot_feature_hexbin/'
dge_base <- '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results/'
features <- c('SOCS1', 'STAT1', 'JAK1', 'JAK2')

rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
dataset_names <- tools::file_path_sans_ext(basename(rds_files))

overwrite_choice <- NULL

# Parse dataset name for directory structure
parse_dataset_info <- function(dataset_name) {
  parts <- strsplit(dataset_name, "_")[[1]]
  location_full <- ifelse(parts[2] == "CO", "Colon", "Ileum")
  return(list(location_full = location_full, cellgroup = parts[3], dataset_num = parts[1]))
}

# Create output directories
create_output_directories <- function(dataset_info) {
  plot_types <- c("dotplot", "featureplot", "hexbin")
  dirs <- list()
  
  for(plot_type in plot_types) {
    dir_path <- file.path(output_base, dataset_info$location_full, dataset_info$cellgroup, plot_type)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    dirs[[plot_type]] <- dir_path
  }
  
  return(dirs)
}

# Function to load DGE summary data
load_dge_summary <- function(dataset_info, comparison) {
  # Construct path to DGE summary file
  dge_path <- file.path(dge_base, comparison, dataset_info$location_full, 
                        dataset_info$cellgroup, paste0("dataset_", dataset_info$dataset_num))
  
  # Find the DEG_summary CSV file
  csv_files <- list.files(dge_path, pattern = "DEG_summary.*\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("Warning: No DGE summary file found in", dge_path, "\n")
    return(NULL)
  }
  
  if (length(csv_files) > 1) {
    cat("Warning: Multiple DGE summary files found in", dge_path, ". Using first one.\n")
  }
  
  tryCatch({
    dge_data <- read.csv(csv_files[1], stringsAsFactors = FALSE)
    return(dge_data)
  }, error = function(e) {
    cat("Error reading DGE file", csv_files[1], ":", e$message, "\n")
    return(NULL)
  })
}

# Function to determine significance asterisks
get_significance_asterisk <- function(padj) {
  if (is.na(padj)) return("")
  if (padj < 0.001) return("***")
  if (padj < 0.01) return("**")
  if (padj < 0.05) return("*")
  return("")
}

# Function to create significance data frame
create_significance_data <- function(plot_data_comp, dataset_info) {
  # Load DGE data for both comparisons
  inactive_dge <- load_dge_summary(dataset_info, "InactiveCD_vs_Control")
  active_dge <- load_dge_summary(dataset_info, "ActiveCD_vs_Control")
  
  # Initialize significance data
  sig_data <- plot_data_comp
  sig_data$significance <- ""
  sig_data$y_pos <- as.numeric(as.factor(sig_data$Celltypes))
  
  # Process InactiveCD significance
  if (!is.null(inactive_dge)) {
    for (gene in features) {
      padj_col <- paste0(gene, "_padj")
      if (padj_col %in% colnames(inactive_dge)) {
        for (i in 1:nrow(inactive_dge)) {
          cell_type <- inactive_dge$cell_type[i]
          padj_value <- inactive_dge[[padj_col]][i]
          asterisk <- get_significance_asterisk(padj_value)
          
          # Find matching rows in sig_data
          mask <- sig_data$Celltypes == cell_type & 
            sig_data$Gene == gene & 
            sig_data$CD_Status == "InactiveCD"
          
          sig_data$significance[mask] <- asterisk
        }
      }
    }
  }
  
  # Process ActiveCD significance
  if (!is.null(active_dge)) {
    for (gene in features) {
      padj_col <- paste0(gene, "_padj")
      if (padj_col %in% colnames(active_dge)) {
        for (i in 1:nrow(active_dge)) {
          cell_type <- active_dge$cell_type[i]
          padj_value <- active_dge[[padj_col]][i]
          asterisk <- get_significance_asterisk(padj_value)
          
          # Find matching rows in sig_data
          mask <- sig_data$Celltypes == cell_type & 
            sig_data$Gene == gene & 
            sig_data$CD_Status == "ActiveCD"
          
          sig_data$significance[mask] <- asterisk
        }
      }
    }
  }
  
  # Filter to only rows with significance asterisks
  sig_data_filtered <- sig_data[sig_data$significance != "" & 
                                  sig_data$CD_Status %in% c("InactiveCD", "ActiveCD"), ]
  
  return(sig_data_filtered)
}

# Check existing files and handle overwrite
check_and_handle_existing_plots <- function(output_dirs, dataset_name) {
  files_to_check <- c(
    file.path(output_dirs[["featureplot"]], paste0(dataset_name, "_umap_features.png")),
    file.path(output_dirs[["hexbin"]], paste0(dataset_name, "_umap_hexbin.png")),
    file.path(output_dirs[["dotplot"]], paste0(dataset_name, "_dotplot.png"))
  )
  
  existing_files <- files_to_check[file.exists(files_to_check)]
  
  if (length(existing_files) > 0) {
    if (is.null(overwrite_choice)) {
      cat("Plot results already exist. Overwrite existing plots? (Y/N): ")
      user_input <- toupper(trimws(readline()))
      overwrite_choice <<- user_input %in% c("Y", "YES")
    }
    
    if (!overwrite_choice) {
      cat("Skipping (plots exist):", dataset_name, "\n")
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Main processing loop
cat("Starting plot generation for", length(dataset_names), "datasets...\n")
success <- 0

for(i in seq_along(rds_files)) {
  dataset_name <- dataset_names[i]
  dataset_info <- parse_dataset_info(dataset_name)
  output_dirs <- create_output_directories(dataset_info)
  
  if (!check_and_handle_existing_plots(output_dirs, dataset_name)) next
  
  cat("Processing:", dataset_name, "\n")
  
  tryCatch({
    seurat_obj <- readRDS(rds_files[i])
    
    # Prepare plot data
    seurat_obj$CD_Status <- factor(seurat_obj$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
    umap_coords <- Embeddings(seurat_obj, reduction = "umap")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    plot_data <- cbind(umap_coords, seurat_obj@meta.data)
    
    # Add feature expression data
    for(feature in features) {
      expr_data <- FetchData(seurat_obj, vars = feature, layer = "data")
      plot_data[[feature]] <- expr_data[, 1]
    }
    
    # UMAP plot
    umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Celltypes)) +
      geom_point(size = 0.1) +
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
    
    # Feature plots
    feature_plots <- list()
    for(feature in features) {
      p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
        geom_point(size = 0.1, alpha = 0.4) +
        scale_color_gradient(low = "grey", high = "#D81B60", name = feature, limits = c(0, 2)) +
        facet_wrap(~ CD_Status, nrow = 1) +
        theme_void() +
        theme(strip.text = element_text(size = 10), legend.position = "right",
              legend.title = element_text(size = 8), legend.text = element_text(size = 6),
              plot.margin = margin(2,2,2,2, "pt"), panel.spacing = unit(0.1, "lines")) +
        coord_fixed()
      feature_plots[[feature]] <- p
    }
    
    # Save feature plots
    full_plot <- wrap_plots(c(list(umap_plot), feature_plots), ncol = 1, heights = c(1.2, rep(1, length(features))))
    setwd(output_dirs[["featureplot"]])
    ggsave(paste0(dataset_name, "_umap_features.png"), plot = full_plot, width = 14, height = 3 + 3 * length(features), units = "in")
    
    # Hexbin plots
    hex_plots <- list()
    for(feature in features) {
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
    
    # Save hexbin plots
    full_hex_plot <- wrap_plots(c(list(umap_plot), hex_plots), ncol = 1, heights = c(1.2, rep(1, length(features))))
    setwd(output_dirs[["hexbin"]])
    ggsave(paste0(dataset_name, "_umap_hexbin.png"), plot = full_hex_plot, width = 14, height = 3 + 3 * length(features), units = "in")
    
    # Dotplot with significance
    all_vars <- c(features, "Celltypes", "CD_Status")
    gene_data <- FetchData(seurat_obj, vars = all_vars)
    
    plot_data_comp <- gene_data %>%
      pivot_longer(cols = all_of(features), names_to = "Gene", values_to = "Expression") %>%
      group_by(Gene, Celltypes, CD_Status) %>%
      summarise(avg_exp = mean(Expression), 
                pct_exp = sum(Expression > 0) / n() * 100, 
                .groups = 'drop')
    
    plot_data_comp$CD_Status <- factor(plot_data_comp$CD_Status, levels = c("Control", "InactiveCD", "ActiveCD"))
    plot_data_comp$Gene <- factor(plot_data_comp$Gene, levels = features)
    
    # Create significance data
    sig_data <- create_significance_data(plot_data_comp, dataset_info)
    
    max_avg_exp <- max(plot_data_comp$avg_exp, na.rm = TRUE)
    max_pct_exp <- max(plot_data_comp$pct_exp, na.rm = TRUE)
    n_cd_status <- length(levels(plot_data_comp$CD_Status))
    n_genes <- length(features)
    
    p <- ggplot(plot_data_comp, aes(x = interaction(CD_Status, Gene), y = Celltypes)) +
      geom_point(aes(size = pct_exp, color = avg_exp)) +
      scale_color_gradient(low = "grey90", high = "blue", name = "Average\nExpression", limits = c(0, max_avg_exp)) +
      scale_size_continuous(range = c(0.5, 8), name = "% Cells\nExpressed", limits = c(0, max_pct_exp)) +
      scale_x_discrete(labels = rep(levels(plot_data_comp$CD_Status), n_genes)) +
      theme_minimal() +
      labs(title = paste("Gene Expression Profile -", dataset_name), 
           x = "", y = "Cell Types") +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 30)),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            legend.position = "right",
            legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
      coord_cartesian(clip = "off")
    
    # Add alternating background
    for(j in seq_along(features)) {
      if(j %% 2 == 1) {
        x_start <- (j - 1) * n_cd_status + 0.5
        x_end <- j * n_cd_status + 0.5
        p <- p + annotate("rect", xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.3)
      }
    }
    
    # Add gene labels
    gene_label_positions <- sapply(seq_along(features), function(j) mean(((j-1) * n_cd_status + 1):(j * n_cd_status)))
    p <- p + annotate("text", x = gene_label_positions, y = Inf, label = features, vjust = -0.5, size = 4, fontface = "bold")
    
    # Re-add points on top
    p <- p + geom_point(aes(size = pct_exp, color = avg_exp))
    
    # Add significance brackets and asterisks
    if (nrow(sig_data) > 0) {
      # Calculate positions for significance elements
      sig_data$x_pos <- NA
      sig_data$control_x_pos <- NA
      
      for (row_idx in 1:nrow(sig_data)) {
        gene_idx <- which(features == sig_data$Gene[row_idx])
        cd_status_idx <- which(levels(plot_data_comp$CD_Status) == sig_data$CD_Status[row_idx])
        control_idx <- which(levels(plot_data_comp$CD_Status) == "Control")
        
        sig_data$x_pos[row_idx] <- (gene_idx - 1) * n_cd_status + cd_status_idx
        sig_data$control_x_pos[row_idx] <- (gene_idx - 1) * n_cd_status + control_idx
      }
      
      # Create bracket data for lines
      bracket_data <- data.frame()
      
      for (row_idx in 1:nrow(sig_data)) {
        y_pos <- sig_data$y_pos[row_idx]
        x_start <- sig_data$control_x_pos[row_idx]
        x_end <- sig_data$x_pos[row_idx]
        cd_status <- sig_data$CD_Status[row_idx]
        
        # Set different heights for InactiveCD vs ActiveCD to prevent overlap
        if (cd_status == "InactiveCD") {
          bracket_height <- 0.15  # Lower height for InactiveCD
          vertical_line_height <- 0.08
        } else if (cd_status == "ActiveCD") {
          bracket_height <- 0.32  # Higher height for ActiveCD to avoid overlap
          vertical_line_height <- 0.08
        }
        
        # Create bracket segments: horizontal line, two short vertical lines
        bracket_segments <- data.frame(
          x = c(x_start, x_start, x_end, x_end),
          xend = c(x_start, x_end, x_end, x_end),
          y = c(y_pos + bracket_height - vertical_line_height, y_pos + bracket_height, y_pos + bracket_height, y_pos + bracket_height - vertical_line_height),
          yend = c(y_pos + bracket_height, y_pos + bracket_height, y_pos + bracket_height - vertical_line_height, y_pos + bracket_height),
          celltype = rep(sig_data$Celltypes[row_idx], 4),
          gene = rep(sig_data$Gene[row_idx], 4),
          significance = rep(sig_data$significance[row_idx], 4),
          cd_status = rep(cd_status, 4)
        )
        
        bracket_data <- rbind(bracket_data, bracket_segments)
      }
      
      # Add bracket lines
      if (nrow(bracket_data) > 0) {
        p <- p + geom_segment(data = bracket_data,
                              aes(x = x, y = y, xend = xend, yend = yend),
                              inherit.aes = FALSE,
                              color = "black",
                              size = 0.5)
      }
      
      # Add significance text at the center of each bracket
      sig_text_data <- sig_data
      sig_text_data$x_center <- (sig_data$control_x_pos + sig_data$x_pos) / 2
      # Use different heights for text based on CD_Status
      sig_text_data$y_text <- ifelse(sig_data$CD_Status == "InactiveCD", 
                                     sig_data$y_pos + 0.15,  # Lower for InactiveCD
                                     sig_data$y_pos + 0.32)  # Higher for ActiveCD
      
      p <- p + geom_text(data = sig_text_data, 
                         aes(x = x_center, y = y_text, label = significance),
                         inherit.aes = FALSE,
                         size = 4, 
                         color = "black", 
                         fontface = "bold",
                         vjust = -0.2)
    }
    
    # Save dotplot with adjusted width
    setwd(output_dirs[["dotplot"]])
    ggsave(paste0(dataset_name, "_dotplot.png"), plot = p, width = 13, height = 12, units = "in", dpi = 300, bg = "white")
    
    success <- success + 1
    cat("Completed:", dataset_name, "\n")
    
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", e$message, "\n")
  })
  
  gc(verbose = FALSE)
}