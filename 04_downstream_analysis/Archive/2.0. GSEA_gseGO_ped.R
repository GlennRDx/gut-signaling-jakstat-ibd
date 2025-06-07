# Required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)

# Define the base directory for DEG results
deg_results_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DEG_Results"

# Function to calculate optimal plot height based on number of GO terms
calculate_plot_height <- function(df, base_height = 4, height_per_term = 0.3, min_height = 6, max_height = 20) {
  # Count the number of GO terms per facet (ont + direction combination)
  terms_per_facet <- df %>%
    group_by(ont, direction) %>%
    summarise(n_terms = n(), .groups = 'drop')
  
  # Find the facet with the most terms (this will determine the height)
  max_terms <- max(terms_per_facet$n_terms)
  
  # Calculate height: base height + (number of terms * height per term)
  calculated_height <- base_height + (max_terms * height_per_term)
  
  # Apply min and max constraints
  final_height <- max(min_height, min(max_height, calculated_height))
  
  return(final_height)
}

# Modified function to run GO enrichment with dynamic height
run_go_enrichment_all_dynamic <- function(base_dir, n_terms = 50) {
  
  # Get list of comparisons
  comparisons <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  
  # Keep track of successful and failed analyses
  success_count <- 0
  failure_count <- 0
  
  for (comparison in comparisons) {
    cat("\n---------------------------------------------\n")
    cat("Processing comparison:", comparison, "\n")
    
    # Get cell types for this comparison
    comparison_dir <- file.path(base_dir, comparison)
    cell_types <- list.dirs(comparison_dir, recursive = FALSE, full.names = FALSE)
    
    for (cell_type in cell_types) {
      cat("\n*** Processing cell type:", cell_type, "\n")
      
      # Use tryCatch to handle any errors and continue processing
      tryCatch({
        # DEG results file path
        deg_file <- file.path(comparison_dir, cell_type, "DEG_results.csv")
        
        # Check if file exists
        if (!file.exists(deg_file)) {
          cat("DEG results file not found:", deg_file, "\n")
          next
        }
        
        # Output plot file path
        output_file <- file.path(comparison_dir, cell_type, "GO_enrichment_plot.pdf")
        
        # Read DEG results
        cat("Reading DEG results file:", deg_file, "\n")
        deg_data <- read.csv(deg_file)
        
        # Rename 'gene' column to 'X' if it exists
        if ("gene" %in% colnames(deg_data)) {
          deg_data$X <- deg_data$gene
        }
        
        # Run GO enrichment
        cat("Running GO enrichment analysis...\n")
        title <- paste0(cell_type, " - ", comparison)
        go_results <- main_workflow(deg_data, n = n_terms, title = title)
        
        # Check if we got results
        if (is.null(go_results)) {
          cat("No GO enrichment results for cell type:", cell_type, "\n")
          failure_count <- failure_count + 1
          next
        }
        
        # Calculate optimal plot height based on number of GO terms
        plot_height <- calculate_plot_height(go_results$results)
        cat("Calculated plot height:", plot_height, "inches\n")
        
        # Save plot to PDF with dynamic height
        cat("Saving GO enrichment plot to:", output_file, "\n")
        
        # Try to save the plot with different approaches
        plot_saved <- FALSE
        
        # First try with ggsave (using dynamic height)
        tryCatch({
          ggsave(output_file, plot = go_results$plot, 
                 width = 16, height = plot_height, limitsize = FALSE)
          plot_saved <- TRUE
          cat("Plot saved successfully with ggsave\n")
        }, error = function(e) {
          cat("Error with ggsave:", e$message, "\n")
        })
        
        # If ggsave failed, try with pdf() and print() (using dynamic height)
        if (!plot_saved) {
          tryCatch({
            pdf(output_file, width = 16, height = plot_height)
            print(go_results$plot)
            dev.off()
            plot_saved <- TRUE
            cat("Plot saved successfully with pdf()\n")
          }, error = function(e) {
            cat("Error with pdf():", e$message, "\n")
          })
        }
        
        # Save results to CSV
        results_file <- file.path(comparison_dir, cell_type, "GO_enrichment_results.csv")
        write.csv(go_results$results, results_file, row.names = FALSE)
        
        cat("GO enrichment analysis completed for cell type:", cell_type, "\n")
        success_count <- success_count + 1
        
      }, error = function(e) {
        cat("Error processing", cell_type, "in", comparison, ":", e$message, "\n")
        failure_count <- failure_count + 1
      })
    }
  }
  
  cat("\n===========================================\n")
  cat("GO enrichment analysis summary:\n")
  cat("Successfully processed:", success_count, "cell types\n")
  cat("Failed to process:", failure_count, "cell types\n")
  cat("===========================================\n")
}

# Alternative: More sophisticated height calculation considering facets
calculate_plot_height_advanced <- function(df, base_height = 4, height_per_term = 0.25, 
                                           facet_padding = 1, min_height = 6, max_height = 25) {
  # Count terms in each facet
  facet_counts <- df %>%
    group_by(ont, direction) %>%
    summarise(n_terms = n(), .groups = 'drop')
  
  # Calculate total height needed
  # Each facet needs space for its terms plus some padding
  total_height <- base_height
  
  for (i in 1:nrow(facet_counts)) {
    facet_height <- facet_counts$n_terms[i] * height_per_term + facet_padding
    total_height <- total_height + facet_height
  }
  
  # Apply constraints
  final_height <- max(min_height, min(max_height, total_height))
  
  return(final_height)
}

# You can also modify the plotting function to adjust font size based on the number of terms
plot_go_enrichment_adaptive <- function(df, base_font_size = 8, legend_size = 12, 
                                        cluster_dot_size = 8, title = NULL) {
  
  # Check for required columns
  required_cols <- c("ID", "Description", "pvalue", "ont", "direction", "gene_ratio", "gene_count", "cluster")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure all required columns are present and non-null
  df <- df %>%
    filter(!is.na(ID) & !is.na(Description) & !is.na(pvalue) &
             !is.na(ont) & !is.na(direction) & !is.na(gene_ratio) &
             !is.na(gene_count) & !is.na(cluster))
  
  # Adaptive font size based on number of terms
  n_terms <- nrow(df)
  if (n_terms > 30) {
    font_size <- base_font_size * 0.8
  } else if (n_terms > 20) {
    font_size <- base_font_size * 0.9
  } else {
    font_size <- base_font_size
  }
  
  # Generate a color palette based on the number of clusters
  n_clusters <- length(unique(df$cluster))
  color_palette <- scales::hue_pal()(n_clusters)
  
  # Create labels with p-value formatting and adaptive description length
  max_chars <- if (n_terms > 25) 80 else 100
  
  df$label <- sapply(1:nrow(df), function(i) {
    pval <- df$pvalue[i]
    pval_formatted <- formatC(pval, format = "e", digits = 2)
    stars <- ""
    if (pval < 0.001) {
      stars <- '***'
    } else if (pval < 0.01) {
      stars <- '**'
    } else if (pval < 0.05) {
      stars <- '*'
    }
    
    # Adaptive description truncation
    description <- df$Description[i]
    if (nchar(description) > max_chars) {
      description <- paste0(substr(description, 1, max_chars), "...")
    }
    
    paste0(description, " - ", df$ID[i], " - p-val: ", pval_formatted, " ", stars)
  })
  
  # Rank by cluster, then by descending NES value
  df <- df %>%
    arrange(cluster, abs(NES)) %>%
    mutate(rank = row_number())
  
  # Create the plot
  tryCatch({
    p <- ggplot(df, aes(x = abs(NES), y = reorder(label, rank))) +
      geom_point(aes(color = factor(cluster), size = gene_count), shape = 21, stroke = 0.9, fill = NA) +
      geom_point(aes(color = factor(cluster), size = gene_count), alpha = 0.7) +
      geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold", show.legend = FALSE) +
      facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
      scale_color_manual(values = color_palette) +
      scale_size_continuous(range = c(2, 10)) +
      scale_x_continuous(limits = c(0, NA)) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = font_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        legend.key.size = unit(cluster_dot_size, "mm"),
        strip.text = element_text(size = font_size + 1)  # Slightly larger facet labels
      ) +
      guides(
        color = guide_legend(override.aes = list(size = cluster_dot_size))
      ) +
      labs(x = "Absolute NES", y = "GO Term", color = "Cluster", size = "Gene Count",
           title = title)
    
    return(p)
  }, error = function(e) {
    cat("Error in plotting function:", e$message, "\n")
    cat("Creating simplified plot...\n")
    
    # Create a very simple plot as fallback
    p_simple <- ggplot(df, aes(x = abs(NES), y = Description)) +
      geom_point(aes(color = factor(cluster), size = gene_count)) +
      facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
      scale_color_manual(values = color_palette) +
      theme_bw() +
      labs(x = "Absolute NES", y = "GO Term", color = "Cluster", size = "Gene Count",
           title = title)
    
    return(p_simple)
  })
}

# Run GO enrichment analysis for all DEG results
run_go_enrichment_all(deg_results_dir, n_terms = 50)