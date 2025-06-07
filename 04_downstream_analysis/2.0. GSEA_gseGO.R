# Required libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Changed from org.Mm.eg.db (mouse) to org.Hs.eg.db (human)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)

# Define the base directory for DEG results
deg_results_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DEG_Results"

# Function to perform GO enrichment analysis using gseGO
perform_go_enrichment <- function(df, ont) {
  # Remove rows with NA gene symbols
  df <- df[!is.na(df$X) & df$X != "", ]
  
  if (nrow(df) == 0) {
    cat("    No valid genes after filtering NAs\n")
    return(NULL)
  }
  
  # Convert gene symbols to Entrez IDs (using human database)
  cat("    Converting", nrow(df), "gene symbols to Entrez IDs...\n")
  entrez_ids <- mapIds(org.Hs.eg.db, keys = df$X, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  # Create a data frame to track mappings
  mapping_df <- data.frame(
    gene_symbol = df$X,
    entrez_id = entrez_ids,
    log2fc = df$log2FoldChange,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA Entrez IDs
  mapping_df <- mapping_df[!is.na(mapping_df$entrez_id), ]
  
  if (nrow(mapping_df) == 0) {
    cat("    No genes could be mapped to Entrez IDs\n")
    return(NULL)
  }
  
  cat("    Successfully mapped", nrow(mapping_df), "genes to Entrez IDs\n")
  
  # Remove duplicates (keep first occurrence for each Entrez ID)
  mapping_df <- mapping_df[!duplicated(mapping_df$entrez_id), ]
  
  # Create a named vector of logFC values
  geneList <- mapping_df$log2fc
  names(geneList) <- mapping_df$entrez_id
  
  # Remove any remaining NAs
  geneList <- geneList[!is.na(geneList) & !is.na(names(geneList))]
  
  if (length(geneList) == 0) {
    cat("    No valid gene list after filtering\n")
    return(NULL)
  }
  
  # Sort geneList in decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  
  cat("    Running GSEA with", length(geneList), "genes for ontology:", ont, "\n")
  
  # Perform GSEA
  set.seed(1)
  tryCatch({
    enrich <- gseGO(geneList = geneList, 
                    OrgDb = org.Hs.eg.db,  # Using human database
                    ont = ont, 
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
    
    if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) {
      cat("    No significant results for ontology:", ont, "\n")
      return(NULL)
    }
    
    cat("    Found", nrow(as.data.frame(enrich)), "significant terms for ontology:", ont, "\n")
    return(enrich)
    
  }, error = function(e) {
    cat("    Error in GSEA for ontology", ont, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to extract top GO terms from the results
extract_top_go_terms <- function(enrich, n) {
  result <- as.data.frame(enrich)
  result <- result %>%
    arrange(pvalue) %>%
    slice_head(n = n)
  
  result$direction <- ifelse(result$NES > 0, "Up", "Down")
  
  return(result)
}

# Function to create a combined dataframe with GO terms, gene ratios, and gene lists
create_combined_df <- function(df, ont_list, n) {
  result_df <- data.frame()
  
  for (ont in ont_list) {
    cat("  Processing ontology:", ont, "\n")
    enrich <- perform_go_enrichment(df, ont)
    
    # Check if enrichment was successful and has results
    if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) {
      cat("    No enrichment results for ontology:", ont, "\n")
      next
    }
    
    top_terms <- extract_top_go_terms(enrich, n)
    
    if (nrow(top_terms) > 0) {
      top_terms$ont <- ont
      top_terms$gene_ratio <- top_terms$setSize / nrow(df)
      top_terms$gene_count <- top_terms$setSize
      
      result_df <- rbind(result_df, top_terms)
    }
  }
  
  return(result_df)
}

# Clustering function using Jaccard distance
cluster_go_terms <- function(go_data, h = 0.9) {
  if (nrow(go_data) == 0) {
    return(data.frame(ID = character(0), cluster = integer(0)))
  }
  
  # Extract GO terms and their associated gene lists
  go_terms <- go_data$ID
  genes_list <- strsplit(go_data$core_enrichment, "/")
  
  # Filter out GO terms with no genes
  valid_indices <- sapply(genes_list, function(gene_data) !is.null(gene_data) && length(gene_data) > 0)
  valid_go_terms <- go_terms[valid_indices]
  valid_genes_list <- genes_list[valid_indices]
  
  # If no valid GO terms or only one term, return simple clustering
  if (length(valid_go_terms) <= 1) {
    return(data.frame(ID = valid_go_terms, cluster = rep(1, length(valid_go_terms))))
  }
  
  # Define Jaccard similarity function
  jaccard_similarity <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) return(0)
    return(intersection / union)
  }
  
  # Initialize similarity matrix
  n <- length(valid_go_terms)
  similarity_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(valid_go_terms, valid_go_terms))
  
  # Calculate pairwise Jaccard similarities
  for (i in 1:n) {
    for (j in i:n) {
      sim <- jaccard_similarity(valid_genes_list[[i]], valid_genes_list[[j]])
      similarity_matrix[i, j] <- sim
      similarity_matrix[j, i] <- sim  # Symmetric matrix
    }
  }
  
  # Convert similarity matrix to distance matrix (1 - similarity)
  distance_matrix <- 1 - similarity_matrix
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(distance_matrix), method = "average")
  
  # Cut the dendrogram into clusters
  clusters <- cutree(hc, h = h)
  
  # Return dataframe with valid GO terms and their clusters
  result_df <- data.frame(ID = valid_go_terms, cluster = clusters)
  
  return(result_df)
}

# Restored plotting function with original styling and p-values
plot_go_enrichment <- function(df, font_size = 8, legend_size = 12, cluster_dot_size = 8, title = NULL) {
  # Check for required columns
  required_cols <- c("ID", "Description", "pvalue", "ont", "direction", "gene_ratio", "gene_count", "cluster", "NES")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    cat("    Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("    Available columns:", paste(colnames(df), collapse = ", "), "\n")
    return(NULL)
  }
  
  # Ensure all required columns are present and non-null
  df <- df %>%
    filter(!is.na(ID) & !is.na(Description) & !is.na(pvalue) & 
             !is.na(ont) & !is.na(direction) & !is.na(gene_ratio) & 
             !is.na(gene_count) & !is.na(cluster) & !is.na(NES))
  
  if (nrow(df) == 0) {
    cat("    No valid data for plotting after filtering.\n")
    return(NULL)
  }
  
  cat("    Creating plot with", nrow(df), "GO terms\n")
  
  # Generate a color palette based on the number of clusters (using original hue palette)
  n_clusters <- length(unique(df$cluster))
  color_palette <- scales::hue_pal()(n_clusters)
  
  # Create labels with p-value formatting and asterisks (keeping truncation)
  df$label <- sapply(1:nrow(df), function(i) {
    # Truncate description to 50 characters like in the second script
    description <- df$Description[i]
    truncated_desc <- ifelse(nchar(description) > 50, 
                             paste0(substr(description, 1, 50), "..."), 
                             description)
    
    pval <- df$pvalue[i]
    pval_formatted <- formatC(pval, format = "e", digits = 2)
    stars <- ""
    if (pval < 0.001) {
      stars <- '<span style="font-family: monospace;">***</span>'
    } else if (pval < 0.01) {
      stars <- '<span style="font-family: monospace;">** </span>'
    } else if (pval < 0.05) {
      stars <- '<span style="font-family: monospace;"> * </span>'
    }
    paste0(truncated_desc, " - ", df$ID[i], " - p-val: ", pval_formatted, " ", stars)
  })
  
  # Rank by cluster, then by DESCENDING absolute NES value (fixed ordering)
  df <- df %>%
    arrange(cluster, desc(abs(NES))) %>%
    mutate(rank = row_number())
  
  # Create the plot with original styling
  tryCatch({
    p <- ggplot(df, aes(x = abs(NES), y = reorder(label, rank))) +
      geom_point(aes(color = factor(cluster), size = gene_count), shape = 21, stroke = 0.9, fill = NA) +  # Outline
      geom_point(aes(color = factor(cluster), size = gene_count), alpha = 0.7) +  # Filled dot
      geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold", show.legend = FALSE) +  # Cluster number
      facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
      scale_color_manual(values = color_palette) +
      scale_size_continuous(range = c(2, 10)) +
      scale_x_continuous(limits = c(0, NA)) +  # Set x-axis to start from 0
      theme_bw() +
      theme(
        axis.text.y = ggtext::element_markdown(size = font_size),  # Using element_markdown for HTML formatting
        axis.text.x = element_text(size = font_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        legend.key.size = unit(cluster_dot_size, "mm"),
        strip.text = element_text(size = font_size + 1),
        plot.title = element_text(size = font_size + 4)
      ) +
      guides(
        color = guide_legend(override.aes = list(size = cluster_dot_size))
      ) +
      labs(x = "Absolute NES", y = "GO Term", color = "Cluster", size = "Gene Count", 
           title = title)
    
    return(p)
    
  }, error = function(e) {
    cat("    Error creating plot:", e$message, "\n")
    return(NULL)
  })
}

# Function to calculate optimal plot height based on number of GO terms
calculate_plot_height <- function(df, base_height = 4, height_per_term = 0.3, min_height = 6, max_height = 20) {
  if (is.null(df) || nrow(df) == 0) {
    return(min_height)
  }
  
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

# Enhanced data validation function
validate_deg_data <- function(df) {
  cat("    Validating DEG data...\n")
  
  # Check if df has the required columns
  if (!"X" %in% colnames(df)) {
    if ("gene" %in% colnames(df)) {
      df$X <- df$gene
      cat("    Using 'gene' column as gene symbols\n")
    } else if ("Gene" %in% colnames(df)) {
      df$X <- df$Gene
      cat("    Using 'Gene' column as gene symbols\n")
    } else if ("symbol" %in% colnames(df)) {
      df$X <- df$symbol
      cat("    Using 'symbol' column as gene symbols\n")
    } else {
      cat("    Error: No gene symbol column found (expected 'X', 'gene', 'Gene', or 'symbol')\n")
      cat("    Available columns:", paste(colnames(df), collapse = ", "), "\n")
      return(NULL)
    }
  }
  
  # Check for log2FoldChange column
  if (!"log2FoldChange" %in% colnames(df)) {
    # Try common alternative column names
    if ("logFC" %in% colnames(df)) {
      df$log2FoldChange <- df$logFC
      cat("    Using 'logFC' column as log2FoldChange\n")
    } else if ("log2FC" %in% colnames(df)) {
      df$log2FoldChange <- df$log2FC
      cat("    Using 'log2FC' column as log2FoldChange\n")
    } else if ("log2_fold_change" %in% colnames(df)) {
      df$log2FoldChange <- df$log2_fold_change
      cat("    Using 'log2_fold_change' column as log2FoldChange\n")
    } else {
      cat("    Error: No log2FoldChange column found\n")
      cat("    Available columns:", paste(colnames(df), collapse = ", "), "\n")
      return(NULL)
    }
  }
  
  # Remove rows with NA or empty gene symbols
  initial_rows <- nrow(df)
  df <- df[!is.na(df$X) & df$X != "" & df$X != "NA", ]
  after_gene_filter <- nrow(df)
  
  if (after_gene_filter < initial_rows) {
    cat("    Removed", initial_rows - after_gene_filter, "rows with missing/invalid gene symbols\n")
  }
  
  # Remove rows with NA log2FoldChange
  df <- df[!is.na(df$log2FoldChange) & is.finite(df$log2FoldChange), ]
  after_fc_filter <- nrow(df)
  
  if (after_fc_filter < after_gene_filter) {
    cat("    Removed", after_gene_filter - after_fc_filter, "rows with missing/invalid log2FoldChange values\n")
  }
  
  if (nrow(df) == 0) {
    cat("    Error: No valid data remaining after filtering\n")
    return(NULL)
  }
  
  cat("    Final dataset:", nrow(df), "genes\n")
  return(df)
}

# Main workflow function with enhanced debugging
main_workflow <- function(df, ont_list = c("BP", "CC", "MF"), n = 20, h = 0.95, title = NULL) {
  # Validate and clean the data
  df <- validate_deg_data(df)
  if (is.null(df)) {
    return(NULL)
  }
  
  # Perform enrichment analysis
  cat("  Starting GO enrichment analysis...\n")
  enrichment_results <- create_combined_df(df, ont_list, n)
  
  if (nrow(enrichment_results) == 0) {
    cat("    No enrichment results found\n")
    return(NULL)
  }
  
  cat("  Found", nrow(enrichment_results), "total enriched terms\n")
  
  # Debug: Check what columns we have
  cat("  Enrichment results columns:", paste(colnames(enrichment_results), collapse = ", "), "\n")
  
  # Perform clustering
  cat("  Clustering GO terms...\n")
  clustered_results <- cluster_go_terms(enrichment_results, h)
  
  if (nrow(clustered_results) == 0) {
    cat("    No clustering results\n")
    return(NULL)
  }
  
  # Merge enrichment and clustering results
  final_results <- merge(enrichment_results, clustered_results, by = "ID", all.x = TRUE)
  
  # Replace any remaining NA clusters with a default cluster
  final_results$cluster[is.na(final_results$cluster)] <- max(final_results$cluster, na.rm = TRUE) + 1
  
  # Debug: Check final results structure
  cat("  Final results columns:", paste(colnames(final_results), collapse = ", "), "\n")
  cat("  Final results rows:", nrow(final_results), "\n")
  
  # Check if we have required columns for plotting
  required_cols <- c("ID", "Description", "pvalue", "ont", "direction", "gene_ratio", "gene_count", "cluster", "NES")
  missing_cols <- setdiff(required_cols, colnames(final_results))
  if (length(missing_cols) > 0) {
    cat("  Missing columns for plotting:", paste(missing_cols, collapse = ", "), "\n")
    return(list(results = final_results, plot = NULL))
  }
  
  # Create the plot
  cat("  Creating visualization...\n")
  go_plot <- plot_go_enrichment(final_results, font_size = 10, title = title)
  
  if (is.null(go_plot)) {
    cat("    Failed to create plot but returning results\n")
    return(list(results = final_results, plot = NULL))
  }
  
  cat("  Successfully created plot\n")
  return(list(results = final_results, plot = go_plot))
}

# Enhanced function to run GO enrichment with detailed progress reporting
run_go_enrichment_all_dynamic <- function(base_dir, n_terms = 50) {
  
  # Get list of comparisons
  comparisons <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  
  # Keep track of successful and failed analyses
  success_count <- 0
  failure_count <- 0
  total_comparisons <- length(comparisons)
  
  cat("===========================================\n")
  cat("STARTING GO ENRICHMENT ANALYSIS\n")
  cat("Total comparisons to process:", total_comparisons, "\n")
  cat("===========================================\n\n")
  
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    cat("---------------------------------------------\n")
    cat("COMPARISON", i, "of", total_comparisons, ":", comparison, "\n")
    cat("---------------------------------------------\n")
    
    # Get cell types for this comparison
    comparison_dir <- file.path(base_dir, comparison)
    cell_types <- list.dirs(comparison_dir, recursive = FALSE, full.names = FALSE)
    total_cell_types <- length(cell_types)
    
    cat("Cell types found:", total_cell_types, "\n")
    if (total_cell_types == 0) {
      cat("No cell types found in this comparison. Skipping...\n\n")
      next
    }
    
    for (j in seq_along(cell_types)) {
      cell_type <- cell_types[j]
      cat("\n  Processing cell type", j, "of", total_cell_types, ":", cell_type, "\n")
      
      # Use tryCatch to handle any errors and continue processing
      tryCatch({
        # DEG results file path
        deg_file <- file.path(comparison_dir, cell_type, "DEG_results.csv")
        
        # Check if file exists
        if (!file.exists(deg_file)) {
          cat("    ❌ DEG results file not found:", deg_file, "\n")
          failure_count <- failure_count + 1
          next
        }
        
        # Output file paths
        output_plot_file <- file.path(comparison_dir, cell_type, "GO_enrichment_plot.pdf")
        output_results_file <- file.path(comparison_dir, cell_type, "GO_enrichment_results.csv")
        
        # Read DEG results
        cat("    📖 Reading DEG results file...\n")
        deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)
        cat("    📖 Found", nrow(deg_data), "genes in DEG results\n")
        
        # Run GO enrichment
        cat("    🔬 Running GO enrichment analysis...\n")
        title <- paste0(cell_type, " - ", comparison)
        go_results <- main_workflow(deg_data, n = n_terms, title = title)
        
        # Check if we got results
        if (is.null(go_results)) {
          cat("    ❌ No GO enrichment results generated\n")
          failure_count <- failure_count + 1
          next
        }
        
        # Check if we have actual results
        if (is.null(go_results$results) || nrow(go_results$results) == 0) {
          cat("    ❌ No enriched GO terms found\n")
          failure_count <- failure_count + 1
          next
        }
        
        cat("    📊 Found", nrow(go_results$results), "enriched GO terms\n")
        
        # Always save the results CSV first
        cat("    💾 Saving GO enrichment results...\n")
        write.csv(go_results$results, output_results_file, row.names = FALSE)
        
        # Check if CSV was saved
        if (file.exists(output_results_file)) {
          cat("    ✅ Results successfully saved:", basename(output_results_file), "\n")
        } else {
          cat("    ❌ Failed to save results CSV\n")
        }
        
        # Only try to save plot if we have one
        if (!is.null(go_results$plot)) {
          
          # Save plot to PDF with better error handling
          cat("    💾 Saving GO enrichment plot...\n")
          
          # Calculate optimal plot height based on number of GO terms
          plot_height <- calculate_plot_height(go_results$results)
          plot_width <- 16
          cat("    📏 Plot dimensions:", plot_width, "x", plot_height, "inches\n")
          
          # Method 1: Try with cairo_pdf for better compatibility
          plot_saved <- FALSE
          tryCatch({
            cairo_pdf(output_plot_file, width = plot_width, height = plot_height)
            print(go_results$plot)
            dev.off()
            plot_saved <- TRUE
            cat("    ✅ Plot saved using cairo_pdf\n")
          }, error = function(e) {
            cat("    ⚠️  Error with cairo_pdf:", e$message, "\n")
            # Make sure device is closed
            if (dev.cur() != 1) dev.off()
          })
          
          # Method 2: Try with regular pdf() if cairo_pdf failed
          if (!plot_saved) {
            tryCatch({
              pdf(output_plot_file, width = plot_width, height = plot_height)
              print(go_results$plot)
              dev.off()
              plot_saved <- TRUE
              cat("    ✅ Plot saved using pdf()\n")
            }, error = function(e) {
              cat("    ⚠️  Error with pdf():", e$message, "\n")
              # Make sure device is closed
              if (dev.cur() != 1) dev.off()
            })
          }
          
          # Method 3: Try with ggsave as last resort
          if (!plot_saved) {
            tryCatch({
              ggsave(output_plot_file, plot = go_results$plot, 
                     width = plot_width, height = plot_height, limitsize = FALSE,
                     device = "pdf")
              plot_saved <- TRUE
              cat("    ✅ Plot saved using ggsave\n")
            }, error = function(e) {
              cat("    ⚠️  Error with ggsave:", e$message, "\n")
            })
          }
          
          # Check if plot file exists and has content
          if (file.exists(output_plot_file)) {
            file_size <- file.info(output_plot_file)$size
            if (file_size > 0) {
              cat("    ✅ Plot successfully saved:", basename(output_plot_file), "(", file_size, "bytes)\n")
            } else {
              cat("    ❌ Plot file is empty (0 bytes)\n")
              # Try to create a simple diagnostic plot
              tryCatch({
                pdf(output_plot_file, width = 8, height = 6)
                plot(1:10, 1:10, main = paste("GO Enrichment Failed for", cell_type))
                text(5, 5, "Plot generation failed\nCheck console for errors", cex = 1.2)
                dev.off()
                cat("    📋 Created diagnostic plot instead\n")
              }, error = function(e) {
                cat("    ❌ Could not create diagnostic plot:", e$message, "\n")
              })
            }
          } else {
            cat("    ⚠️  No plot generated, but results were saved\n")
          }
        }
        
        cat("    ✅ COMPLETED:", cell_type, "analysis\n")
        success_count <- success_count + 1
        
      }, error = function(e) {
        cat("    ❌ ERROR processing", cell_type, ":", e$message, "\n")
        failure_count <- failure_count + 1
      })
    }
    
    cat("\nCompleted comparison:", comparison, "\n")
  }
  
  cat("\n===========================================\n")
  cat("GO ENRICHMENT ANALYSIS SUMMARY\n")
  cat("===========================================\n")
  cat("✅ Successfully processed:", success_count, "cell types\n")
  cat("❌ Failed to process:", failure_count, "cell types\n")
  cat("📊 Total processed:", success_count + failure_count, "cell types\n")
  cat("📈 Success rate:", round(success_count / (success_count + failure_count) * 100, 1), "%\n")
  cat("===========================================\n")
}
  
  # Run GO enrichment analysis for all DEG results
  cat("Starting GO enrichment analysis for all comparisons and cell types...\n\n")
  run_go_enrichment_all_dynamic(deg_results_dir, n_terms = 50)