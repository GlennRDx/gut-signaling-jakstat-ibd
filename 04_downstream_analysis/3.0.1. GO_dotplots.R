# Required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)
library(digest)  # Added for hash-based caching
library(ggplot2)
library(scales)

# Define the base directory for DEG results
deg_results_dir <- "/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results"

# Cache for gene mappings to avoid repeated database queries
gene_mapping_cache <- new.env()

# Global variable to store user's overwrite choice
overwrite_choice <- NULL

# Efficient gene mapping with caching using hash-based keys
map_genes_cached <- function(gene_symbols) {
  # Create hash-based cache key to avoid variable name length limits
  cache_key <- digest(paste(sort(gene_symbols), collapse = "|"), algo = "md5")
  
  # Check cache first
  if (exists(cache_key, envir = gene_mapping_cache)) {
    return(get(cache_key, envir = gene_mapping_cache))
  }
  
  # Query all genes (simplified approach - cache individual results)
  gene_mappings <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", 
                          keytype = "SYMBOL", multiVals = "first")
  
  # Cache the result for this combination
  assign(cache_key, gene_mappings, envir = gene_mapping_cache)
  return(gene_mappings)
}

# Perform GO enrichment analysis for all three ontologies
perform_go_enrichment <- function(df, n = 50) {
  # Clean data once
  df <- df[!is.na(df$X) & df$X != "" & !is.na(df$log2FoldChange), ]
  if (nrow(df) == 0) return(NULL)
  
  # Single gene mapping call with caching
  entrez_ids <- map_genes_cached(df$X)
  mapping_df <- data.frame(
    gene_symbol = df$X,
    entrez_id = entrez_ids,
    log2fc = df$log2FoldChange,
    stringsAsFactors = FALSE
  )
  
  mapping_df <- mapping_df[!is.na(mapping_df$entrez_id) & !duplicated(mapping_df$entrez_id), ]
  if (nrow(mapping_df) == 0) return(NULL)
  
  geneList <- sort(setNames(mapping_df$log2fc, mapping_df$entrez_id), decreasing = TRUE)
  
  # Run GO enrichment analysis for all three ontologies
  result_df <- data.frame()
  ontologies <- c("BP", "MF", "CC")  # Biological Process, Molecular Function, Cellular Component
  
  for (ont in ontologies) {
    tryCatch({
      go_enrich <- gseGO(geneList = geneList,
                         OrgDb = org.Hs.eg.db,
                         ont = ont,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         verbose = FALSE)
      
      if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
        top_terms <- as.data.frame(go_enrich) %>%
          arrange(pvalue) %>%
          slice_head(n = n) %>%
          mutate(
            direction = ifelse(NES > 0, "Up", "Down"),
            ontology = ont
          )
        
        result_df <- rbind(result_df, top_terms)
      }
    }, error = function(e) {
      cat("Warning: Failed to process GO", ont, "enrichment:", e$message, "\n")
    })
  }
  
  if (nrow(result_df) > 0) {
    return(result_df)
  }
  
  return(NULL)
}

# Optimized main workflow for GO analysis
run_go_analysis <- function(deg_file, output_dir, cell_type, cell_group, tissue_dir, comparison, dataset) {
  if (!file.exists(deg_file)) return(list(success = FALSE, results = NULL))
  
  # Check if GO results already exist
  csv_output <- file.path(output_dir, "GO.csv")
  
  if (file.exists(csv_output)) {
    # If overwrite choice hasn't been made yet, prompt user
    if (is.null(overwrite_choice)) {
      cat("GO results already exist. Overwrite existing results? (Y/N): ")
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
    
    # If user chose not to overwrite, skip this analysis but return existing results
    if (!overwrite_choice) {
      cat("Skipping (results exist):", deg_file, "\n")
      # Read existing results for summary
      existing_results <- tryCatch({
        read.csv(csv_output, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(existing_results)) {
        # Add metadata columns if they don't exist
        if (!"cell_type" %in% colnames(existing_results)) {
          existing_results$cell_type <- cell_type
          existing_results$cell_group <- cell_group
          existing_results$tissue <- tissue_dir
          existing_results$comparison <- comparison
          existing_results$dataset <- dataset
        }
      }
      
      return(list(success = TRUE, results = existing_results))
    }
  }
  
  cat("Processing:", deg_file, "\n")
  
  # Single file read with column standardization
  df <- read.csv(deg_file, stringsAsFactors = FALSE)
  if (!"X" %in% colnames(df)) {
    if ("gene" %in% colnames(df)) df$X <- df$gene
    else {
      cat("Warning: No gene column found in", deg_file, "\n")
      return(list(success = FALSE, results = NULL))
    }
  }
  if (!"log2FoldChange" %in% colnames(df)) {
    if ("logFC" %in% colnames(df)) df$log2FoldChange <- df$logFC
    else {
      cat("Warning: No log2FoldChange column found in", deg_file, "\n")
      return(list(success = FALSE, results = NULL))
    }
  }
  
  # Perform GO enrichment analysis
  results <- perform_go_enrichment(df)
  if (is.null(results) || nrow(results) == 0) {
    cat("No GO enrichment results for", deg_file, "\n")
    return(list(success = FALSE, results = NULL))
  }
  
  # Add metadata columns for summary
  results$cell_type <- cell_type
  results$cell_group <- cell_group
  results$tissue <- tissue_dir
  results$comparison <- comparison
  results$dataset <- dataset
  
  # Save results
  write.csv(results, csv_output, row.names = FALSE)
  
  cat("Completed:", deg_file, "\n")
  return(list(success = TRUE, results = results))
}

# Function to create GO dotplot from summary data
create_go_dotplot <- function(summary_file, output_dir, max_pathways = 30, pvalue_cutoff = 0.05) {
  if (!file.exists(summary_file)) {
    cat("Summary file not found:", summary_file, "\n")
    return(FALSE)
  }
  
  # Read summary data
  go_data <- read.csv(summary_file, stringsAsFactors = FALSE)
  
  if (nrow(go_data) == 0) {
    cat("No data in summary file:", summary_file, "\n")
    return(FALSE)
  }
  
  # Filter by p-value and prepare data
  go_filtered <- go_data %>%
    filter(pvalue <= pvalue_cutoff) %>%
    mutate(
      # Create a combined identifier for GO terms
      go_term_id = paste(ID, Description, sep = " - "),
      # Transform p-value for size (smaller p-value = larger dot)
      neg_log_pval = -log10(pvalue),
      # Cap NES values for better visualization
      NES_capped = pmax(pmin(NES, 3), -3),
      # Format ontology names
      ontology = case_when(
        ontology == "BP" ~ "Biological Process",
        ontology == "MF" ~ "Molecular Function",
        ontology == "CC" ~ "Cellular Component",
        TRUE ~ ontology
      )
    )
  
  if (nrow(go_filtered) == 0) {
    cat("No significant GO terms found (p <=", pvalue_cutoff, ") in:", summary_file, "\n")
    return(FALSE)
  }
  
  # Create separate plots for each ontology
  ontologies <- unique(go_filtered$ontology)
  
  for (ont in ontologies) {
    ont_data <- go_filtered %>%
      filter(ontology == ont)
    
    if (nrow(ont_data) == 0) next
    
    # Select top pathways by significance (lowest p-values) for this ontology
    top_terms <- ont_data %>%
      group_by(go_term_id) %>%
      summarise(min_pval = min(pvalue), .groups = 'drop') %>%
      arrange(min_pval) %>%
      slice_head(n = max_pathways) %>%
      pull(go_term_id)
    
    # Filter data to top terms
    plot_data <- ont_data %>%
      filter(go_term_id %in% top_terms) %>%
      mutate(
        go_term_id = factor(go_term_id, levels = rev(top_terms)),
        cell_type = factor(cell_type)
      )
    
    # Create the dotplot
    p <- ggplot(plot_data, aes(x = cell_type, y = go_term_id)) +
      geom_point(aes(size = neg_log_pval, color = NES_capped), alpha = 0.8) +
      scale_size_continuous(
        name = "-log10(p-value)",
        range = c(1, 8),
        breaks = c(1, 2, 3, 4, 5),
        labels = c("0.1", "0.01", "0.001", "0.0001", "0.00001")
      ) +
      scale_color_gradient2(
        name = "NES",
        low = "red", 
        mid = "white", 
        high = "green",
        midpoint = 0,
        limits = c(-3, 3),
        breaks = c(-3, -1.5, 0, 1.5, 3),
        labels = c("≤-3", "-1.5", "0", "1.5", "≥3")
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Cell Type",
        y = "GO Term",
        title = paste("GO", ont, "Enrichment Across Cell Types"),
        subtitle = paste("Top", length(top_terms), "terms (p ≤", pvalue_cutoff, ")")
      )
    
    # Save the plot
    ont_safe <- gsub(" ", "_", ont)
    plot_file <- file.path(output_dir, paste0("GO_", ont_safe, "_dotplot.pdf"))
    ggsave(plot_file, p, width = max(8, length(unique(plot_data$cell_type)) * 0.8), 
           height = max(6, length(top_terms) * 0.3), limitsize = FALSE)
    
    cat("Created GO", ont, "dotplot:", plot_file, "\n")
  }
  
  # Create combined plot with all ontologies
  if (length(ontologies) > 1) {
    # Select fewer terms per ontology for combined plot
    terms_per_ont <- max(10, floor(max_pathways / length(ontologies)))
    
    combined_data <- go_filtered %>%
      group_by(ontology) %>%
      arrange(pvalue) %>%
      slice_head(n = terms_per_ont) %>%
      ungroup() %>%
      mutate(
        go_term_id = factor(go_term_id, levels = rev(unique(go_term_id))),
        cell_type = factor(cell_type)
      )
    
    p_combined <- ggplot(combined_data, aes(x = cell_type, y = go_term_id)) +
      geom_point(aes(size = neg_log_pval, color = NES_capped), alpha = 0.8) +
      facet_wrap(~ ontology, scales = "free_y", ncol = 1) +
      scale_size_continuous(
        name = "-log10(p-value)",
        range = c(1, 6),
        breaks = c(1, 2, 3, 4, 5),
        labels = c("0.1", "0.01", "0.001", "0.0001", "0.00001")
      ) +
      scale_color_gradient2(
        name = "NES",
        low = "red", 
        mid = "white", 
        high = "green",
        midpoint = 0,
        limits = c(-3, 3)
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Cell Type",
        y = "GO Term",
        title = "GO Enrichment Across Cell Types (All Ontologies)",
        subtitle = paste("Top", terms_per_ont, "terms per ontology (p ≤", pvalue_cutoff, ")")
      )
    
    plot_file_combined <- file.path(output_dir, "GO_Combined_dotplot.pdf")
    ggsave(plot_file_combined, p_combined, 
           width = max(10, length(unique(combined_data$cell_type)) * 0.8), 
           height = max(12, nrow(combined_data) * 0.15), limitsize = FALSE)
    
    cat("Created combined GO dotplot:", plot_file_combined, "\n")
  }
  
  return(TRUE)
}

# Function to create GO summary for each dataset
create_go_summary <- function(dataset_path) {
  summary_file <- file.path(dataset_path, "GO_Summary.csv")
  
  # Check if summary already exists
  if (file.exists(summary_file) && !is.null(overwrite_choice) && !overwrite_choice) {
    cat("GO summary already exists, skipping:", summary_file, "\n")
    return(TRUE)
  }
  
  # Find all GO.csv files in subdirectories
  go_files <- list.files(dataset_path, pattern = "GO\\.csv$", recursive = TRUE, full.names = TRUE)
  
  if (length(go_files) == 0) {
    cat("No GO.csv files found in", dataset_path, "\n")
    return(FALSE)
  }
  
  # Read and combine all GO results
  all_results <- data.frame()
  
  for (go_file in go_files) {
    tryCatch({
      go_data <- read.csv(go_file, stringsAsFactors = FALSE)
      if (nrow(go_data) > 0) {
        all_results <- rbind(all_results, go_data)
      }
    }, error = function(e) {
      cat("Warning: Could not read", go_file, ":", e$message, "\n")
    })
  }
  
  if (nrow(all_results) > 0) {
    # Sort by p-value
    all_results <- all_results %>% arrange(pvalue)
    write.csv(all_results, summary_file, row.names = FALSE)
    cat("Created GO summary:", summary_file, "(", nrow(all_results), "GO terms)\n")
    
    # Create dotplot
    create_go_dotplot(summary_file, dataset_path)
    
    return(TRUE)
  } else {
    cat("No GO results to summarize in", dataset_path, "\n")
    return(FALSE)
  }
}

# Updated directory processing for new structure
process_all_comparisons <- function(base_dir) {
  cat("Starting GO enrichment analysis...\n")
  
  comparisons <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  success <- 0
  total <- 0
  summary_success <- 0
  total_datasets <- 0
  
  for (comparison in comparisons) {
    comparison_dir <- file.path(base_dir, comparison)
    tissue_dirs <- list.dirs(comparison_dir, recursive = FALSE, full.names = FALSE)
    
    for (tissue_dir in tissue_dirs) {
      tissue_path <- file.path(comparison_dir, tissue_dir)
      cell_group_dirs <- list.dirs(tissue_path, recursive = FALSE, full.names = FALSE)
      
      for (cell_group in cell_group_dirs) {
        cell_group_path <- file.path(tissue_path, cell_group)
        dataset_dirs <- list.dirs(cell_group_path, recursive = FALSE, full.names = FALSE)
        
        for (dataset in dataset_dirs) {
          dataset_path <- file.path(cell_group_path, dataset)
          cell_type_dirs <- list.dirs(dataset_path, recursive = FALSE, full.names = FALSE)
          
          # Process all cell types in this dataset
          for (cell_type in cell_type_dirs) {
            cell_type_path <- file.path(dataset_path, cell_type)
            deg_files <- list.files(cell_type_path, pattern = "_DEG_results\\.csv$", full.names = TRUE)
            
            for (deg_file in deg_files) {
              total <- total + 1
              result <- run_go_analysis(deg_file, cell_type_path, cell_type, cell_group, tissue_dir, comparison, dataset)
              if (result$success) success <- success + 1
            }
          }
          
          # Create summary immediately after all cell types in this dataset are processed
          total_datasets <- total_datasets + 1
          cat("Creating GO summary for dataset:", dataset, "\n")
          if (create_go_summary(dataset_path)) {
            summary_success <- summary_success + 1
          }
        }
      }
    }
  }
  
  cat("GO enrichment completed:", success, "of", total, "analyses successful\n")
  cat("GO summaries created:", summary_success, "of", total_datasets, "datasets\n")
  cat("Gene mapping cache contains", length(ls(gene_mapping_cache)), "entries\n")
}

# Run analysis
process_all_comparisons(deg_results_dir)