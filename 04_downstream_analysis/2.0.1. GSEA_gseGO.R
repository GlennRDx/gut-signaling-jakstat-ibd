# Required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)
library(AnnotationDbi)
library(digest)  # Added for hash-based caching

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

# Perform GO enrichment for all ontologies in one batch
perform_go_enrichment_batch <- function(df, n = 50) {
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
  
  # Run all ontologies in batch
  ont_list <- c("BP", "CC", "MF")
  result_df <- data.frame()
  
  # Process each ontology individually (more reliable than compareCluster)
  for (ont in ont_list) {
    tryCatch({
      enrich <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = ont, 
                      minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
      
      if (!is.null(enrich) && nrow(as.data.frame(enrich)) > 0) {
        top_terms <- as.data.frame(enrich) %>%
          arrange(pvalue) %>%
          slice_head(n = n) %>%
          mutate(ont = ont, direction = ifelse(NES > 0, "Up", "Down"))
        
        result_df <- rbind(result_df, top_terms)
      }
    }, error = function(e) {
      cat("Warning: Failed to process ontology", ont, ":", e$message, "\n")
    })
  }
  
  if (nrow(result_df) > 0) {
    return(result_df)
  }
  
  return(NULL)
}

# Function to filter GO terms containing specific gene IDs
filter_go_terms <- function(go_data, target_genes = c("8651", "3716", "3717", "6772")) {
  if (nrow(go_data) == 0) return(go_data)
  
  # Check which rows contain any of the target genes
  contains_target <- sapply(go_data$core_enrichment, function(genes) {
    gene_list <- strsplit(genes, "/")[[1]]
    any(target_genes %in% gene_list)
  })
  
  return(go_data[contains_target, ])
}

# Vectorized clustering for efficiency
cluster_go_terms_fast <- function(go_data, h = 0.9) {
  if (nrow(go_data) == 0) return(data.frame(ID = character(0), cluster = integer(0)))
  
  genes_list <- strsplit(go_data$core_enrichment, "/")
  n <- length(genes_list)
  
  if (n <= 1) return(data.frame(ID = go_data$ID, cluster = 1))
  
  # Vectorized Jaccard similarity calculation
  sim_matrix <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      intersection <- length(intersect(genes_list[[i]], genes_list[[j]]))
      union <- length(union(genes_list[[i]], genes_list[[j]]))
      sim <- if (union == 0) 0 else intersection / union
      sim_matrix[i, j] <- sim_matrix[j, i] <- sim
    }
  }
  diag(sim_matrix) <- 1
  
  # Fast clustering
  dist_matrix <- 1 - sim_matrix
  hc <- hclust(as.dist(dist_matrix), method = "average")
  clusters <- cutree(hc, h = h)
  
  return(data.frame(ID = go_data$ID, cluster = clusters))
}

# Streamlined plotting
plot_go_enrichment_fast <- function(df, title = NULL) {
  if (nrow(df) == 0) return(NULL)
  
  # Pre-calculate all formatting to avoid repeated operations
  df$label <- paste0(
    ifelse(nchar(df$Description) > 50, 
           paste0(substr(df$Description, 1, 50), "..."), 
           df$Description),
    " - ", df$ID, " - p: ", 
    formatC(df$pvalue, format = "e", digits = 2),
    ifelse(df$pvalue < 0.001, " ***", 
           ifelse(df$pvalue < 0.01, " **", 
                  ifelse(df$pvalue < 0.05, " *", "")))
  )
  
  df <- df %>% 
    arrange(cluster, desc(abs(NES))) %>% 
    mutate(rank = row_number())
  
  # Single color palette generation
  color_palette <- scales::hue_pal()(length(unique(df$cluster)))
  
  ggplot(df, aes(x = abs(NES), y = reorder(label, rank))) +
    geom_point(aes(color = factor(cluster), size = setSize), alpha = 0.7) +
    geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold") +
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), strip.text = element_text(size = 9)) +
    labs(x = "Absolute NES", y = "GO Term", color = "Cluster", size = "Gene Count", title = title)
}

# Optimized main workflow
run_go_analysis <- function(deg_file, output_dir, title = NULL) {
  if (!file.exists(deg_file)) return(FALSE)
  
  # Check if GO results already exist
  csv_output_full <- file.path(output_dir, "GO_enrichment_results_full.csv")
  pdf_output_full <- file.path(output_dir, "GO_enrichment_plot_full.pdf")
  csv_output_filt <- file.path(output_dir, "GO_enrichment_results_filt.csv")
  pdf_output_filt <- file.path(output_dir, "GO_enrichment_plot_filt.pdf")
  
  if (file.exists(csv_output_full) || file.exists(pdf_output_full) || 
      file.exists(csv_output_filt) || file.exists(pdf_output_filt)) {
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
    
    # If user chose not to overwrite, skip this analysis
    if (!overwrite_choice) {
      cat("Skipping (results exist):", deg_file, "\n")
      return(TRUE)  # Return TRUE to indicate successful handling
    }
  }
  
  cat("Processing:", deg_file, "\n")
  
  # Single file read with column standardization
  df <- read.csv(deg_file, stringsAsFactors = FALSE)
  if (!"X" %in% colnames(df)) {
    if ("gene" %in% colnames(df)) df$X <- df$gene
    else {
      cat("Warning: No gene column found in", deg_file, "\n")
      return(FALSE)
    }
  }
  if (!"log2FoldChange" %in% colnames(df)) {
    if ("logFC" %in% colnames(df)) df$log2FoldChange <- df$logFC
    else {
      cat("Warning: No log2FoldChange column found in", deg_file, "\n")
      return(FALSE)
    }
  }
  
  # Batch enrichment analysis
  results <- perform_go_enrichment_batch(df)
  if (is.null(results) || nrow(results) == 0) {
    cat("No enrichment results for", deg_file, "\n")
    return(FALSE)
  }
  
  # Fast clustering and plotting for full results
  clusters_full <- cluster_go_terms_fast(results)
  final_results_full <- merge(results, clusters_full, by = "ID", all.x = TRUE)
  plot_full <- plot_go_enrichment_fast(final_results_full, title = title)
  
  # Filter results and process filtered data
  results_filt <- filter_go_terms(results)
  if (nrow(results_filt) > 0) {
    clusters_filt <- cluster_go_terms_fast(results_filt)
    final_results_filt <- merge(results_filt, clusters_filt, by = "ID", all.x = TRUE)
    plot_filt <- plot_go_enrichment_fast(final_results_filt, title = title)
  } else {
    final_results_filt <- data.frame()
    plot_filt <- NULL
  }
  
  # Batch file operations
  write.csv(final_results_full, csv_output_full, row.names = FALSE)
  if (!is.null(plot_full)) {
    ggsave(pdf_output_full, plot_full, width = 16, height = 12, limitsize = FALSE)
  }
  
  write.csv(final_results_filt, csv_output_filt, row.names = FALSE)
  if (!is.null(plot_filt)) {
    ggsave(pdf_output_filt, plot_filt, width = 16, height = 12, limitsize = FALSE)
  }
  
  cat("Completed:", deg_file, "\n")
  return(TRUE)
}

# Updated directory processing for new structure
process_all_comparisons <- function(base_dir) {
  cat("Starting GO enrichment analysis...\n")
  
  comparisons <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  success <- 0
  total <- 0
  
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
          
          for (cell_type in cell_type_dirs) {
            cell_type_path <- file.path(dataset_path, cell_type)
            deg_files <- list.files(cell_type_path, pattern = "_DEG_results\\.csv$", full.names = TRUE)
            
            for (deg_file in deg_files) {
              total <- total + 1
              title <- paste(cell_type, cell_group, tissue_dir, comparison, sep = " - ")
              if (run_go_analysis(deg_file, cell_type_path, title)) success <- success + 1
            }
          }
        }
      }
    }
  }
  
  cat("GO enrichment completed:", success, "of", total, "analyses successful\n")
  cat("Gene mapping cache contains", length(ls(gene_mapping_cache)), "entries\n")
}

# Run analysis
process_all_comparisons(deg_results_dir)