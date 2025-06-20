# Required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)
library(AnnotationDbi)
library(digest)  # Added for hash-based caching
library(GO.db)   # Added for GO term filtering

# Define the base directory for DEG results
deg_results_dir <- "/home/glennrdx/Documents/gut-signaling-jakstat-ibd/05_results_repository/DGE_Results"

# Cache for gene mappings to avoid repeated database queries
gene_mapping_cache <- new.env()

# Cache for GO term filtering results
go_filter_cache <- new.env()

# Target genes for JAK-STAT pathway filtering
TARGET_GENES <- c("JAK1", "JAK2", "STAT1", "SOCS1")

# Function to get GO terms containing specific genes
get_go_terms_with_genes <- function(target_genes = TARGET_GENES, ontologies = c("BP", "CC", "MF")) {
  cache_key <- paste("go_filter", paste(target_genes, collapse = "_"), paste(ontologies, collapse = "_"), sep = "_")
  
  # Check cache first
  if (exists(cache_key, envir = go_filter_cache)) {
    cat("Using cached GO term filter results...\n")
    return(get(cache_key, envir = go_filter_cache))
  }
  
  cat("Filtering GO terms containing target genes:", paste(target_genes, collapse = ", "), "...\n")
  
  # Convert gene symbols to Entrez IDs
  target_entrez <- mapIds(org.Hs.eg.db, keys = target_genes, column = "ENTREZID", 
                          keytype = "SYMBOL", multiVals = "first")
  target_entrez <- target_entrez[!is.na(target_entrez)]
  
  if (length(target_entrez) == 0) {
    cat("Warning: No Entrez IDs found for target genes\n")
    return(list())
  }
  
  cat("Found Entrez IDs:", paste(target_entrez, collapse = ", "), "\n")
  
  # Get all GO terms for each ontology
  filtered_terms <- list()
  
  for (ont in ontologies) {
    cat("Processing ontology:", ont, "...\n")
    
    # Get all GO annotations for the organism
    go_annotations <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), 
                                            columns = c("ENTREZID", "GO"), keytype = "ENTREZID")
    
    # Filter for specific ontology
    go_terms_ont <- AnnotationDbi::select(GO.db, keys = go_annotations$GO, 
                                          columns = c("GOID", "ONTOLOGY"), keytype = "GOID")
    go_terms_ont <- go_terms_ont[go_terms_ont$ONTOLOGY == ont & !is.na(go_terms_ont$ONTOLOGY), ]
    
    # Merge to get annotations for this ontology
    ont_annotations <- merge(go_annotations, go_terms_ont, by.x = "GO", by.y = "GOID")
    ont_annotations <- ont_annotations[!is.na(ont_annotations$ENTREZID), ]
    
    # Find GO terms containing any of our target genes
    target_go_terms <- ont_annotations[ont_annotations$ENTREZID %in% target_entrez, "GO"]
    target_go_terms <- unique(target_go_terms[!is.na(target_go_terms)])
    
    cat("Found", length(target_go_terms), "GO terms in", ont, "containing target genes\n")
    
    if (length(target_go_terms) > 0) {
      # Get term descriptions
      term_info <- AnnotationDbi::select(GO.db, keys = target_go_terms, 
                                         columns = c("GOID", "TERM"), keytype = "GOID")
      term_info <- term_info[!is.na(term_info$TERM), ]
      
      filtered_terms[[ont]] <- list(
        go_ids = target_go_terms,
        term_info = term_info
      )
      
      # Print some examples
      if (nrow(term_info) > 0) {
        cat("Example terms for", ont, ":\n")
        print(head(term_info, 3))
      }
    }
  }
  
  # Cache the results
  assign(cache_key, filtered_terms, envir = go_filter_cache)
  
  cat("GO term filtering completed. Found terms in", length(filtered_terms), "ontologies\n")
  return(filtered_terms)
}

# Modified function to create custom GO term database
create_custom_go_background <- function(filtered_terms, ont) {
  if (!ont %in% names(filtered_terms) || length(filtered_terms[[ont]]$go_ids) == 0) {
    return(NULL)
  }
  
  return(filtered_terms[[ont]]$go_ids)
}

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

# Function to convert Entrez IDs back to gene symbols
entrez_to_symbol <- function(entrez_ids) {
  entrez_vector <- unlist(strsplit(entrez_ids, "/"))
  symbols <- mapIds(org.Hs.eg.db, keys = entrez_vector, column = "SYMBOL", 
                    keytype = "ENTREZID", multiVals = "first")
  # Remove NAs and return as comma-separated string
  symbols <- symbols[!is.na(symbols)]
  return(paste(symbols, collapse = "/"))
}

# Modified GO enrichment function with pre-filtering
perform_go_enrichment_batch_filtered <- function(df, filtered_terms, n = 50) {
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
  
  # Run all ontologies with pre-filtered terms
  ont_list <- c("BP", "CC", "MF")
  result_df <- data.frame()
  
  # Process each ontology individually with custom background
  for (ont in ont_list) {
    custom_go_terms <- create_custom_go_background(filtered_terms, ont)
    
    if (is.null(custom_go_terms) || length(custom_go_terms) == 0) {
      cat("No filtered GO terms available for ontology:", ont, "\n")
      next
    }
    
    cat("Running GSEA for", ont, "with", length(custom_go_terms), "pre-filtered terms...\n")
    
    tryCatch({
      # Create custom TERM2GENE mapping for filtered terms
      all_annotations <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), 
                                               columns = c("ENTREZID", "GO"), keytype = "ENTREZID")
      
      # Filter annotations to only include our custom GO terms
      custom_annotations <- all_annotations[all_annotations$GO %in% custom_go_terms, ]
      custom_annotations <- custom_annotations[!is.na(custom_annotations$ENTREZID), ]
      
      if (nrow(custom_annotations) == 0) {
        cat("No gene annotations found for filtered GO terms in", ont, "\n")
        next
      }
      
      # Create TERM2GENE data frame
      term2gene <- custom_annotations[, c("GO", "ENTREZID")]
      
      # Get term descriptions
      term2name <- AnnotationDbi::select(GO.db, keys = custom_go_terms, 
                                         columns = c("GOID", "TERM"), keytype = "GOID")
      term2name <- term2name[!is.na(term2name$TERM), c("GOID", "TERM")]
      colnames(term2name) <- c("GO", "Description")
      
      # Run GSEA with custom terms
      enrich <- GSEA(geneList = geneList, 
                     TERM2GENE = term2gene,
                     TERM2NAME = term2name,
                     minGSSize = 5,  # Reduced minimum size for focused analysis
                     maxGSSize = 500, 
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
      
      if (!is.null(enrich) && nrow(as.data.frame(enrich)) > 0) {
        top_terms <- as.data.frame(enrich) %>%
          arrange(pvalue) %>%
          slice_head(n = n) %>%
          mutate(ont = ont, direction = ifelse(NES > 0, "Up", "Down"))
        
        # Add gene symbols column
        if ("core_enrichment" %in% colnames(top_terms)) {
          cat("Converting Entrez IDs to gene symbols for", ont, "...\n")
          top_terms$gene_symbols <- sapply(top_terms$core_enrichment, entrez_to_symbol)
        }
        
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

# Enhanced plotting with JAK-STAT focus
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
  
  p <- ggplot(df, aes(x = abs(NES), y = reorder(label, rank))) +
    geom_point(aes(color = factor(cluster), size = setSize), alpha = 0.7) +
    geom_text(aes(label = cluster), size = 2.5, color = "white", fontface = "bold") +
    facet_grid(ont ~ direction, scales = "free_y", space = "free_y") +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), 
          strip.text = element_text(size = 9),
          plot.title = element_text(size = 12, face = "bold")) +
    labs(x = "Absolute NES", 
         y = "GO Term", 
         color = "Cluster", 
         size = "Gene Count", 
         title = paste("JAK-STAT Focused GO Enrichment:", title),
         subtitle = paste("Filtered for GO terms containing:", paste(TARGET_GENES, collapse = ", ")))
  
  return(p)
}

# Modified main workflow with pre-filtering
run_go_analysis_filtered <- function(deg_file, output_dir, filtered_terms, title = NULL) {
  if (!file.exists(deg_file)) return(FALSE)
  
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
  
  # Batch enrichment analysis with pre-filtering
  results <- perform_go_enrichment_batch_filtered(df, filtered_terms)
  if (is.null(results) || nrow(results) == 0) {
    cat("No enrichment results for", deg_file, "\n")
    return(FALSE)
  }
  
  # Fast clustering and plotting
  clusters <- cluster_go_terms_fast(results)
  final_results <- merge(results, clusters, by = "ID", all.x = TRUE)
  plot <- plot_go_enrichment_fast(final_results, title = title)
  
  # Batch file operations
  write.csv(final_results, file.path(output_dir, "GO_enrichment_results_JAK_STAT_filtered.csv"), row.names = FALSE)
  if (!is.null(plot)) {
    ggsave(file.path(output_dir, "GO_enrichment_plot_JAK_STAT_filtered.pdf"), plot, width = 16, height = 12, limitsize = FALSE)
  }
  
  cat("Completed:", deg_file, "\n")
  return(TRUE)
}

# Enhanced directory processing with pre-filtering
process_all_comparisons_filtered <- function(base_dir) {
  cat("Starting JAK-STAT focused GO enrichment analysis...\n")
  cat("Target genes:", paste(TARGET_GENES, collapse = ", "), "\n")
  
  # Pre-filter GO terms (this happens once for all analyses)
  filtered_terms <- get_go_terms_with_genes()
  
  if (length(filtered_terms) == 0) {
    cat("Error: No GO terms found containing target genes. Exiting.\n")
    return()
  }
  
  # Print summary of filtered terms
  total_terms <- sum(sapply(filtered_terms, function(x) length(x$go_ids)))
  cat("Total filtered GO terms across all ontologies:", total_terms, "\n")
  for (ont in names(filtered_terms)) {
    cat("  ", ont, ":", length(filtered_terms[[ont]]$go_ids), "terms\n")
  }
  
  comparisons <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  success <- 0
  total <- 0
  
  for (comparison in comparisons) {
    comparison_dir <- file.path(base_dir, comparison)
    tissue_dirs <- list.dirs(comparison_dir, recursive = FALSE, full.names = FALSE)
    
    for (tissue_dir in tissue_dirs) {
      tissue_path <- file.path(comparison_dir, tissue_dir)
      cell_type_dirs <- list.dirs(tissue_path, recursive = FALSE, full.names = FALSE)
      
      if (length(cell_type_dirs) == 0) {
        deg_file <- file.path(tissue_path, "DEG_results.csv")
        if (file.exists(deg_file)) {
          total <- total + 1
          title <- paste(tissue_dir, comparison, sep = " - ")
          if (run_go_analysis_filtered(deg_file, tissue_path, filtered_terms, title)) success <- success + 1
        }
      } else {
        for (cell_type in cell_type_dirs) {
          cell_type_path <- file.path(tissue_path, cell_type)
          deg_file <- file.path(cell_type_path, "DEG_results.csv")
          if (file.exists(deg_file)) {
            total <- total + 1
            title <- paste(cell_type, tissue_dir, comparison, sep = " - ")
            if (run_go_analysis_filtered(deg_file, cell_type_path, filtered_terms, title)) success <- success + 1
          }
        }
      }
    }
  }
  
  cat("JAK-STAT focused GO enrichment completed:", success, "of", total, "analyses successful\n")
  cat("Gene mapping cache contains", length(ls(gene_mapping_cache)), "entries\n")
  cat("GO filter cache contains", length(ls(go_filter_cache)), "entries\n")
}

# Run the enhanced analysis
cat("=== JAK-STAT Focused GO Enrichment Analysis ===\n")
cat("This analysis will only test GO terms that contain JAK1, JAK2, STAT1, or SOCS1\n\n")

process_all_comparisons_filtered(deg_results_dir)