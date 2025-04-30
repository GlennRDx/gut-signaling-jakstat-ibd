# Load necessary libraries
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(ggplot2)
library(Matrix)
library(MAST)


# Define the path to your Seurat object and where to save results
file_path <- "C:/Users/student/Desktop/aaliya/integrated_analysis/harm/harmseurat_TI_Str_integrated_harm1.rds"
output_path <- "C:/Users/student/Desktop/aaliya/integrated_analysis/harm/"

# Load your Seurat object
seurat_TI_Str <- readRDS(file_path)
range(seurat_TI_Str@assays$RNA@layers$data)
# Check metadata columns
print(colnames(seurat_obj@meta.data))

# Ensure the 'Celltype' column is correct and contains cell types
cell_types <- unique(seurat_obj$Celltype)
print(cell_types)

# Loop through each cell type
for (cell_type in cell_types) {
  
  # Subset the Seurat object by the current cell type
  cell_type_obj <- subset(seurat_obj, subset = Celltype == cell_type)
  
  # Check if the subset has cells
  if (ncol(cell_type_obj) == 0) {
    warning(paste("No cells found for cell type:", cell_type))
    next
  }
  
  # Print the number of cells in the subset
  cat("Number of cells in cell type", cell_type, ":", ncol(cell_type_obj), "\n")
  
  # Define the file name and save the subset as an RDS file
  saveRDS(cell_type_obj, file = paste0(output_path, "Seurat_", gsub(" ", "_", cell_type), ".rds"))
  
  # Print confirmation
  cat("Saved subset for cell type", cell_type, "to", paste0(output_path, "Seurat_", gsub(" ", "_", cell_type), ".rds"), "\n")
}



# Define the path to the RDS file
file_path <- "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/Seurat_Endothelial_cells_LTC4S_SEMA3G.rds"

# Read the RDS file
seurat_obj <- readRDS(file_path)

# Convert the metadata column "Type" to a factor 
seurat_obj$Type <- as.factor(seurat_obj$Type)

# Set identities based on metadata column "Type"
Idents(seurat_obj) <- seurat_obj@meta.data$Type

# Get the metadata and expression data
expr_data <- GetAssayData(seurat_obj, slot = "counts")
metadata <- seurat_obj@meta.data

# Define disease groups
disease_groups <- unique(metadata$Type)


filter_genes <- function(expr_data, metadata, disease_groups, min_cells_per_group = 0.1) {
  # List to store genes to keep
  genes_to_keep <- list()
  
  # Loop through each disease group
  for (disease in disease_groups) {
    # Subset expression data for the current disease group
    group_cells <- which(metadata$Type == disease)
    group_expr_data <- expr_data[, group_cells]
    
    # Check if there are any cells in this subset
    if (ncol(group_expr_data) == 0) {
      warning(paste("No cells found for group:", disease))
      next
    }
    
    # Calculate the fraction of cells with non-zero expression for each gene
    gene_expr_fraction <- rowSums(group_expr_data > 0) / ncol(group_expr_data)
    
    # Filter genes that are expressed in at least min_cells_per_group fraction of cells
    filtered_genes <- names(gene_expr_fraction[gene_expr_fraction >= min_cells_per_group])
    
    # Add the filtered genes to the list
    genes_to_keep[[disease]] <- filtered_genes
  }
  
  # Find the union of filtered genes across all disease groups
  genes_to_keep_final <- Reduce(union, genes_to_keep)
  
  # Print the number of genes in the final list
  cat("Total genes to keep after filtering across all groups:", length(genes_to_keep_final), "\n")
  
  return(genes_to_keep_final)
}

# Apply the function to filter genes
genes_to_keep <- filter_genes(expr_data, metadata, disease_groups, min_cells_per_group = 0.1)

# Define the genes you want to check
genes_of_interest <- c("ADM", "CALCRL", "RAMP2", "RAMP3","ACKR3","ICAM1","VCAM1")

# Check if these genes are in the `genes_to_keep` vector
genes_of_interest %in% degs_noninfl_infl

# Subset the Seurat object to keep only the filtered genes
seurat_obj <- subset(seurat_obj, features = genes_to_keep)



# Differential expression using MAST for each cell type
degs_healthy_vs_noninflamed <- FindMarkers(seurat_obj, ident.1 = "NonI", ident.2 = "Heal", test.use = "MAST", min.cells.group = 10)
degs_healthy_vs_inflamed <- FindMarkers(seurat_obj, ident.1 = "Infl", ident.2 = "Heal", test.use = "MAST", min.cells.group = 10)
degs_noninfl_infl <-  FindMarkers(seurat_obj, ident.1 = "NonI", ident.2 = "Infl", test.use = "MAST", min.cells.group = 10)

# Save differential expression results as CSV files
write.csv(degs_healthy_vs_noninflamed, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/degs_healthy_vs_noninflamed.csv", row.names = TRUE)
write.csv(degs_healthy_vs_inflamed, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/degs_healthy_vs_inflamed.csv", row.names = TRUE)
write.csv(degs_noninfl_infl, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/degs_noninfl_infl.csv", row.names = TRUE)


# Filter based on p-value threshold
filtered_degs_noninflamed <- degs_healthy_vs_noninflamed[degs_healthy_vs_noninflamed$p_val_adj <= 0.05, ]
filtered_degs_inflamed <- degs_healthy_vs_inflamed[degs_healthy_vs_inflamed$p_val_adj <= 0.05, ]
filtered_degs_noninfl_infl <- degs_noninfl_infl[degs_noninfl_infl$p_val_adj <= 0.05, ]

# Convert gene names to Entrez IDs
filtered_degs_noninflamed$entrez <- mapIds(org.Hs.eg.db, keys=row.names(filtered_degs_noninflamed), column="ENTREZID", keytype="SYMBOL", multiVals="first")
filtered_degs_inflamed$entrez <- mapIds(org.Hs.eg.db, keys=row.names(filtered_degs_inflamed), column="ENTREZID", keytype="SYMBOL", multiVals="first")
filtered_degs_noninfl_infl$entrez <- mapIds(org.Hs.eg.db, keys=row.names(filtered_degs_noninfl_infl), column="ENTREZID", keytype="SYMBOL", multiVals="first")

# Filter out NA values
filtered_degs_noninflamed <- filtered_degs_noninflamed[!is.na(filtered_degs_noninflamed$entrez), ]
filtered_degs_inflamed <- filtered_degs_inflamed[!is.na(filtered_degs_inflamed$entrez), ]
filtered_degs_noninfl_infl <- filtered_degs_noninfl_infl[!is.na(filtered_degs_noninfl_infl$entrez), ]


# Function to separate upregulated and downregulated genes
separate_genes <- function(degs) {
  upregulated <- degs[degs$avg_log2FC > 0, ]
  downregulated <- degs[degs$avg_log2FC < 0, ]
  return(list(upregulated = upregulated, downregulated = downregulated))
}

# Apply the function to each DEG list
degs_noninflamed_separated <- separate_genes(filtered_degs_noninflamed)
degs_inflamed_separated <- separate_genes(filtered_degs_inflamed)
degs_noninfl_infl_separated <- separate_genes(filtered_degs_noninfl_infl)

# Function to prepare gene lists
get_gene_list <- function(degs) {
  list(
    upregulated = as.character(degs$entrez[degs$avg_log2FC > 0]),
    downregulated = as.character(degs$entrez[degs$avg_log2FC < 0])
  )
}

# Prepare gene lists for each condition
gene_lists_noninflamed <- get_gene_list(filtered_degs_noninflamed)
gene_lists_inflamed <- get_gene_list(filtered_degs_inflamed)
gene_lists_noninfl_infl <- get_gene_list(filtered_degs_noninfl_infl)


library(clusterProfiler)

# Function to perform KEGG enrichment analysis
perform_kegg_analysis <- function(gene_list, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05) {
  kegg_up <- enrichKEGG(gene = gene_list$upregulated, organism = organism, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  kegg_down <- enrichKEGG(gene = gene_list$downregulated, organism = organism, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  return(list(upregulated = kegg_up, downregulated = kegg_down))
}

# Perform KEGG enrichment analysis for each condition
kegg_results_noninflamed <- perform_kegg_analysis(gene_lists_noninflamed)
kegg_results_inflamed <- perform_kegg_analysis(gene_lists_inflamed)
kegg_results_noninfl_infl <- perform_kegg_analysis(gene_lists_noninfl_infl)



# Function to combine upregulated and downregulated KEGG results
combine_kegg_results <- function(up, down) {
  up_df <- as.data.frame(up)
  down_df <- as.data.frame(down)
  
  up_df$regulation <- "Upregulated"
  down_df$regulation <- "Downregulated"
  
  combined_df <- bind_rows(up_df, down_df)
  return(combined_df)
}

# Combine results for each condition
combined_results_noninflamed <- combine_kegg_results(kegg_results_noninflamed$upregulated, kegg_results_noninflamed$downregulated)
combined_results_inflamed <- combine_kegg_results(kegg_results_inflamed$upregulated, kegg_results_inflamed$downregulated)
combined_results_noninfl_infl <- combine_kegg_results(kegg_results_noninfl_infl$upregulated, kegg_results_noninfl_infl$downregulated)

# Save combined results to CSV files
write.csv(combined_results_noninflamed, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/combined_results_noninflamed.csv", row.names = FALSE)
write.csv(combined_results_inflamed, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/combined_results_inflamed.csv", row.names = FALSE)
write.csv(combined_results_noninfl_infl, "D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/LTC4/combined_results_noninfl_infl.csv", row.names = FALSE)


library(ggplot2)

# Function to create a plot for KEGG results without size legend
plot_kegg_results <- function(combined_results, title) {
  ggplot(combined_results, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), color = regulation)) +
    geom_point(aes(size = GeneRatio), alpha = 0.7) +
    coord_flip() +
    scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red")) +
    labs(title = title, x = "Pathway", y = "-log10(p.adjust)", color = "Regulation") +
    guides(size = "none") + # Remove the size legend
    theme_minimal()
}

# Plot results for each condition
plot_noninflamed <- plot_kegg_results(combined_results_noninflamed, "KEGG Pathways - Noninflamed vs Healthy")
plot_inflamed <- plot_kegg_results(combined_results_inflamed, "KEGG Pathways - Inflamed vs Healthy")
plot_noninfl_infl <- plot_kegg_results(combined_results_noninfl_infl, "KEGG Pathways - Noninflamed vs Inflamed")

# Print plots
print(plot_noninflamed)
print(plot_inflamed)
print(plot_noninfl_infl)


#KEGG PATHWAY VIEWER FOR ENDO CA4CD36 NONINFLAMED VS INFLAMED
degs_noninfl_infl <- read.csv("D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/ca4_cd36/degs_noninfl_infl.csv")
filtered_degs_noninfl_infl <- degs_noninfl_infl[degs_noninfl_infl$p_val_adj <= 0.05, ]

# Ensure you're using the correct column ("X") for the gene symbols
filtered_degs_noninfl_infl$entrez <- mapIds(
  org.Hs.eg.db, 
  keys = filtered_degs_noninfl_infl$X,   # Use the "X" column for gene symbols
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

# Create a named vector with Entrez IDs and log fold change for pathview
gene_data <- setNames(filtered_degs_noninfl_infl$avg_log2FC, filtered_degs_noninfl_infl$entrez)


# List of pathways you're interested in
pathways <- c("hsa04010", "hsa04014","hsa04015","hsa04022",
              "hsa04066","hsa04151","hsa04210",
              "hsa04370","hsa04510","hsa04520",
              "hsa04520","hsa04530","hsa04611",
              "hsa04962","hsa05418","hsa04512","hsa04270")  # Add more KEGG pathway IDs here


BiocManager::install("pathview")
library(pathview)
getwd()
setwd("D:/APC_MICROBIOME/aaliya/integrated_analysis/harm/ca4_cd36/pathview")

# Loop through and visualize each pathway
for (pathway in pathways) {
  pathview(gene.data = gene_data,
           pathway.id = pathway,
           species = "hsa",
           out.suffix = paste0("noninflamed_vs_inflamed_", pathway),
           kegg.native = TRUE)
}
















