# Load necessary libraries for data manipulation, visualization, and analysis
library(Matrix)        # For handling sparse matrices
library(RSpectra)      # For spectral decomposition (e.g., PCA)
library(RcppEigen)     # Interface to the Eigen C++ library for linear algebra
library(SeuratObject)  # Core Seurat data structures
library(sctransform)   # For normalization and variance stabilization
library(uwot)          # For UMAP dimensionality reduction
library(spatstat.core) # For spatial statistics (not directly used here)
library(Seurat)        # For single-cell RNA-seq data analysis
library(ggplot2)       # For data visualization

# Read in the counts matrix from a Matrix Market file (sparse format)
counts <- readMM("C:/Users/student/Desktop/aaliya/TI_STR.scp.matrix.mtx")

# Read in the barcodes and features (genes) associated with the counts matrix
barcodes <- read.table("C:/Users/student/Desktop/aaliya/TI_STR.scp.barcodes.tsv", stringsAsFactors = FALSE)[,1]
features <- read.csv("C:/Users/student/Desktop/aaliya/TI_STR.scp.features.tsv", stringsAsFactors = FALSE, sep = "\t", header = FALSE)

# Assign row names (genes) and column names (cells) to the counts matrix
rownames(counts) <- make.unique(features[,2])  # Ensure unique gene names
colnames(counts) <- barcodes                   # Cell barcodes as column names

# Create a Seurat object with the counts data for downstream analysis
seurat_TI_Str <- CreateSeuratObject(counts, project = "TI_Str")

# Store the initial number of cells for reference
initial_cell_count <- ncol(seurat_TI_Str)

# Quality control: Calculate mitochondrial gene percentage
seurat_TI_Str[["percent.mt"]] <- PercentageFeatureSet(seurat_TI_Str, pattern = "^MT[-\\.]")

# Visualize QC metrics to identify potential outliers
VlnPlot(seurat_TI_Str, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat_TI_Str, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Scatter plots to explore relationships between QC metrics
library(patchwork)  # For combining multiple plots
plot1 <- FeatureScatter(seurat_TI_Str, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_TI_Str, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1  # Display first scatter plot
plot2  # Display second scatter plot

# Define assay and layer names for data manipulation
assay_name <- "RNA"
layer_name <- "data"

# Retrieve the counts layer matrix from the RNA assay
counts_matrix <- seurat_TI_Str[[assay_name]]@layers$counts

# Add a new layer named 'data' with the same values as the 'counts' layer
LayerData(seurat_TI_Str[[assay_name]], layer = layer_name) <- counts_matrix

# List all layers in the RNA assay to confirm the addition
Layers(seurat_TI_Str[[assay_name]])

# Identify highly variable features (genes) to focus on for downstream analysis
seurat_TI_Str <- FindVariableFeatures(seurat_TI_Str, selection.method = "vst", nfeatures = 2000)

# Plot the variable features to visualize their dispersion
VariableFeaturePlot(seurat_TI_Str)

# Highlight the top 20 variable features on the plot
top_features <- head(VariableFeatures(seurat_TI_Str), 20)
plot1 <- VariableFeaturePlot(seurat_TI_Str)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2  # Combine plots to show labels

# Scale the data to normalize expression values across cells
seurat_TI_Str <- ScaleData(seurat_TI_Str, features = rownames(seurat_TI_Str))

# Perform Principal Component Analysis (PCA) on the variable features
seurat_TI_Str <- RunPCA(seurat_TI_Str, features = VariableFeatures(object = seurat_TI_Str))

# Plot an elbow plot to determine the number of principal components to use
ElbowPlot(seurat_TI_Str)

# Read in the metadata file containing additional information about cells
metadata <- read.table("C:/Users/student/Desktop/aaliya/scp_metadata_combined.v2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Set the first row as column names for the metadata dataframe
colnames(metadata) <- metadata[1, ]

# Remove the first row as it's now the header
metadata <- metadata[-1, ]

# Convert specific columns to numeric data type if needed
metadata$n_genes <- as.numeric(metadata$n_genes)
metadata$n_counts <- as.numeric(metadata$n_counts)

# View the cleaned metadata to verify correctness
head(metadata)

# Extract relevant columns from metadata for integration
metadata_cols <- metadata[, c("NAME", "Celltype", "Type", "Chem", "Site", "library_preparation_protocol", "biosample_id")]

# Get the cell names from the Seurat object
cell_names <- colnames(seurat_TI_Str)

# Identify matching cells between the Seurat object and the metadata
matching_cells <- intersect(cell_names, metadata_cols$NAME)

# Subset the metadata to include only matching cells
metadata_to_add <- metadata[metadata$NAME %in% matching_cells, ]

# Set the cell names as row names for the metadata
rownames(metadata_to_add) <- metadata_to_add$NAME

# Remove the 'NAME' column as it's now redundant
metadata_to_add$NAME <- NULL

# Add the metadata to the Seurat object
seurat_TI_Str <- AddMetaData(object = seurat_TI_Str, metadata = metadata_to_add)

# Verify that the metadata has been added correctly
head(seurat_TI_Str@meta.data)

# Perform neighborhood graph construction and clustering
seurat_TI_Str <- FindNeighbors(seurat_TI_Str, dims = 1:10)
seurat_TI_Str <- FindClusters(seurat_TI_Str, resolution = 0.5)

# Run UMAP for dimensionality reduction and visualization
seurat_TI_Str <- RunUMAP(seurat_TI_Str, dims = 1:10)

# Plot UMAP to visualize clusters without cell type labels
DimPlot(seurat_TI_Str, reduction = "umap")

# Load ggplot2 for advanced plotting
library(ggplot2)

# Create a UMAP plot grouped by 'Celltype' with labels
p <- DimPlot(seurat_TI_Str, reduction = "umap", group.by = "Celltype", label = TRUE, repel = TRUE)

# Define custom colors for each cell type for consistency
celltype_colors <- c(
  "Fibroblasts ADAMDEC1" = "#9370DB",    # MediumPurple
  "Fibroblasts SMOC2 PTGIS" = "#66CDAA", # MediumAquamarine
  "Fibroblasts SFRP2 SLPI" = "#87CEFA",  # LightSkyBlue
  "Fibroblasts NPY SLITRK6" = "#4682B4", # SteelBlue
  "Fibroblasts KCNN3 LY6H" = "#4169E1",  # RoyalBlue
  "Endothelial cells CD36" = "#FFB6C1",  # LightPink
  "Endothelial cells DARC" = "#FF7F50",  # Coral
  "Endothelial cells CA4 CD36" = "#FFD700", # Gold
  "Endothelial cells LTC4S SEMA3G" = "#FFA500", # Orange
  "Myofibroblasts HHIP NPNT" = "#98FB98", # PaleGreen
  "Myofibroblasts GREM1 GREM2" = "darkgreen", # DarkGreen
  "Glial cells" = "#FF69B4",              # HotPink
  "Pericytes HIGD1B STEAP4" = "#D2B48C",  # Tan
  "Pericytes RERGL NTRK2" = "#F4A460",    # SandyBrown
  "Activated fibroblasts CCL19 ADAMADEC1" = "#FFDAB9", # PeachPuff
  "Lymphatics" = "#40E0D0"                # Turquoise
)

# Add a title and apply the custom colors to the UMAP plot
p <- p + 
  ggtitle("Terminal Ileum - Stromal cells") +
  scale_color_manual(values = celltype_colors)

# Display the customized UMAP plot
print(p)

# Plot UMAP grouped by 'Type' (e.g., Healthy, Inflamed)
DimPlot(seurat_TI_Str, reduction = "umap", group.by = "Type")

# Set the cell identities to 'Celltype' for further analysis
Idents(seurat_TI_Str) <- 'Celltype'

# Calculate the proportions of each cell type across different conditions
proportions <- prop.table(table(Idents(seurat_TI_Str), seurat_TI_Str$Type), margin = 2)

# Convert the proportions table to a dataframe for plotting
prop_df <- as.data.frame(proportions)

# Save the proportion dataframe to a CSV file
write.csv(prop_df, "C:/Users/student/Desktop/aaliya/latest/prop_df.csv")

# Load the proportion dataframe (if needed)
prop_df <- read.csv("C:/Users/student/Desktop/aaliya/latest/prop_df.csv")

# Rename columns for clarity
colnames(prop_df) <- c("CellType", "Type", "Proportion")

# Plot the proportions using ggplot2 with custom colors
ggplot(prop_df, aes(x = Proportion, y = Type, fill = CellType)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colors) + # Apply custom cell type colors
  theme_minimal() +
  labs(
    x = "Proportion", 
    y = "Type", 
    fill = "Cell Type", 
    title = "Terminal Ileum - Stromal Cells"
  ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) # Adjust y-axis text

# Visualize expression of marker genes using FeaturePlot
FeaturePlot(seurat_TI_Str, features = c("ADM", "CALCRL"), split.by = "Type")

# Additional FeaturePlots for genes of interest
FeaturePlot(seurat_TI_Str, features = c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3"))
FeaturePlot(seurat_TI_Str, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"))
FeaturePlot(seurat_TI_Str, features = c("TNF", "IL6", "ICAM1", "VCAM1", "SELE", "PECAM1"))

# Define the desired order of cell types for consistent plotting
desired_order <- c(
  "Fibroblasts ADAMDEC1", "Fibroblasts SMOC2 PTGIS", "Fibroblasts SFRP2 SLPI", 
  "Fibroblasts NPY SLITRK6", "Fibroblasts KCNN3 LY6H", "Activated fibroblasts CCL19 ADAMADEC1",
  "Endothelial cells CD36", "Endothelial cells DARC", "Endothelial cells CA4 CD36", 
  "Endothelial cells LTC4S SEMA3G", 
  "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2",
  "Glial cells", "Pericytes HIGD1B STEAP4", "Pericytes RERGL NTRK2", 
  "Lymphatics"
)

# Reorder the 'Celltype' factor levels in the metadata
seurat_TI_Str@meta.data$Celltype <- factor(seurat_TI_Str@meta.data$Celltype, levels = desired_order)

# Verify that the factor levels are correctly set
print(levels(seurat_TI_Str@meta.data$Celltype))

# Set cell identities to 'Celltype' again if necessary
Idents(seurat_TI_Str) <- 'Celltype'

# Customize violin plot for a specific gene (e.g., RAMP3)
VlnPlot(seurat_TI_Str, features = "RAMP3", pt.size = 0.1) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Expression of RAMP3 by Cell Types") + 
  xlab("Cell Types") + 
  ylab("Expression Level") +
  scale_x_discrete(limits = desired_order) # Ensure the order is applied

# Create violin plots for genes of interest grouped by 'Type' (condition)
VlnPlot(seurat_TI_Str, features = c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3"), group.by = "Type")
VlnPlot(seurat_TI_Str, features = c("ICAM1", "VCAM1", "SELE", "TNF", "IL6"), group.by = "Type")
VlnPlot(seurat_TI_Str, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"), group.by = "Type")

# Convert the 'Type' column to a factor with the desired order
seurat_TI_Str@meta.data$Type <- factor(seurat_TI_Str@meta.data$Type, levels = c("Heal", "NonI", "Infl"))

# Create a violin plot with 'ACKR3' expression split by 'Type'
VlnPlot(seurat_TI_Str, features = "ACKR3", pt.size = 0.1, split.by = "Type") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Expression Levels of ACKR3 by Cell Type and Diagnosis") +
  xlab("Cell Types") +
  ylab("Expression Level") +
  scale_x_discrete(limits = desired_order) + # Ensure the order is applied
  scale_fill_manual(values = c("Heal" = "yellow", "Infl" = "red", "NonI" = "cyan"))

# Load dplyr for data manipulation
library(dplyr)

# Extract metadata from the Seurat object
metadata2 <- seurat_TI_Str@meta.data

# Group by 'Type' and 'Celltype' and count the number of cells in each group
cell_counts <- metadata2 %>%
  group_by(Type, Celltype) %>%
  summarise(Count = n()) %>%
  arrange(Type, Celltype)

# Print the cell counts
print(cell_counts)

# Save the cell counts to a CSV file
write.csv(cell_counts, file = "C:/Users/student/Desktop/aaliya/latest/cell_counts.csv", row.names = FALSE)

# Plotting cell counts with ggplot2
library(ggplot2)

ggplot(cell_counts, aes(x = Celltype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 3, angle = 90) +
  labs(title = "Terminal Ileum - Stromal Cells",
       x = "Cell Type",
       y = "Total Cell Count",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# List of genes of interest for expression analysis
genes_of_interest <- c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3", "ICAM1", "VCAM1")

# Extract the normalized data matrix for the genes of interest
data_matrix <- GetAssayData(seurat_TI_Str, slot = "data")[genes_of_interest, ]

# Ensure that metadata and data matrix have matching cells
matching_cells <- intersect(rownames(metadata2), colnames(data_matrix))

# Subset and order both metadata and data matrix to match
metadata2 <- metadata2[matching_cells, , drop = FALSE]
data_matrix <- data_matrix[, matching_cells, drop = FALSE]

# Combine the metadata with the expression data
combined_data <- cbind(metadata2, t(data_matrix))

# Define a threshold for considering gene expression as "expressed"
expression_threshold <- 0

# Function to count the number of cells expressing each gene
count_expressing_cells <- function(data, gene) {
  data %>%
    filter(get(gene) > expression_threshold) %>%
    group_by(Type, Celltype) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(Gene = gene) %>%
    select(Gene, Type, Celltype, Count)
}

# Apply the function to all genes of interest and combine the results
expressing_cells_list <- lapply(genes_of_interest, function(gene) count_expressing_cells(combined_data, gene))
expressing_cells_df <- do.call(rbind, expressing_cells_list)

# Print the number of expressing cells for each gene
print(expressing_cells_df)

# Save the expressing cells data to a CSV file
write.csv(expressing_cells_df, file = "C:/Users/student/Desktop/aaliya/expressing_cells_df.csv", row.names = FALSE)

# Merge expressing cells data with total cell counts
combined_results <- left_join(expressing_cells_df, cell_counts, by = c("Type", "Celltype"))

# Rename columns for clarity
combined_results <- combined_results %>%
  rename(
    Expressed_counts = Count.x,
    Total_cell_counts = Count.y
  )

# Print the combined results
print(combined_results)

# Calculate the percentage of cells expressing each gene
combined_results <- combined_results %>%
  mutate(Percentage = (Expressed_counts / Total_cell_counts) * 100)

# Print the updated results with percentages
print(combined_results)

# Save the combined results to a CSV file
write.csv(combined_results, file = "C:/Users/student/Desktop/aaliya/latest/combined_results.csv", row.names = FALSE)

# Subset the Seurat object for different conditions
seurat_healthy <- subset(seurat_TI_Str, subset = Type == "Heal")
seurat_crohns_infl <- subset(seurat_TI_Str, subset = Type == "Infl")
seurat_crohns_noninf <- subset(seurat_TI_Str, subset = Type == "NonI")

# List of genes of interest for dot plots
genes_of_interest <- c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3")

# Create dot plots for each condition
dot_plot_healthy <- DotPlot(seurat_healthy, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Healthy Condition", x = "Genes", y = "Cell Types")

dot_plot_crohns_infl <- DotPlot(seurat_crohns_infl, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Crohn's Disease Inflamed Condition", x = "Genes", y = "Cell Types")

dot_plot_crohns_noninfl <- DotPlot(seurat_crohns_noninf, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Crohn's Disease Non-Inflamed Condition", x = "Genes", y = "Cell Types")

# Extract data from dot plots
dot_plot_healthy_data <- dot_plot_healthy$data
dot_plot_crohns_infl_data <- dot_plot_crohns_infl$data 
dot_plot_crohns_noninfl_data <- dot_plot_crohns_noninfl$data 

# Save dot plot data to CSV files
write.csv(dot_plot_healthy_data, file = "C:/Users/student/Desktop/aaliya/dot_plot_healthy_data1.csv", row.names = FALSE)
write.csv(dot_plot_crohns_infl_data, file = "C:/Users/student/Desktop/aaliya/dot_plot_crohns_infl_data1.csv", row.names = FALSE)
write.csv(dot_plot_crohns_noninfl_data, file = "C:/Users/student/Desktop/aaliya/dot_plot_crohns_noninfl_data1.csv", row.names = FALSE)

# Modify gene names to include the condition for combined plotting
dot_plot_healthy_data$features.plot <- paste(dot_plot_healthy_data$features.plot, "Healthy", sep = "_")
dot_plot_crohns_infl_data$features.plot <- paste(dot_plot_crohns_infl_data$features.plot, "Crohn_Infl", sep = "_")
dot_plot_crohns_noninfl_data$features.plot <- paste(dot_plot_crohns_noninfl_data$features.plot, "Crohn_NonInf", sep = "_")

# Combine the dot plot data into a single dataframe
combined_dot_plot_data <- rbind(dot_plot_healthy_data, dot_plot_crohns_infl_data, dot_plot_crohns_noninfl_data)

# Set the order of the features.plot factor
ordered_features <- unlist(lapply(genes_of_interest, function(gene) {
  paste(gene, c("Healthy", "Crohn_NonInf", "Crohn_Infl"), sep = "_")
}))
combined_dot_plot_data$features.plot <- factor(combined_dot_plot_data$features.plot, levels = ordered_features)

# Define the desired order of cell types for plotting
desired_order <- rev(c(
  "Pericytes RERGL NTRK2", "Pericytes HIGD1B STEAP4",
  "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2",
  "Glial cells", "Fibroblasts SMOC2 PTGIS", "Fibroblasts NPY SLITRK6",
  "Fibroblasts SFRP2 SLPI", "Fibroblasts KCNN3 LY6H",
  "Fibroblasts ADAMDEC1", "Activated fibroblasts CCL19 ADAMADEC1",
  "Lymphatics", "Endothelial cells LTC4S SEMA3G", 
  "Endothelial cells DARC", "Endothelial cells CD36", 
  "Endothelial cells CA4 CD36"
))

# Ensure that the 'id' column matches the desired levels exactly
combined_dot_plot_data$id <- factor(combined_dot_plot_data$id, levels = desired_order)

# Create the combined dot plot
combined_dot_plot <- ggplot(combined_dot_plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = c("lightblue", "blue", "darkblue"), limits = c(min(combined_dot_plot_data$avg.exp.scaled), max(combined_dot_plot_data$avg.exp.scaled))) +
  scale_size_area(max_size = 10, limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Healthy and Crohn's Disease Conditions (Terminal Ileum - Stromal Cells)", x = "Genes", y = "Cell Types")

# Display the combined dot plot
print(combined_dot_plot)

# Save the combined dot plot data to a CSV file
write.csv(combined_dot_plot_data, file = "C:/Users/student/Desktop/aaliya/combined_dot_plot_data1.csv", row.names = FALSE)

# Save the subsetted Seurat objects for each condition
saveRDS(seurat_crohns_infl, file = "/home/aaliya09/Documents/APC_RNASEq/TI_STR/seurat_crohns_infl.rds")
saveRDS(seurat_healthy, file = "/home/aaliya09/Documents/APC_RNASEq/TI_STR/seurat_healthy.rds")
saveRDS(seurat_crohns_noninf, file = "/home/aaliya09/Documents/APC_RNASEq/TI_STR/seurat_crohns_noninf.rds")

# Define the cell types to subset
cell_types <- c(
  "Fibroblasts ADAMDEC1",
  "Endothelial cells CD36",
  "Myofibroblasts HHIP NPNT",
  "Fibroblasts SMOC2 PTGIS",
  "Endothelial cells DARC",
  "Fibroblasts NPY SLITRK6",
  "Myofibroblasts GREM1 GREM2",
  "Endothelial cells CA4 CD36",
  "Glial cells",
  "Fibroblasts SFRP2 SLPI",
  "Endothelial cells LTC4S SEMA3G",
  "Pericytes HIGD1B STEAP4",
  "Activated fibroblasts CCL19 ADAMADEC1",
  "Lymphatics",
  "Fibroblasts KCNN3 LY6H",
  "Pericytes RERGL NTRK2"
)

# Loop through each cell type and create separate Seurat objects
for (cell_type in cell_types) {
  # Create a valid variable name by replacing spaces and special characters
  var_name <- gsub(" ", "_", cell_type)  # Replace spaces with underscores
  var_name <- make.names(var_name)       # Ensure valid R variable name
  
  # Subset the Seurat object for the specific cell type
  assign(var_name, subset(seurat_TI_Str, subset = Celltype == cell_type))
}

# Save each subsetted cell type as an RDS file
saveRDS(Activated_fibroblasts_CCL19_ADAMADEC1, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Activated_fibroblasts_CCL19_ADAMADEC1.rds")
saveRDS(Fibroblasts_ADAMDEC1, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Fibroblasts_ADAMDEC1.rds")
saveRDS(Endothelial_cells_CD36, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Endothelial_cells_CD36.rds")
saveRDS(Myofibroblasts_HHIP_NPNT, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Myofibroblasts_HHIP_NPNT.rds")
saveRDS(Fibroblasts_SMOC2_PTGIS, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Fibroblasts_SMOC2_PTGIS.rds")
saveRDS(Endothelial_cells_DARC, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Endothelial_cells_DARC.rds")
saveRDS(Fibroblasts_NPY_SLITRK6, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Fibroblasts_NPY_SLITRK6.rds")
saveRDS(Myofibroblasts_GREM1_GREM2, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Myofibroblasts_GREM1_GREM2.rds")
saveRDS(Endothelial_cells_CA4_CD36, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Endothelial_cells_CA4_CD36.rds")
saveRDS(Pericytes_HIGD1B_STEAP4, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Pericytes_HIGD1B_STEAP4.rds")
saveRDS(Lymphatics, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Lymphatics.rds")
saveRDS(Fibroblasts_KCNN3_LY6H, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Fibroblasts_KCNN3_LY6H.rds")
saveRDS(Pericytes_RERGL_NTRK2, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Pericytes_RERGL_NTRK2.rds")
saveRDS(Glial_cells, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Glial_cells.rds")
saveRDS(Fibroblasts_SFRP2_SLPI, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Fibroblasts_SFRP2_SLPI.rds")
saveRDS(Endothelial_cells_LTC4S_SEMA3G, file = "/media/aaliya09/Expansion/123112981_scRNASEq/TI_STROMAL/DATA/Endothelial_cells_LTC4S_SEMA3G.rds")
