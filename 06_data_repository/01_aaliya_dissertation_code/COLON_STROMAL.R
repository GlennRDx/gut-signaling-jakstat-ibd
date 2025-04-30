# Load necessary libraries
library(Matrix)
library(RSpectra)
library(RcppEigen)
library(SeuratObject)
library(sctransform)
library(uwot)
library(spatstat.core)
library(Seurat)
library(ggplot2)

# Read the matrix market file (sparse matrix format)
counts <- readMM("D:/APC_MICROBIOME/Colon_Str/CO_STR.scp.matrix.mtx")

# Read the barcode and feature files
barcodes <- read.table("D:/APC_MICROBIOME/Colon_Str/CO_STR.scp.barcodes.tsv", stringsAsFactors = FALSE)[, 1]
features <- read.csv("D:/APC_MICROBIOME/Colon_Str/CO_STR.scp.features.tsv", stringsAsFactors = FALSE, sep = "\t", header = FALSE)

# Assign row and column names to the counts matrix
rownames(counts) <- make.unique(features[, 2])  # Make sure gene names are unique
colnames(counts) <- barcodes                    # Assign barcodes (cell IDs) to column names

# Create a Seurat object from the counts matrix
seurat_Colon_Str <- CreateSeuratObject(counts, project = "Col_Str")

# Store the initial number of cells in the Seurat object for reference
initial_cell_count <- ncol(seurat_Colon_Str)

# Add mitochondrial percentage data for quality control
seurat_Colon_Str[["percent.mt"]] <- PercentageFeatureSet(seurat_Colon_Str, pattern = "^MT[-\\.]")

# Plot violin plots for quality control metrics (number of features, number of counts, and mitochondrial percentage)
VlnPlot(seurat_Colon_Str, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Plot without point overlay
VlnPlot(seurat_Colon_Str, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Scatter plots for quality control visualization
library(patchwork)
plot1 <- FeatureScatter(seurat_Colon_Str, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_Colon_Str, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1  # Display plot1
plot2  # Display plot2

# Define assay and layer names for data manipulation
assay_name <- "RNA"
layer_name <- "data"

# Retrieve the counts matrix from the RNA assay
counts_matrix <- seurat_Colon_Str[[assay_name]]@layers$counts

# Add a new layer named 'data' to the Seurat object with the same values as the 'counts' layer
LayerData(seurat_Colon_Str[[assay_name]], layer = layer_name) <- counts_matrix

# List all layers to confirm the addition of the 'data' layer
Layers(seurat_Colon_Str[[assay_name]])

# Find the top 2000 variable features in the data using 'vst' method
seurat_Colon_Str <- FindVariableFeatures(seurat_Colon_Str, selection.method = "vst", nfeatures = 2000)

# Plot the variable features (genes)
VariableFeaturePlot(seurat_Colon_Str)

# Highlight the top 20 variable features
top_features <- head(VariableFeatures(seurat_Colon_Str), 20)
plot1 <- VariableFeaturePlot(seurat_Colon_Str)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2  # Combine both plots

# Scale the data (normalization)
seurat_Colon_Str <- ScaleData(seurat_Colon_Str, features = rownames(seurat_Colon_Str))

# Perform Principal Component Analysis (PCA) on the variable features
seurat_Colon_Str <- RunPCA(seurat_Colon_Str, features = VariableFeatures(object = seurat_Colon_Str))

# Plot the elbow plot to determine the number of PCs to use
ElbowPlot(seurat_Colon_Str)

# Read the metadata file containing additional information about the cells
metadata <- read.table("D:/APC_MICROBIOME/scp_metadata_combined.v2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Set the first row as column names and remove it from the data
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,]

# Convert numeric columns in the metadata to numeric format
metadata$n_genes <- as.numeric(metadata$n_genes)
metadata$n_counts <- as.numeric(metadata$n_counts)

# Preview the cleaned metadata
head(metadata)

# Match and integrate metadata with the Seurat object
metadata_cols <- metadata[, c("NAME", "Celltype", "Type", "Chem", "Site", "library_preparation_protocol", "biosample_id")]
cell_names <- colnames(seurat_Colon_Str)
matching_cells <- intersect(cell_names, metadata_cols$NAME)
metadata_to_add <- metadata[metadata$NAME %in% matching_cells,]
rownames(metadata_to_add) <- metadata_to_add$NAME
metadata_to_add$NAME <- NULL

# Add metadata to the Seurat object
seurat_Colon_Str <- AddMetaData(object = seurat_Colon_Str, metadata = metadata_to_add)

# Verify that metadata has been added
head(seurat_Colon_Str@meta.data)

# Perform neighborhood graph and clustering
seurat_Colon_Str <- FindNeighbors(seurat_Colon_Str, dims = 1:10)
seurat_Colon_Str <- FindClusters(seurat_Colon_Str, resolution = 0.5)

# Run UMAP for visualization
seurat_Colon_Str <- RunUMAP(seurat_Colon_Str, dims = 1:10)

# Plot UMAP of the clusters
DimPlot(seurat_Colon_Str, reduction = "umap")

# Define custom colors for the cell types
stromal_cell_colors <- c(
  "Fibroblasts ADAMDEC1" = "#9370DB",   # MediumPurple
  "Endothelial cells CD36" = "#FFB6C1", # LightPink
  "Fibroblasts KCNN3 LY6H" = "#4169E1", # RoyalBlue
  "Lymphatics" = "#40E0D0",             # Turquoise
  "Endothelial cells LTC4S SEMA3G" = "#FFA500", # Orange
  "Endothelial cells DARC" = "#FF7F50", # Coral
  "Pericytes HIGD1B STEAP4" = "#D2B48C",# Tan
  "Glial cells" = "#FF69B4",            # HotPink
  "Myofibroblasts HHIP NPNT" = "#98FB98", # PaleGreen
  "Stromal Cycling cells" = "#B0C4DE",  # LightSteelBlue
  "Inflammatory fibroblasts IL11 CHI3L1" = "#A52A2A", # Brown
  "Fibroblasts NPY SLITRK6" = "#4682B4", # SteelBlue
  "Fibroblasts SFRP2 SLPI" = "#87CEFA",  # LightSkyBlue
  "Fibroblasts SMOC2 PTGIS" = "#66CDAA", # MediumAquamarine
  "Pericytes RERGL NTRK2" = "#F4A460",   # SandyBrown
  "Activated fibroblasts CCL19 ADAMADEC1" = "#008080", # teal
  "Myofibroblasts GREM1 GREM2" = "darkgreen" # DarkGreen
)

# Verify unique cell types
unique(seurat_Colon_Str$Celltype)

# Plot UMAP with corrected color mapping
p <- DimPlot(
  seurat_Colon_Str,
  reduction = "umap",
  group.by = "Celltype",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + 
  scale_color_manual(values = stromal_cell_colors) + 
  ggtitle("Colon - Stromal Cells")

# Print the UMAP plot
print(p)

# Set cell identities to 'Celltype'
Idents(seurat_Colon_Str) <- "Celltype"

# Visualize marker genes on UMAP
FeaturePlot(seurat_Colon_Str, features = c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3"))
FeaturePlot(seurat_Colon_Str, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"))
FeaturePlot(seurat_Colon_Str, features = c("TNF", "IL6", "ICAM1", "VCAM1", "SELE", "PECAM1"))

# Define a desired order for cell types and reorder the 'Celltype' factor
desired_order <- c(
  "Fibroblasts ADAMDEC1", "Fibroblasts KCNN3 LY6H", "Fibroblasts NPY SLITRK6", 
  "Fibroblasts SFRP2 SLPI", "Fibroblasts SMOC2 PTGIS", "Inflammatory fibroblasts IL11 CHI3L1", 
  "Activated fibroblasts CCL19 ADAMADEC1", "Endothelial cells CD36", "Endothelial cells LTC4S SEMA3G", 
  "Endothelial cells DARC", "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2", 
  "Pericytes HIGD1B STEAP4", "Pericytes RERGL NTRK2", "Lymphatics", "Glial cells", "Stromal Cycling cells"
)

# Reorder the Celltype factor levels in the metadata
# Ensures the Celltype column has a specified order for better plotting
seurat_Colon_Str@meta.data$Celltype <- factor(seurat_Colon_Str@meta.data$Celltype, levels = desired_order)

# Verify that the factor levels for Celltype are correctly set in the metadata
print(levels(seurat_Colon_Str@meta.data$Celltype))

# Customize the violin plot for ACKR3 gene expression by Celltype
VlnPlot(seurat_Colon_Str, features = "ACKR3", pt.size = 0.1) + 
  theme(legend.position = "none",                               # Remove the legend
        axis.text.x = element_text(angle = 45, hjust = 1)) +    # Rotate the x-axis text for better visibility
  ggtitle("Expression of ACKR3 by Cell Types") +                # Add plot title
  xlab("Cell Types") +                                          # X-axis label
  ylab("Expression Level") +                                    # Y-axis label
  scale_x_discrete(limits = desired_order)                      # Ensure the cell types are in the correct order

# Violin plots for a general set of genes across different conditions (Type variable)
VlnPlot(seurat_Colon_Str, features = c("ADM", "CALCRL", "RAMP2", "RAMP3", "ACKR3"), group.by = "Type")

# Violin plots for another set of genes between conditions (grouped by Type)
VlnPlot(seurat_Colon_Str, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"), group.by = "Type")

# Violin plots for inflammation-related genes between two conditions (grouped by Type)
VlnPlot(seurat_Colon_Str, features = c("TNF", "IL6", "ICAM1", "VCAM1", "SELE", "PECAM1"), group.by = "Type")

# Reorder the 'Type' factor levels (conditions: Heal, NonI, Infl)
seurat_Colon_Str@meta.data$Type <- factor(seurat_Colon_Str@meta.data$Type, levels = c("Heal", "NonI", "Infl"))

# Violin plot with ACKR3 expression split by condition (Type)
VlnPlot(seurat_Colon_Str, features = "ACKR3", pt.size = 0.1, split.by = "Type") +
  theme(legend.position = "top",                                # Move legend to the top
        axis.text.x = element_text(angle = 45, hjust = 1)) +    # Rotate the x-axis text
  ggtitle("Expression Levels of ACKR3 by Cell Type and Diagnosis") + # Add plot title
  xlab("Cell Types") +                                          # X-axis label
  ylab("Expression Level") +                                    # Y-axis label
  scale_x_discrete(limits = desired_order) +                    # Ensure cell types are in the correct order
  scale_fill_manual(values = c("Heal" = "yellow",               # Set custom colors for each condition
                               "Infl" = "red", 
                               "NonI" = "cyan"))


# Load the necessary library for data manipulation
library(dplyr)

# Extract metadata from the Seurat object for further analysis
metadata2 <- seurat_Colon_Str@meta.data

# Group data by 'Type' (condition) and 'Celltype', and count the number of cells in each group
cell_counts <- metadata2 %>%
  group_by(Type, Celltype) %>%
  summarise(Count = n()) %>%
  arrange(Type, Celltype)                                      # Sort by condition and cell type

# Print the grouped cell count data
print(cell_counts)

# Save the cell count data to a CSV file for external analysis
write.csv(cell_counts, file = "D:/APC_MICROBIOME/Colon_Str/latest/cell_counts.csv", row.names = FALSE)

# Load cell counts from CSV file
cell_counts <- read.csv("D:/APC_MICROBIOME/Colon_Str/cell_counts.csv")

# Define the desired order of cell types for consistent plotting
desired_order <- c("Pericytes RERGL NTRK2", "Pericytes HIGD1B STEAP4",
                   "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2",
                   "Glial cells", "Fibroblasts SMOC2 PTGIS", "Fibroblasts NPY SLITRK6",
                   "Fibroblasts SFRP2 SLPI", "Fibroblasts KCNN3 LY6H",
                   "Fibroblasts ADAMDEC1", "Activated fibroblasts CCL19 ADAMADEC1",
                   "Lymphatics", "Endothelial cells LTC4S SEMA3G", 
                   "Endothelial cells DARC", "Endothelial cells CD36", 
                   "Stromal Cycling cells","Inflammatory fibroblasts IL11 CHI3L1")

# Convert the 'Celltype' column to a factor with the desired order
cell_counts$Celltype <- factor(cell_counts$Celltype, levels = desired_order)

# Define the desired order for conditions (Type variable)
desired_type_order <- c("Heal", "NonI", "Infl")

# Convert the 'Type' column to a factor with the specified order
cell_counts$Type <- factor(cell_counts$Type, levels = desired_type_order)

# Load the ggplot2 library for data visualization
library(ggplot2)

# Plotting cell count data across cell types and conditions with bar plot
ggplot(cell_counts, aes(x = Celltype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +            # Bar plot with separate bars for each condition
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 3, angle = 90) + # Add cell count labels
  labs(title = "Colon - Stromal Cells",                        # Add plot title
       x = "Cell Type",                                        # X-axis label
       y = "Total Cell Count",                                 # Y-axis label
       fill = "Condition") +                                   # Legend title
  theme_minimal() +                                            # Minimal theme for a clean plot
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # Rotate x-axis labels for better readability

# List of genes of interest for expression analysis
genes_of_interest <- c("ADM", "CALCRL","RAMP2","RAMP3","ACKR3","TNF", "IL6", "ICAM1", "VCAM1", "SELE")

# Extract the normalized gene expression data from the Seurat object
data_matrix <- GetAssayData(seurat_Colon_Str, slot = "data")[genes_of_interest, ]

# Ensure the row names of metadata match the column names of the data matrix
matching_cells <- intersect(rownames(metadata2), colnames(data_matrix))

# Subset and order both metadata and data matrix to ensure they match
metadata2 <- metadata2[matching_cells, , drop = FALSE]
data_matrix <- data_matrix[, matching_cells, drop = FALSE]

# Combine metadata with gene expression data for further analysis
combined_data <- cbind(metadata2, t(data_matrix))

# Define a threshold for gene expression (e.g., any expression level greater than 0)
expression_threshold <- 0

# Function to count the number of cells expressing each gene of interest
count_expressing_cells <- function(data, gene) {
  data %>%
    filter(get(gene) > expression_threshold) %>%                 # Filter cells expressing the gene above the threshold
    group_by(Type, Celltype) %>%                                 # Group by condition and cell type
    summarise(Count = n(), .groups = 'drop') %>%                 # Count cells expressing the gene
    mutate(Gene = gene) %>%                                      # Add gene name for reference
    select(Gene, Type, Celltype, Count)                          # Select relevant columns
}

# Apply the counting function to each gene and combine the results into a single dataframe
expressing_cells_list <- lapply(genes_of_interest, function(gene) count_expressing_cells(combined_data, gene))
expressing_cells_df <- do.call(rbind, expressing_cells_list)

# Print the number of expressing cells per gene, type, and cell type
print(expressing_cells_df)

# Join the expressing cell counts with total cell counts for each condition and cell type
combined_results <- left_join(expressing_cells_df, cell_counts, by = c("Type", "Celltype"))

# Rename columns to better describe the data
combined_results <- combined_results %>%
  rename(
    Expressed_counts = Count.x,      # Number of cells expressing the gene
    Total_cell_counts = Count.y      # Total number of cells
  )

# Print the combined results for validation
print(combined_results)

# Calculate the percentage of cells expressing each gene per condition and cell type
combined_results <- combined_results %>%
  mutate(Percentage = (Expressed_counts / Total_cell_counts) * 100)

# Print the final dataframe with percentages
print(combined_results)

# Save the final results as a CSV file
write.csv(combined_results, file = "D:/APC_MICROBIOME/Colon_Str/latest/combined_results.csv", row.names = FALSE)



# Set the identity of the Seurat object to "Celltype" for proportion calculation
Idents(seurat_Colon_Str) <- "Celltype"

# Calculate the proportions of each cell type across different conditions (Heal, NonI, Infl)
proportion_data <- prop.table(table(Idents(seurat_Colon_Str), seurat_Colon_Str$Type), margin = 2)

# Convert the proportion data to a dataframe for further plotting
prop_df <- as.data.frame(proportion_data)

# Save the proportion dataframe to a CSV file
write.csv(prop_df, 'D:/APC_MICROBIOME/Colon_Str/latest/colon_str_prop_df.csv')

# Load the proportion dataframe from the CSV file
prop_df <- read.csv("D:/APC_MICROBIOME/Colon_Str/latest/colon_str_prop_df.csv")

# Rename columns to match the cell type and condition names
colnames(prop_df) <- c("CellType", "Type", "Proportion")

# Step 1: Define the desired order of cell types for visualization (reversed for better plotting)
desired_order <- rev(c("Pericytes RERGL NTRK2", "Pericytes HIGD1B STEAP4",
                       "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2",
                       "Glial cells", "Fibroblasts SMOC2 PTGIS", "Fibroblasts NPY SLITRK6",
                       "Fibroblasts SFRP2 SLPI", "Fibroblasts KCNN3 LY6H",
                       "Fibroblasts ADAMDEC1", "Activated fibroblasts CCL19 ADAMDEC1",
                       "Lymphatics", "Endothelial cells LTC4S SEMA3G", 
                       "Endothelial cells DARC", "Endothelial cells CD36", 
                       "Stromal Cycling cells","Inflammatory fibroblasts IL11 CHI3L1"))

# Step 2: Convert the 'CellType' column to a factor with the specified order
prop_df$CellType <- factor(prop_df$CellType, levels = desired_order)

# Step 1: Define the desired order of conditions (Type variable)
desired_type_order <- c("Infl", "NonI", "Heal")

# Step 2: Convert the 'Type' column to a factor with the specified order
prop_df$Type <- factor(prop_df$Type, levels = desired_type_order)

# Define custom colors for each cell type
celltype_colors <- c(
  "Pericytes RERGL NTRK2" = "#F4A460",                  # SandyBrown
  "Pericytes HIGD1B STEAP4" = "#D2B48C",                 # Tan
  "Myofibroblasts HHIP NPNT" = "#98FB98",                # PaleGreen
  "Myofibroblasts GREM1 GREM2" = "#006400",              # DarkGreen
  "Glial cells" = "#FF69B4",                             # HotPink
  "Fibroblasts SMOC2 PTGIS" = "#66CDAA",                 # MediumAquamarine
  "Fibroblasts NPY SLITRK6" = "#4682B4",                 # SteelBlue
  "Fibroblasts SFRP2 SLPI" = "#87CEFA",                  # LightSkyBlue
  "Fibroblasts KCNN3 LY6H" = "#4169E1",                  # RoyalBlue
  "Fibroblasts ADAMDEC1" = "#9370DB",                    # MediumPurple
  "Activated fibroblasts CCL19 ADAMDEC1" = "#FFFF00",    # Yellow
  "Lymphatics" = "#40E0D0",                              # Turquoise
  "Endothelial cells LTC4S SEMA3G" = "#FFA500",          # Orange
  "Endothelial cells DARC" = "#FF7F50",                  # Coral
  "Endothelial cells CD36" = "#FFB6C1",                  # LightPink
  "Stromal Cycling cells" = "#B0C4DE",                   # LightSteelBlue
  "Inflammatory fibroblasts IL11 CHI3L1" = "#A52A2A"     # Brown
)

# Plot the proportions of cell types across conditions using ggplot2
ggplot(prop_df, aes(x = Proportion, y = Type, fill = CellType)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = celltype_colors) +  # Apply custom colors
  theme_minimal() +
  labs(
    x = "Proportion", 
    y = "Condition", 
    fill = "Cell Types", 
    title = "Colon - Stromal Cells"
  ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Customize the y-axis text
  guides(fill = guide_legend(reverse = TRUE))  # Reverse legend order if needed


# Subset Seurat object for different conditions (Healthy, Crohn's Inflamed, Crohn's Non-Inflamed)
seurat_healthy <- subset(seurat_Colon_Str, subset = Type == "Heal")
seurat_crohns_infl <- subset(seurat_Colon_Str, subset = Type == "Infl")
seurat_crohns_noninf <- subset(seurat_Colon_Str, subset = Type == "NonI")

# List of genes of interest for dot plot
genes_of_interest <- c("ADM", "CALCRL","RAMP2","RAMP3","ACKR3")

# Create dot plot for healthy condition
dot_plot_healthy <- DotPlot(seurat_healthy, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Healthy Condition", x = "Genes", y = "Cell Types")

# Create dot plot for Crohn's disease inflamed condition
dot_plot_crohns_infl <- DotPlot(seurat_crohns_infl, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Crohn's Disease Inflamed Condition", x = "Genes", y = "Cell Types")

# Create dot plot for Crohn's disease non-inflamed condition
dot_plot_crohns_noninfl <- DotPlot(seurat_crohns_noninf, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Crohn's Disease Non-Inflamed Condition", x = "Genes", y = "Cell Types")


# Extract data from the dot plots for exporting and further analysis
dot_plot_healthy_data <- dot_plot_healthy$data
dot_plot_crohns_infl_data <- dot_plot_crohns_infl$data 
dot_plot_crohns_noninfl_data <- dot_plot_crohns_noninfl$data 

# Write the extracted dot plot data to CSV files for each condition
write.csv(dot_plot_healthy_data, file = "/home/aaliya09/Documents/APC_RNASEq/CO_IMMUNE/dot_plot_healthy_data.csv", row.names = FALSE)
write.csv(dot_plot_crohns_infl_data, file = "/home/aaliya09/Documents/APC_RNASEq/CO_IMMUNE/dot_plot_crohns_infl_data.csv", row.names = FALSE)
write.csv(dot_plot_crohns_noninfl_data, file = "/home/aaliya09/Documents/APC_RNASEq/CO_IMMUNE/dot_plot_crohns_noninfl_data.csv", row.names = FALSE)

# Modify the gene names to include the condition (for combined plotting)
dot_plot_healthy_data$features.plot <- paste(dot_plot_healthy_data$features.plot, "Healthy", sep = "_")
dot_plot_crohns_infl_data$features.plot <- paste(dot_plot_crohns_infl_data$features.plot, "Crohn_Infl", sep = "_")
dot_plot_crohns_noninfl_data$features.plot <- paste(dot_plot_crohns_noninfl_data$features.plot, "Crohn_NonInf", sep = "_")

# Combine all the dot plot data into a single dataframe for comprehensive visualization
combined_dot_plot_data <- rbind(dot_plot_healthy_data, dot_plot_crohns_infl_data, dot_plot_crohns_noninfl_data)

# Set the order of the features.plot factor for proper ordering in the plot
ordered_features <- unlist(lapply(genes_of_interest, function(gene) {
  paste(gene, c("Healthy", "Crohn_NonInf", "Crohn_Infl"), sep = "_")
}))
combined_dot_plot_data$features.plot <- factor(combined_dot_plot_data$features.plot, levels = ordered_features)

# Ensure the 'id' column matches the desired levels for cell types
combined_dot_plot_data$id <- factor(combined_dot_plot_data$id, levels = desired_order)


# Create the combined dot plot with size and color representing expression levels and percentage
combined_dot_plot <- ggplot(combined_dot_plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradientn(colors = c("lightblue", "blue", "darkblue"), 
                        limits = c(min(combined_dot_plot_data$avg.exp.scaled), max(combined_dot_plot_data$avg.exp.scaled))) +
  scale_size_area(max_size = 10, limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression in Healthy and Crohn's Disease Conditions (Colon - Stromal)", x = "Genes", y = "Cell Types")

# Display the combined dot plot
print(combined_dot_plot)

# Save the combined dot plot data to a CSV file
write.csv(combined_dot_plot_data, file = "D:/APC_MICROBIOME/Colon_Str/latest/combined_dot_plot_data.csv", row.names = FALSE)
