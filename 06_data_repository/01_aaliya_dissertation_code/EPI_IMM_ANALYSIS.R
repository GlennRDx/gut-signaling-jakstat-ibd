# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# Load data: counts, barcodes, and features
counts <- readMM("/home/aaliya09/Documents/APC_RNASEq/TI_EPI/TI_EPI.scp.matrix.mtx")
barcodes <- read.table("/home/aaliya09/Documents/APC_RNASEq/TI_EPI/TI_EPI.scp.barcodes.tsv", stringsAsFactors=F)[,1]
features <- read.csv("/home/aaliya09/Documents/APC_RNASEq/TI_EPI/TI_EPI.scp.features.tsv", stringsAsFactors=F, sep="\t", header=F)

# Assign unique row and column names to the counts matrix
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes

# Create a Seurat object with counts data
seurat_TI_Epi <- CreateSeuratObject(counts, project="TI_Epi")
initial_cell_count <- ncol(seurat_TI_Epi)  # Store the initial number of cells

# Quality control: Calculate the percentage of mitochondrial genes
seurat_TI_Epi[["percent.mt"]] <- PercentageFeatureSet(seurat_TI_Epi, pattern = "^MT[-\\.]")

# Plot violin plots for quality control metrics
VlnPlot(seurat_TI_Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat_TI_Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# Feature scatter plots to check relationships between key metrics
FeatureScatter(seurat_TI_Epi, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_TI_Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Find highly variable features for downstream analysis
seurat_TI_Epi <- FindVariableFeatures(seurat_TI_Epi, selection.method = "vst", nfeatures = 2000)
 
# Visualize the top variable features
top_features <- head(VariableFeatures(seurat_TI_Epi), 20)
plot1 <- VariableFeaturePlot(seurat_TI_Epi)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

# Load metadata and ensure cell names in the Seurat object match the metadata
metadata <- read.table("/home/aaliya09/Documents/APC_RNASEq/scp_metadata_combined.v2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata) <- metadata[1, ]  # Set the first row as column names
metadata <- metadata[-1, ]  # Remove the first row
metadata$n_genes <- as.numeric(metadata$n_genes)
metadata$n_counts <- as.numeric(metadata$n_counts)

# Ensure metadata and Seurat object cell names match
metadata_cols <- metadata[, c("NAME", "Celltype", "Type","Chem","Site", "library_preparation_protocol", "biosample_id")]
matching_cells <- intersect(colnames(seurat_TI_Epi), metadata_cols$NAME)
metadata_to_add <- metadata[metadata$NAME %in% matching_cells, ]
rownames(metadata_to_add) <- metadata_to_add$NAME
metadata_to_add$NAME <- NULL  # Remove the NAME column

# Add metadata to the Seurat object
seurat_TI_Epi <- AddMetaData(object = seurat_TI_Epi, metadata = metadata_to_add)

# Load UMAP coordinates and integrate with Seurat object
UMAP_coord_file_path <- "/home/aaliya09/Documents/APC_RNASEq/TI_EPI/TI_EPI.scp.X_umap.coords.txt"
UMAP_coordinates <- read.table(UMAP_coord_file_path, stringsAsFactors = F, header = T, row.names = 1)
colnames(UMAP_coordinates) <- paste0("UMAP_", 1:ncol(UMAP_coordinates))  # Rename columns to match UMAP format

# Subset Seurat object and UMAP to include only common cells
common_cells <- intersect(rownames(UMAP_coordinates), colnames(seurat_TI_Epi))
seurat_TI_Epi <- subset(seurat_TI_Epi, cells = common_cells)
UMAP_coordinates <- UMAP_coordinates[common_cells, ]
UMAP_coordinates <- as.data.frame(lapply(UMAP_coordinates, as.numeric))

# Add UMAP coordinates to Seurat object
seurat_TI_Epi[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(UMAP_coordinates), key = "UMAP_", assay = DefaultAssay(seurat_TI_Epi))

# Plot UMAP with cell type labels
DimPlot(seurat_TI_Epi, reduction = "umap", group.by = 'Celltype', label = TRUE, raster = FALSE, repel = TRUE)

# Set cell type identities for downstream analysis
Idents(seurat_TI_Epi) <- "Celltype"

# Visualize known markers on UMAP using FeaturePlot
FeaturePlot(seurat_TI_Epi, features = c("ADM", "CALCRL", "RAMP2", "RAMP3"),raster = FALSE)
FeaturePlot(seurat_TI_Epi, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"),raster = FALSE)
FeaturePlot(seurat_TI_Epi, features = c("TNF", "IL6", "SELE"),raster = FALSE)
FeaturePlot(seurat_TI_Epi, features = c("ICAM1", "VCAM1"),raster = FALSE)

# Reorder cell type factor levels for plotting
desired_order <- c("Stem cells OLFM4 LGR5", "Stem cells OLFM4", "Stem cells OLFM4 GSTA1", "Stem cells OLFM4 PCNA",
                   "Epithelial Cycling cells", "Epithelial cells METTL12 MAFB", "Epithelial HBB HBA",
                   "Enterocytes TMIGD1 MEP1A", "Enterocytes BEST4", "Enterocytes TMIGD1 MEP1A GSTA1",
                   "Goblet cells MUC2 TFF1", "Goblet cells MUC2 TFF1-", "Goblet cells SPINK4",
                   "Paneth cells", "Tuft cells", "Enterochromaffin cells", "L cells")

seurat_TI_Epi@meta.data$Celltype <- factor(seurat_TI_Epi@meta.data$Celltype, levels = desired_order)
Idents(seurat_TI_Epi) <- seurat_TI_Epi@meta.data$Celltype

# Create a violin plot for VCAM1 expression across cell types
VlnPlot(seurat_TI_Epi, features = "VCAM1", pt.size = 0.1) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Expression of VCAM1 by Cell Type") + 
  xlab("Cell Types") + 
  ylab("Expression Level")

# Create violin plots for genes of interest between different conditions
VlnPlot(seurat_TI_Epi, features = c("ADM", "CALCRL", "RAMP2", "RAMP3"), group.by = "Type")
VlnPlot(seurat_TI_Epi, features = c("ICAM1", "VCAM1", "PECAM1","SELE"), group.by = "Type")
VlnPlot(seurat_TI_Epi, features = c("ACKR3", "CXCR7", "CRCP", "GPR182"), group.by = "Type")

# Define condition order and plot GPR182 expression split by 'Type'
seurat_TI_Epi@meta.data$Type <- factor(seurat_TI_Epi@meta.data$Type, levels = c("Heal", "NonI", "Infl"))
VlnPlot(seurat_TI_Epi, features = "GPR182", pt.size = 0.1, split.by = "Type") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Expression Levels of GPR182 by Cell Type and Diagnosis") +
  xlab("Cell Types") +
  ylab("Expression Level") +
  scale_fill_manual(values = c("Heal" = "yellow", "NonI" = "cyan", "Infl" = "red"))


# Load the dplyr library for data manipulation
library(dplyr)

# Extract metadata from the Seurat object
metadata2 <- seurat_TI_Epi@meta.data

# Group the metadata by 'Type' (condition) and 'Celltype' and count the number of cells in each group
cell_counts <- metadata2 %>%
  group_by(Type, Celltype) %>%       # Group by Type and Celltype
  summarise(Count = n()) %>%         # Count the number of occurrences in each group
  arrange(Type, Celltype)             # Arrange the results by Type and Celltype

# Print the resulting cell counts
print(cell_counts)

# Write the cell counts to a CSV file for external analysis
write.csv(cell_counts, file = "/home/aaliya09/Documents/APC_RNASEq/TI_EPI/FINAL_CORRECT_ANALYSIS/cell_counts.csv", row.names = FALSE)

# Load ggplot2 for plotting
library(ggplot2)

# Create a bar plot to visualize the cell counts by Type and Celltype
ggplot(cell_counts, aes(x = Celltype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Create dodged bars for each Type
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 3, angle = 90) +  # Add count labels above bars
  labs(title = "Terminal Ileum - Epithelial cells",  # Add title and labels
       x = "Cell Type",
       y = "Total Cell Count",
       fill = "Condition") +
  theme_minimal() +  # Use a minimal theme for the plot
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels for better readability

# Set the cell identities to 'Celltype' for further analysis
Idents(seurat_TI_Epi) <- 'Celltype'

# Calculate the proportions of each cell type across different conditions (Type)
proportions <- prop.table(table(Idents(seurat_TI_Epi), seurat_TI_Epi$Type), margin = 2)  # Calculate proportions by Type

# Convert the proportions table to a dataframe for plotting
prop_df <- as.data.frame(proportions)

# Rename columns for clarity 
colnames(prop_df) <- c("CellType", "Type", "Proportion")

# Step 2: Convert the 'CellType' column to a factor with reversed levels based on desired order
prop_df$CellType <- factor(prop_df$CellType, levels = desired_order)

# Step 3: Define the desired order for the 'Type' variable
desired_type_order <- c("Infl", "NonI", "Heal")

# Convert the 'Type' column to a factor with the specified order
prop_df$Type <- factor(prop_df$Type, levels = desired_type_order)

# Define colors for each cell type for visualization
celltype_colors <- c(
  "Paneth cells" = "#9370DB",                # MediumPurple
  "Goblet cells MUC2 TFF1" = "#66CDAA",      # MediumAquamarine
  "Goblet cells MUC2 TFF1-" = "#87CEFA",     # LightSkyBlue
  "Goblet cells SPINK4" = "#4682B4",         # SteelBlue
  "Enterocytes TMIGD1 MEP1A" = "#4169E1",    # RoyalBlue
  "Enterocytes BEST4" = "#FFB6C1",           # LightPink
  "Epithelial Cycling cells" = "#FF7F50",    # Coral
  "Stem cells OLFM4 LGR5" = "#FFD700",       # Gold
  "Stem cells OLFM4 PCNA" = "#FFA500",       # Orange
  "Stem cells OLFM4" = "#98FB98",            # PaleGreen
  "Tuft cells" = "#006400",                  # DarkGreen
  "Stem cells OLFM4 GSTA1" = "#8A2BE2",      # BlueViolet
  "Epithelial cells METTL12 MAFB" = "#FF6347", # Tomato
  "Epithelial HBB HBA" = "#FF4500",         # OrangeRed
  "Enterocytes TMIGD1 MEP1A GSTA1" = "#2E8B57", # SeaGreen
  "Enterochromaffin cells" = "#D2691E",      # Chocolate
  "L cells" = "#A52A2A"                      # Brown
)

# Create a stacked bar plot to visualize the proportions of cell types across conditions
ggplot(prop_df, aes(x = Proportion, y = Type, fill = CellType)) + 
  geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot to show proportions
  scale_fill_manual(values = celltype_colors) +  # Apply custom colors defined earlier
  theme_minimal() +  # Use a minimal theme for the plot
  labs(
    x = "Proportion", 
    y = "Condition", 
    fill = "Cell Types", 
    title = "Terminal Ileum - Epithelial Cells"
  ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # Adjust y-axis text for better readability
  guides(fill = guide_legend(reverse = TRUE))  # Reverse legend order if needed

# Subset the Seurat object for each condition (healthy and Crohn's disease)
seurat_healthy <- subset(seurat_TI_Epi, subset = Type == "Heal")          # Healthy condition
seurat_crohns_noninf <- subset(seurat_TI_Epi, subset = Type == "NonI")   # Crohn's disease non-inflamed condition
seurat_crohns_infl <- subset(seurat_TI_Epi, subset = Type == "Infl")      # Crohn's disease inflamed condition

# List of genes of interest for analysis
genes_of_interest <- c("ADM", "CALCRL", "RAMP2", "RAMP3", "VCAM1", "ICAM1", "SELE", "ACKR3", "CRCP", "IL6", "TNF")

# Create a dot plot for the healthy condition
dot_plot_healthy <- DotPlot(seurat_healthy, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  labs(title = "Gene Expression in Healthy Condition", x = "Genes", y = "Cell Types")  # Add title and labels

# Create a dot plot for Crohn's disease inflamed condition
dot_plot_crohns_infl <- DotPlot(seurat_crohns_infl, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Gene Expression in Crohn's Disease Inflamed Condition", x = "Genes", y = "Cell Types")  # Add title and labels

# Create a dot plot for Crohn's disease non-inflamed condition
dot_plot_crohns_noninfl <- DotPlot(seurat_crohns_noninf, features = genes_of_interest, group.by = "Celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Gene Expression in Crohn's Disease Non-Inflamed Condition", x = "Genes", y = "Cell Types")  # Add title and labels

# Extract dot plot data for saving to CSV files
dot_plot_healthy_data <- dot_plot_healthy$data
dot_plot_crohns_infl_data <- dot_plot_crohns_infl$data 
dot_plot_crohns_noninfl_data <- dot_plot_crohns_noninfl$data 

# Write the dot plot data to CSV files for external analysis
write.csv(dot_plot_healthy_data, file = "/home/aaliya09/Documents/APC_RNASEq/TI_EPI/FINAL_CORRECT_ANALYSIS/dot_plot_healthy_data.csv", row.names = FALSE)
write.csv(dot_plot_crohns_infl_data, file = "/home/aaliya09/Documents/APC_RNASEq/TI_EPI/FINAL_CORRECT_ANALYSIS/dot_plot_crohns_infl_data.csv", row.names = FALSE)
write.csv(dot_plot_crohns_noninfl_data, file = "/home/aaliya09/Documents/APC_RNASEq/TI_EPI/FINAL_CORRECT_ANALYSIS/dot_plot_crohns_noninfl_data.csv", row.names = FALSE)

# Modify gene names to include condition for combined analysis
dot_plot_healthy_data$features.plot <- paste(dot_plot_healthy_data$features.plot, "Healthy", sep = "_")
dot_plot_crohns_noninfl_data$features.plot <- paste(dot_plot_crohns_noninfl_data$features.plot, "Crohn_NonInf", sep = "_")
dot_plot_crohns_infl_data$features.plot <- paste(dot_plot_crohns_infl_data$features.plot, "Crohn_Infl", sep = "_")

# Combine the dot plot data into a single dataframe for visualization
combined_dot_plot_data <- rbind(dot_plot_healthy_data, dot_plot_crohns_noninfl_data, dot_plot_crohns_infl_data)

# Set the order of the features.plot factor for consistent display
ordered_features <- unlist(lapply(genes_of_interest, function(gene) {
  paste(gene, c("Healthy", "Crohn_NonInf", "Crohn_Infl"), sep = "_")  # Create ordered feature names
}))
combined_dot_plot_data$features.plot <- factor(combined_dot_plot_data$features.plot, levels = ordered_features)  # Set factor levels

# Create the combined dot plot to visualize gene expression across conditions
combined_dot_plot <- ggplot(combined_dot_plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +  # Use point size for percentage expression and color for average scaled expression
  scale_color_gradientn(colors = c("lightblue", "blue", "darkblue"), limits = c(min(combined_dot_plot_data$avg.exp.scaled), max(combined_dot_plot_data$avg.exp.scaled))) +  # Apply a gradient color scale
  scale_size_area(max_size = 10, limits = c(0, 100)) +  # Set size limits for points
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(title = "Gene Expression in Healthy and Crohn's Disease Conditions", x = "Genes", y = "Cell Types")  # Add title and axis labels

# Print the combined dot plot
print(combined_dot_plot)

