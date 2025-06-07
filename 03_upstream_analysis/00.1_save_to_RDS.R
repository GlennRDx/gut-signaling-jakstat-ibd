library(Seurat)

setwd('/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/02_Pediatric_CD/Archive')

count_data <- Read10X(data.dir = 'matrix_files')

metadata <- read.csv('metadata.csv')

rownames(metadata) <- metadata$index

sobj <- CreateSeuratObject(counts = count_data, meta.data = metadata)

sobj

saveRDS(sobj, file = 'norm_pediatric_chron_ileum.rds')