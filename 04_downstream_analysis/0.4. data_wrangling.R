sobj <- readRDS('/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/04_li/terminal_ileum_li.rds')

# PatientID
sobj$PatientID <- sobj$Sample.name
sobj$Sample.name <- NULL  # Remove the old column

# CD_Status
sobj$CD_Status <- sobj$inferred.state
sobj$inferred.state <- NULL  # Remove the old column

# Celltypes
sobj$Celltypes <- sobj$annotation_V2
sobj$annotation_V2 <- NULL  # Remove the old column

sobj$CD_Status[sobj$CD_Status == "Active_CD"] <- "ActiveCD"
sobj$CD_Status[sobj$CD_Status == "Inactive_CD"] <- "InactiveCD"
sobj$CD_Status[sobj$CD_Status == "Control"] <- "Control"

# Overwrite the original file
saveRDS(sobj, '/home/glennrossdolan/Documents/gut-signaling-jakstat-ibd/06_data_repository/00_rds_datasets/09_TI_LI.rds')
