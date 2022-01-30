options(stringsAsFactors = FALSE)

Path_Code <- "./01_Bulk_RNASeq/01_Code/"
Path_Data <- "./01_Bulk_RNASeq/02_Bulk_RNASeq_Data/"

setwd(Path_Data)

source(paste0(Path_Code, "01_Normal_And_Tumor_boxplot.R"))
source(paste0(Path_Code, "02_Survival.R"))
source(paste0(Path_Code, "03_Clinical.R"))
source(paste0(Path_Code, "04_DiffGene"))
source(paste0(Path_Code, "05_WGCNA.R"))
source(paste0(Path_Code, "06_KeyGeneExp.R"))
source(paste0(Path_Code, "07_GO.R"))







