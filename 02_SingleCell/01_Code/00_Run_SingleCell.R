options(stringsAsFactors = FALSE)

Path_Code <- "./02_SingleCell/01_Code/"
Path_Data <- "./02_SingleCell/02_SingleCell_Data/"

setwd(Path_Data)

library("Seurat")
library("dplyr")
library("ggsci")
library("Matrix")
library("cowplot")
library("gridExtra")
library("ape")
library("celldex")
library("SingleR")
library("scater")
library("SummarizedExperiment")
library("psych")
library("qgraph")
library("igraph")
library("tidyverse")

source(paste0(Path_Code, "01_Seurat.R"))
source(paste0(Path_Code, "02_SingleR.R"))
source(paste0(Path_Code, "03_CellPhonedb.R"))
source(paste0(Path_Code, "04_Tissue_stem_cells"))
source(paste0(Path_Code, "05_Endothelial_cells.R"))
source(paste0(Path_Code, "06_Smooth_muscle_cells.R"))





