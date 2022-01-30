###################################
options(stringsAsFactors = FALSE)
###################################

load("All_Cell_Cluster.RData")
load("cluster4_8.markers.RData")
load("pbmc_all_assays_data.RData")

# 读取关键基因
Key_N <- read.table("N_tan.txt")[, 1]
Key_P <- read.table("P_greenyellow.txt")[, 1]
Key_NP <- c(Key_P, Key_N)

library("WGCNA")
#####################
# 获取内皮细胞位置
#####################
# Tissue_stem_cells
Tissue_stem_cells_position <- which(All_Cell_Cluster$Cluster_number == 4)

# 获取内皮细胞样本,并去除其他细胞样本
Tissue_stem_Sample <- All_Cell_Cluster$Cell_id[Tissue_stem_cells_position]
Tissue_stem_Sample <- Tissue_stem_Sample[-which(All_Cell_Cluster$Cell_Type[Tissue_stem_cells_position] != "Tissue_stem_cells")]

# 获取所有单细胞基因和样本
all_gene <- rownames(pbmc_all_assays_data)
all_sample <- colnames(pbmc_all_assays_data)

# 获取所有差异基因
diff4 <- rownames(cluster4_8.markers)
# 获取上调和下调基因
diff4_up <- diff8[which(cluster4_8.markers$avg_log2FC > 0.5)]
diff4_0 <- diff8[which(cluster4_8.markers$avg_log2FC == 0)]
diff4_down <- diff8[which(cluster4_8.markers$avg_log2FC < -0.5)]

# 获取内皮细胞以基差异基因的单细胞表达矩阵
diff4_gene_position <- which((all_gene %in% diff4) == TRUE)
Tissue_stem_Sample_position <- which(all_sample %in% Tissue_stem_Sample)

Tissue_stem_assays_data <- pbmc_all_assays_data[diff4_gene_position, Tissue_stem_Sample_position]

# 计算差异基因之间的相关性
diff4_gene_corAndPvalue <- corAndPvalue(t(Tissue_stem_assays_data))

diff4_gene_cor <- diff4_gene_corAndPvalue$cor
diff4_gene_p <- diff4_gene_corAndPvalue$p

diff4_gene_cor[which(is.na(diff4_gene_cor) == TRUE)] <- 0
diff4_gene_p[which(is.na(diff4_gene_p) == TRUE)] <- 1

# 获取关键基因的相关性矩阵
KeyNP_in_diff4 <- which((Key_NP %in% diff4) ==TRUE)
Key_assays_data <- diff4_gene_cor[, which((Key_NP %in% diff4) ==TRUE)]
Key_assays_data <- Key_assays_data[-KeyNP_in_diff4, ]  # 去除自身的基因


# 获取与关键基因相关的高相关基因----------定义为高风险基因
Higth_Risk_gene <- c()
for (i in 1:dim(Key_assays_data)[2]) {
  n_1 <- which(Key_assays_data[, i] > 0.5)
  Higth_Risk_gene <- append(Higth_Risk_gene, n_1)
  
}
Higth_Risk_gene <- unique(Higth_Risk_gene)
Higth_Risk_gene <- rownames(Key_assays_data)[Higth_Risk_gene]

Higth_Risk_gene <- Higth_Risk_gene[which((Higth_Risk_gene %in% diff8_up) == TRUE)]
Higth_Risk_gene <- unique(Higth_Risk_gene)




save(Higth_Risk_gene, file = "Higth_Risk_gene.RData")
write.table(Higth_Risk_gene, file = "Higth_Risk_gene.txt", row.names = FALSE, col.names = FALSE)


# 富集分析
setwd("./0902_GO/")
source(paste0(path_code, "0902_Smooth_muscle_HighGene_GO.R"))

Smooth_muscle_GO <- kk

write.csv(Smooth_muscle_GO, file = "Smooth_muscle_GO.csv")
setwd("..")

# 生存分析

source("1002_Smooth_muscle_cells_survival.R")








