
# path_GSM5276938 <- "E:/WZH/Single/GEO/Test_Seurat/Data/GSM5276938/"
# setwd(path_GSM5276938)
# BiocManager::install("Seurat")
# library("Seurat")
# library("dplyr")
# library("ggsci")  # 我想改变一下配色
# library("Matrix")
# library("cowplot")
# library("celldex")
# 
# library("SingleR")
# library("scater")
# library("SummarizedExperiment")
# packageVersion("Seurat")
########
# 注意 ################################################################################
# 参考集里面是logcounts矩阵，后面对于单细胞数据集也要做类似的处理。 导入UMI count矩阵 #
#######################################################################################
test.count=as.data.frame(pbmc[["RNA"]]@counts)

#####################################################
# 导入参考集，保证两个数据集的基因相同，然后log处理 #
###########################################################################################
load(file="./HPCA/hpca.se.RData")
common_hpca <- intersect(rownames(test.count), rownames(hpca.se))
hpca.se <- hpca.se[common_hpca,]
test.count_forhpca <- test.count[common_hpca,]
test.count_forhpca.se <- SummarizedExperiment(assays=list(counts=test.count_forhpca))
test.count_forhpca.se <- logNormCounts(test.count_forhpca.se)
###########################################################################################

######################################################
# 接下来是注释步骤，在这一步里，我只用了main大类注释 #############
# 还有一个fine小类注释，这里没有演示，因为我觉得小类注释不太准 #
################################################################
###main
pred.main.hpca <- SingleR(test = test.count_forhpca.se, ref = hpca.se, labels = hpca.se$label.main)
result_main_hpca <- as.data.frame(pred.main.hpca$labels)
result_main_hpca$CB <- rownames(pred.main.hpca)
colnames(result_main_hpca) <- c('HPCA_Main', 'CB')
###########################################
# 得到的结果是这样的，每个CB都有一个label #
# result_main_hpca <- read.table("HPCA_Main.txt")
# colnames(result_main_hpca) <- result_main_hpca[1, ]
# result_main_hpca <- result_main_hpca[-1, ]
head(result_main_hpca)
###########################################

####################################################
# 我们接下来要把这个结果整合到tpbmc@meta.data中 #
# 然后画tsne/umap展示一下 ##########################
##########################
pbmc@meta.data$CB=rownames(pbmc@meta.data)
pbmc@meta.data=merge(pbmc@meta.data,result_main_hpca,by="CB")
rownames(pbmc@meta.data)=pbmc@meta.data$CB



my_cols <- c('3'='#F68282',"midnightblue",'5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             "greenyellow","red","blue","brown")

DimPlot(pbmc, reduction = "umap", group.by = "HPCA_Main", pt.size=1,
        cols = as.character(my_cols))+theme(
          axis.line = element_blank(),
          axis.ticks = element_blank(),axis.text = element_blank()
        )


p5 <- DimPlot(pbmc, reduction = "umap", group.by = "HPCA_Main", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank())

p6 <- DimPlot(pbmc, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
fig_umap <- plot_grid(p6, p5, labels = c('ident','HPCA_Main'),rel_widths = c(2,3))

setwd(Path_Data)
ggsave(filename = paste0(SC_FileName[n_file], "_umap4.png")
       , plot = fig_umap, device = 'png', width = 36, height = 12, units = 'cm')




