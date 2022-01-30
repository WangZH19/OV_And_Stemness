options(stringsAsFactors = FALSE)

Path_Code <- "E:/2022_GitHub/OV_and_Stemness/02_SingleCell/01_Code/"
Path_Data <- "E:/2022_GitHub/OV_and_Stemness/02_SingleCell/02_Data/"


setwd(Path_Data)

Color_N_tan <- unlist(read.table("N_tan.txt", stringsAsFactors = FALSE))
Color_P_greenyellow <- unlist(read.table("P_greenyellow.txt", stringsAsFactors = FALSE))

Color_gene <- list()
Color_gene[[1]] <- Color_P_greenyellow
Color_gene[[2]] <- Color_N_tan
names(Color_gene) <- c("P", "N")

Color_geneDist <- list()
for (n_1 in 1:length(Color_gene)) {########## for n_1
  Color_geneDist[[n_1]] <- list()
  n_int <- length(Color_gene[[n_1]]) %/% 9
  n_rem <- length(Color_gene[[n_1]]) %% 9
  genePosition_9 <- list()
  x <- 1
  y <- 9
  if (n_rem == 0) {number <- n_int} else {number <- n_int + 1}
  for (i in 1:number) {########## for i
    position <- seq(x, y, 1)
    if (n_rem != 0) {########## if n_rem
      if (i == (number - 1)) {######### if i
        x <- x + 9
        y <- y + 1
      }else {
        x <- x + 9
        y <- y + 9
      }
    }########## if i
    else {######### else n_rem
      x <- x + 9
      y <- y + 9
    }######### else n_rem
    print(position)
    Color_geneDist[[n_1]][[i]] <- Color_gene[[n_1]][position]
  }########## for i
}########## for n_1
names(Color_geneDist) <- c("Up", "Down")
####################################################################################################
####################################################################################################


##################################
# 设置读取单细胞数据的文件的位置 #
##################################
####################
# 导入需要的依赖包 #
#######################################
# BiocManager::install('Seurat')

# 单细胞数据文件夹
SC_FileName <- c("GSM5276940")
n_file <- 1

for (n_file in 1:length(SC_FileName)) {
  # for (n_file in 1:1) {
  print(paste0("runing file---", SC_FileName[n_file]))
  
  
  ############
  # 数据读取 #########################################
  # 参考博客：https://www.jianshu.com/p/03b94b2034d5 ############################
  ###############################################################################
  #单细胞读取数据，因为涉及到了稀疏矩阵，所以要与稀疏矩阵联动
  pbmc.counts <- Read10X(data.dir = paste0("./", SC_FileName[n_file], "/RawData/"))
  
  
  ###############################################################################
  
  ################
  # 设置存储路径 ####################################
  ###################################################
  #######################################
  # 创建存储基因在细胞亚群中的分布情况
  #######################################
  # setwd(path_picture)
  dir.create("./DistPicture/")
  dir.create("./DistPicture/N_picture")
  dir.create("./DistPicture/P_picture")
  #######################################
  #######################################
  
  
  ############
  # 创建对象 #
  #########################################################################################################
  # pbmc <- CreateSeuratObject(counts = pbmc.counts)
  pbmc <- CreateSeuratObject(counts = pbmc.counts, project = "pbmc3k", min.cells = 3, min.features = 200)
  
  
  ##########################
  # 线粒体细胞和红细胞比例 #
  #低质量或者垂死细胞会有线粒体污染，所以要进行线粒体基因占比计算 #
  ##################################################################################################
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  #人类血液常见红细胞基因
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
  #统计干细胞与数据中的RNA匹配数量
  HB_m <- match(HB.genes_total,rownames(pbmc@assays$RNA))
  HB.genes <- rownames(pbmc@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc,features=HB.genes)
  ##################################################################################################
  # head(pbmc@meta.data)[,c(2,3,4,5)]
  #########################################################################################################
  
  ##################
  # 均一化以标准化 #
  ###########################################################################################
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  ###########################################################################################
  
  #######################
  # 特征选择：高变基因  #
  ###################################################################################
  # 识别高度可变的特征(特征选择)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  ###################################################################################
  
  ########################
  # 找出10个最易变的基因 #
  ###########################################################################
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(pbmc), 10)
  ###########################################################################
  
  ###########
  # 标准化  #
  ###################################################
  # Scaling the data 
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  ###################################################
  
  ############
  # 数据降维 #
  #############################################################
  # RunPCA   
  # RunUMAP：进行群间分类较好，更加紧实
  # RunTSNE：进行群内分类较好，比较稀疏
  # 选取其中一个降维算法即可，因为他们的核心算法都是一样的，
  #############################################################
  
  #######
  # PCA #
  #######
  # 特征提取：PCA降维 #
  #Perform linear dimensional reductionPerform linear dimensional reduction
  #############################################################################
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  #############################################################################
  ############################
  # 查看一下 PCA 都有那些结果
  ############################
  # 每个细胞在PC轴上的坐标
  head(pbmc@reductions$pca@cell.embeddings)
  # 每个基因对每个PC轴的贡献度（loading值）
  head(pbmc@reductions$pca@feature.loadings)
  
  
  head(pbmc@reductions)
  head(pbmc@reductions$pca@key)
  
  
  ######################
  # 对Loading值一番研究
  ######################
  # Get the feature loadings for a given DimReduc
  t(Loadings(object = pbmc[["pca"]])[1:5,1:5])
  # Get the feature loadings for a specified DimReduc in a Seurat object
  t(Loadings(object = pbmc, reduction = "pca")[1:5,1:5])
  # Set the feature loadings for a given DimReduc
  new.loadings <- Loadings(object = pbmc[["pca"]])
  new.loadings <- new.loadings + 0.01
  Loadings(object = pbmc[["pca"]]) <- new.loadings
  
  ##########################################################
  # Examine and visualize PCA results a few different ways
  # 通过几种不同的方法检查和可视化PCA结果
  ##########################################################
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  
  
  
  # DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
  
  
  #########################################
  # 选择多少个维度进行下一阶段的分析呢？
  # 基于PAC有以下几种方法可以探索。
  #########################################
  # Determine the  "dimensionality" of the dataset
  # NOTE: This process can take a long time for big datasets, comment out for expediency.
  # More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
  # 可以看出大概在PC为某一个值的时候，每个轴是有区分意义的。
  pbmc <- JackStraw(pbmc, num.replicate = 100)
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
  
  # 每个轴显著相关的基因，这个也可以作为后面分析选择基因的一个参考。
  # ?PCASigGenes
  head(PCASigGenes(pbmc,pcs.use=2,pval.cut = 0.7))  
  
  ########
  # 聚类 #
  ########
  # 最后一步
  # Seurat采用的是graph-based聚类方法，k-means方法在V3中已经不存在了
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  
  
  table(pbmc@active.ident) # 查看每一类有多少个细胞
  
  # 提取某一类细胞。
  head(subset(as.data.frame(pbmc@active.ident),pbmc@active.ident=="2"))
  
  # 提取某一cluster细胞
  pbmc
  subpbmc<-subset(x = pbmc,idents="2")
  subpbmc
  
  # ?WhichCells
  head(WhichCells(pbmc,idents="2"))
  # Look at cluster IDs of the first 5 cells
  head(Idents(pbmc), 5)
  
  #提取部分细胞
  pbmc
  head(colnames(pbmc@assays$RNA@counts)[1:30])
  subbset<-subset(x=pbmc,cells=colnames(pbmc@assays$RNA@counts)[1:30])
  subbset
  
  #########################################################################################################
  #########################################################################################################
  ########### 当然，我们用的基本都是默认参数，建议？FindClusters一下
  ########### 看看具体的参数设置，比如虽然是图聚类，但是却有不同的算法，这个要看相应的文献了。
  #########################################################################################################
  # ?FindClusters
  # Algorithm for modularity optimization (1 = original Louvain algorithm;
  #                                        2 = Louvain algorithm with multilevel refinement;
  #                                        3 = SLM algorithm;
  #                                        4 = Leiden algorithm). Leiden requires the leidenalg python.
  #########################################################################################################
  #########################################################################################################
  
  ################
  # 系统发育分析 #
  #################################################
  # （Phylogenetic Analysis of Identity Classes）
  #################################################
  # ?BuildClusterTree
  # BiocManager::install('ape')
  # library(ape)
  pbmc<-BuildClusterTree(pbmc)
  Tool(object = pbmc, slot = 'BuildClusterTree')
  
  
  
  # ?CalculateBarcodeInflections
  pbmc<-CalculateBarcodeInflections(pbmc)
  SubsetByBarcodeInflections(pbmc)
  
  ##############
  # 可视化降维 #
  ####################################################################################################
  #############
  # UMAP 降维 #
  ####################################
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  
  head(pbmc@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
  
  #############
  # TSNE 降维 #
  ####################################
  pbmc <- RunTSNE(pbmc, dims = 1:10)
  head(pbmc@reductions$tsne@cell.embeddings)
  
  ############################
  # 比较一下两个可视化的结果 #
  # 两者的降维可视化的结构是一致的,
  # UMAP方法更加紧凑-在降维图上,同一cluster离得更近，不同cluster离得更远
  # 作为一种后来的算法有一定的优点，但是t-SNE结构也能很好地反映cluster的空间结构。
  ############################
  plot1<-DimPlot(pbmc, reduction = "umap",label = TRUE)+scale_color_npg()
  plot2<-DimPlot(pbmc, reduction = "tsne",label = TRUE)+scale_color_npg()
  CombinePlots(plots = list(plot1, plot2),legend="bottom")
  
  
  ######
  # 这里之前画出细胞分群的图
  
  ############
  # 差异分析 #
  ####################################################################################################
  # Finding differentially expressed features (cluster biomarkers) 
  # find all markers of cluster 1
  ##########################################
  # 这是一种one-others的差异分析方法，
  # 就是cluster1与其余的cluster来做比较
  # 当然这个是可以指定的,参数就是ident.2。
  ##########################################
  cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
  head(cluster1.markers, n = 5)
  # ?FindMarkers
  # find all markers distinguishing cluster 5 from clusters 0 and 3
  cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
  head(cluster5.markers, n = 5)
  
  

  
  cluster0_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(1, 2, 3, 4, 5, 6, 7, 8, 9), min.pct = 0.25)
  cluster1_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 2, 3, 4, 5, 6, 7, 8, 9), min.pct = 0.25)
  cluster2_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 3, 4, 5, 6, 7, 8, 9), min.pct = 0.25)
  cluster3_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 4, 5, 6, 7, 8, 9), min.pct = 0.25)
  cluster4_8.markers <- FindMarkers(pbmc, ident.1 = 4, ident.2 = c(0, 1, 2, 3, 5, 6, 7, 8, 9), min.pct = 0.25)
  cluster5_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 3, 4, 6, 7, 8, 9), min.pct = 0.25)
  cluster6_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 3, 4, 5, 7, 8, 9), min.pct = 0.25)
  cluster7_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 3, 4, 5, 6, 8, 9), min.pct = 0.25)
  cluster8_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 3, 4, 5, 6, 7, 9), min.pct = 0.25)
  cluster9_8.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = c(0, 1, 2, 3, 4, 5, 6, 7, 8), min.pct = 0.25)
  
  save(cluster0_8.markers, file = "cluster0_8.markers.RData")
  save(cluster1_8.markers, file = "cluster1_8.markers.RData")
  save(cluster2_8.markers, file = "cluster2_8.markers.RData")
  save(cluster3_8.markers, file = "cluster3_8.markers.RData")
  save(cluster4_8.markers, file = "cluster4_8.markers.RData")
  save(cluster5_8.markers, file = "cluster5_8.markers.RData")
  save(cluster6_8.markers, file = "cluster6_8.markers.RData")
  save(cluster7_8.markers, file = "cluster7_8.markers.RData")
  save(cluster8_8.markers, file = "cluster8_8.markers.RData")
  save(cluster9_8.markers, file = "cluster9_8.markers.RData")
  
  cluster357_8.markers <- FindMarkers(pbmc, ident.1 = c(3, 5, 7), ident.2 = c(0, 1, 2, 4, 6, 8, 9), min.pct = 0.25)
  save(cluster357_8.markers, file = "cluster357_8.markers.RData")
  
  
  logFCfilter <- 1  # 大于为上调，小于上调
  
  # cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name  
  # 如果是直接写出基因名的加这个参数，转化为基因ID的话不用加。
  
  
  
  # 看看输出结果都是什么。
  # ?FindMarkers
  
  # 我们还可以输出一个总表
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  # 这里可以注意一下only.pos 参数，可以指定返回positive markers 基因
  # test.use可以指定检验方法，可选择的有：wilcox，bimod，roc，t，negbinom，poisson，LR，MAST，DESeq2
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  head(pbmc.markers)
  
  ####################################################################################
  # pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)  # avg_logFC
  pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)  # avg_log2FC
  ####################################################################################
  # ?top_n
  
  ####################################
  # cluster间保守conserved marker基因.
  ####################################
  #Finds markers that are conserved between the groups
  #构建一个分组方式：
  pbmc[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = pbmc), replace = TRUE)
  head(FindConservedMarkers(pbmc, ident.1 = 0, ident.2 = 1, grouping.var = "groups"))
  
  
  # source("E:/WZH/01_Single/GEO/Test_SingleR/ttt.R")
  
  
  ####################################################################################################
  ####################################################################################################
  ########## 重要，在细胞群中看及关键基因的分布
  ####################################################################################################
  ####################################################################################################
  # setwd(paste0(path_AllFile, SC_FileName[n_file]))
  # 
  # for (n_1 in 1:length(Color_geneDist)) {########## for n_1
  #   if (n_1 == 1) {setwd(paste0(Path_Data, "/DistPicture/P_picture/"))}
  #   if (n_1 == 2) {setwd(paste0(Path_Data, "/DistPicture/N_picture/"))}
  #   
  #   for (n_2 in 1:length(Color_geneDist[[n_1]])) {########## for n_2
  #     
  #     file = paste0(names(Color_geneDist)[n_1], "_", n_2, ".png")
  #     png(file, width = 1920, height = 1080)
  #     # plot1 <-FeaturePlot(pbmc, features = c("USP18","SAMD12","EHF","UBALD2","BCL2L14","LRP8","SLC44A2","RTP4"),min.cutoff = 0, max.cutoff = 4)
  #     plot1 <- FeaturePlot(pbmc, features = c(as.character(na.omit(Color_geneDist[[n_1]][[n_2]]))),min.cutoff = 0, max.cutoff = 4)
  #     # plot1
  #     print(plot1)
  #     dev.off()
  #     
  #     plot1 <- FeaturePlot(pbmc, features = c(as.character(na.omit(Color_geneDist[[n_1]][[n_2]]))),min.cutoff = 0, max.cutoff = 4)
  #     top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  #     # 绘制热图
  #     plot2<-DoHeatmap(pbmc, features = top10$gene) + NoLegend()+scale_color_npg()
  #     # library(gridExtra)
  #     grid.arrange(plot1,plot2,ncol = 2, nrow = 1)
  #     # plot1<-FeaturePlot(pbmc, features = a,min.cutoff = 0, max.cutoff = 4)
  #     
  #   }########## for n_2
  #   # setwd(paste0(path_AllFile, SC_FileName[n_file]))
  # }########## for n_1
  # 
  # # setwd(paste0(path_AllFile, SC_FileName[n_file]))
  
}
getwd()
##############################################








##########################################################################################################################################
##########################################################################################################################################
# 基因的全分布图
##########################################################################################################################################
##########################################################################################################################################
KeyGene_Endometrioid <- c("CDH11", "CTHRC1", "VCAM1", "GFPT2", "COL6A1", "TMEM158", "TPM1", "AXL", "MMP2", "COL6A3", 
                          "TIMP2", "COL6A3", "MRC2")

KeyGene_HighGradeSerous <- c("FOXP1", "CDH11",  "CTHRC1",  "CLMP",  "GFPT2",  "TMEM158",  "SERPINF1",  "MMP2",  "FBN1", "COL6A2", 
                             "COL6A3",  "MRC2",  "BICC1", "SRPX2", "COL6A1", "DKK3", "LUM", "BGN", "PRRX1", "SPON2",
                             "PRKG1", "FAM19A5", "PDLIM3", "ADAMTS2", "TAGLN", "AEBP1")

# a <- c("CTHRC1", "")


# FeaturePlot(pbmc, features = KeyGene_Endometrioid,min.cutoff = 0, max.cutoff = 4)
# ?FeaturePlot


mat1 <- as.data.frame(pbmc@assays[["RNA"]]@data[KeyGene_HighGradeSerous[1], ])
colnames(mat1)="exp"
mat2=Embeddings(pbmc,"umap")
mat3=merge(mat2,mat1,by="row.names")

aa <- mat3
for (i in 2:length(KeyGene_HighGradeSerous)) {
  mat1 <- as.data.frame(pbmc@assays[["RNA"]]@data[KeyGene_HighGradeSerous[i], ])
  colnames(mat1)="exp"
  mat2=Embeddings(pbmc,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  
  aa$exp <- aa$exp + mat3$exp
}
aa%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp))+
  scale_color_gradient(low = "grey",high = "purple")+theme_bw()

DimPlot(pbmc, reduction = "umap",label = TRUE)+scale_color_npg()
# mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp))+
#   scale_color_gradient(low = "grey",high = "purple")+theme_bw()

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


################################################################################
################################################################################
# 获取目标群基因 
# UMAP_1 (0, 10)
# UMAP_2 (-9, 0)
################################################################################
################################################################################

################################################################################
# 通过使用限定坐标的方法获取指定类别中的细胞
# 与获取 pbmc@meta.data$seurat_clusters 中各个闭包分类对应的数值是一样的
################################################################################
# UMAP <- pbmc@reductions$umap@cell.embeddings %>% 
#   as.data.frame() %>% cbind(tx = pbmc@meta.data$seurat_clusters)
# 
# # 坐标获取
# T_cell.embeddings <- pbmc@reductions[["umap"]]@cell.embeddings
# 
# Cluster_4_t <- intersect(intersect(intersect(which(T_cell.embeddings[, "UMAP_1"] > 0),
#                                              which(T_cell.embeddings[, "UMAP_1"] < 5)),
#                                    which(T_cell.embeddings[, "UMAP_2"] > -9)),
#                          which(T_cell.embeddings[, "UMAP_2"] < 0))
################################################################################
################################################################################

# 获取各个细胞的类别

# 提取每个细胞类型对应的坐标系中的位置，与类别，数字直接就是代表的某一类
UMAP <- pbmc@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(tx = pbmc@meta.data$seurat_clusters)

pbmc_assays_data <- as.data.frame(pbmc@assays[["RNA"]]@data, strinsAsFactor = FALSE)

All_Cell <- colnames(pbmc_assays_data)
All_genes <- rownames(Cluster_4_assays)

# 获取类别 4 的细胞
Cluster_4_position <- which(UMAP$tx == 4)
Cluster_4_Cell <- rownames(UMAP)[Cluster_4_position]



Cell_4_positioin <- which((All_Cell %in% Cluster_4_Cell) == TRUE)

Cluster_4_assays <- pbmc_assays_data[, Cell_4_positioin]


Cluster_4_genes <- c()
for (n_row in 1:dim(Cluster_4_assays)[1]) {
  if (!(all(Cluster_4_assays[n_row, ] == 0))) {
    Cluster_4_genes <- append(Cluster_4_genes, All_genes[n_row])
    
  }
}


aaaaaaa <- 1
################################################################################
################################################################################



















