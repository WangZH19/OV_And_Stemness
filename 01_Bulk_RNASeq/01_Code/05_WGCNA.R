######
# 05 #
####################################################################################################
####################################################################################################
################################################ WGCNA #############################################
####################################################################################################
####################################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
#install.packages("WGCNA")
library("WGCNA")

load("Sample_HighLow.RData")
data=read.table("diffGeneExp.txt",sep="\t",header=T,check.names=F,row.names=1)
load("OV_mRNAsi.RData")
OV_mRNAsi <- na.omit(OV_mRNAsi)

All_Sample <- substring(colnames(data), 1, 12)

data=data[, which((All_Sample %in% Sample_HighLow[[2]]) == TRUE)]
datExpr0=t(data)

################
###检查缺失值
################
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr0_sample <- substring(rownames(datExpr0), 1, 12)
rownames(datExpr0) <- datExpr0_sample
#################################################################################################################

##############
# 样品聚类
# 删除离散的样品
##############
sampleTree = hclust(dist(datExpr0), method = "average")

clust = cutreeStatic(sampleTree, 37000, minSize = 5)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]


traitData = OV_mRNAsi
fpkmSamples = rownames(datExpr0) # 表达的文件
traitSamples =rownames(traitData) # 临床的文件

traitSamples <- substring(traitSamples, 1, 12)

sameSample=intersect(fpkmSamples,traitSamples)
datExpr0=datExpr0[sameSample,] # 保留交集的样品，为了后面向相关性的分析
datTraits=traitData[sameSample,]
datTraits <- datTraits[, -2]
##############################################################################


sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
######################################################################

enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
###########################################################################################################

##################
###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
# softPower <- 5
adjacency = adjacency(datExpr0, power = softPower)

###TOM矩阵【计算两个基因之间的距离】
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM  # 相当于邻接矩阵
###########################################################################################################


###############
###基因聚类
###############
geneTree = hclust(as.dist(dissTOM), method = "average");
##########################################################################################

#########################
###动态剪切模块识别
##################################################################################################
# 模块基因数目【就是每个模块最少要有50个基因，如果模块的基因少于50，就只能与其他的模块进行合并】
# 这个是可以自己进行修改，如果得到的模块特别多数值放大点，太少的话就放小点
minModuleSize = 80

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
########################################################################

#########################################################################
# 相似模块聚类：对模块进行剪切，合并相似的模块
# 剪切的高度：需要我们结合画出来的图进行确定，进行合并
# 合并低于剪切先一下的模块
#########################################################################
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")

######################################################################################


######################################## 
# 相似模块合并：为了得到更少的模块
########################################
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

####################################################################
# 就是模块和临床的相关性：这里是模块和干细胞指数的额相关性
# 模块与性状数据热图
# 后续分析我们可以只拿出来两个模块进行分析
# 最负相关的模块和最正相关的模块
####################################################################
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# pdf(file="8_Module_trait.pdf",width=6,height=6)
pdf(file="8_Module_trait.pdf",width=10,height=10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
# par(mar = c(6, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               # xLabels = names(datTraits),
               xLabels = "mRNAsi",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

######################################
###计算MM和GS值
# MM：基因和模块的相关性；
# GS：基于和干细胞指数的相关性
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###输出GS_MM数据
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)

GS_MM <- geneInfo
save(GS_MM, file = "GS_MM.RData")
####################################################################################################
####################################################################################################
