######
# 07 #
####################################################################################################
####################################################################################################
################################################# GO ###############################################
####################################################################################################
####################################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db")

NP_gene <- as.character(read.table("NP_gene.txt", stringsAsFactors = FALSE)[, 1])

entrezIDs <- mget(NP_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=data.frame(NP_gene,entrezID=entrezIDs)
out <- out[-which(out$entrezID == "NA"), ]
colnames(out)=c("symbol","entrezID")
write.table(out,file="NP_id.txt",sep="\t",quote=F,row.names=F)    #输出结果


path_GO <- "E:/WZH/01_2021_important/01_Project_new_20211215backups/20_GO/"
setwd(path_GO)                 #设置工作目录
# rt=read.table("N_id.txt",sep="\t",header=T,check.names=F)           #读取id.txt文件
rt=read.table("NP_id.txt",sep="\t",header=T,check.names=F)
load("diffSub_up.RData")
load("diffSub_down.RData")

NP_up <- rt[which((rt$symbol %in% diffSub_up$gene) == TRUE), ]
NP_down <- rt[which((rt$symbol %in% diffSub_down$gene) == TRUE), ]


rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

#GO富集分析
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =1,
               # qvalueCutoff = 0.05,
               ont="all",
               readable =T)

write.table(kk,file="NP_GO.xls",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.csv(kk,file="NP_GO.csv")


NP_up=NP_up[is.na(NP_up[,"entrezID"])==F,]                                 #去除基因id为NA的基因
gene=NP_up$entrezID
NP_GO_up <- enrichGO(gene = gene,
                     OrgDb = org.Hs.eg.db, 
                     pvalueCutoff =1,
                     # qvalueCutoff = 0.05,
                     ont="all",
                     readable =T)

NP_KEGG_up <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)   #富集分析


NP_down=NP_down[is.na(NP_down[,"entrezID"])==F,] 
gene=NP_down$entrezID
NP_GO_down <- enrichGO(gene = gene,
                       OrgDb = org.Hs.eg.db, 
                       pvalueCutoff =1,
                       # qvalueCutoff = 0.05,
                       ont="all",
                       readable =T)
NP_KEGG_down <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)   #富集分析


write.table(NP_GO_up,file="NP_GO_up.xls",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.csv(NP_GO_up,file="NP_GO_up.csv")
write.table(NP_GO_down,file="NP_GO_down.xls",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.csv(NP_GO_down,file="NP_GO_down.csv")

write.table(NP_KEGG_up,file="NP_KEGG_up.xls",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.csv(NP_KEGG_up,file="NP_KEGG_up.csv")
write.table(NP_KEGG_down,file="NP_KEGG_down.xls",sep="\t",quote=F,row.names = F)                 #保存富集结果
write.csv(NP_KEGG_down,file="NP_KEGG_down.csv")