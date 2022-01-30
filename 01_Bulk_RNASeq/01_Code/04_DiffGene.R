######
# 04 #
####################################################################################################
####################################################################################################
############################################## DiffGene ############################################
####################################################################################################
####################################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
# install.packages("pheatmap")
library("limma")  # 就是为了将重复的基因去均值

fdrFilter = 0.05                                                    #fdr临界值
logFCfilter = 1                                                     #logFC临界值
conNum = 88                                                         #normal组样品数目
treatNum = 416                                                      #tumor组样品数目

load("OV_tcga_gtex_tpm_list.RData")
rt1 <- data.frame(OV_tcga_gtex_tpm_list[[1]],
                  OV_tcga_gtex_tpm_list[[2]],
                  OV_tcga_gtex_tpm_list[[3]],
                  OV_tcga_gtex_tpm_list[[4]],
                  OV_tcga_gtex_tpm_list[[5]])

t_name <- gsub("\\.", "-",colnames(rt1))
colnames(rt1) <- t_name
rt <- rt1

#读取输入文件
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
# rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
# rownames(rt)=rt[,1]
# exp=rt[,2:ncol(rt)]
exp <- rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)  # 将重复的基因去均值
data=data[rowMeans(data)>0.2,]  # 去除低表达的数据

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]  # 提取基因的名字
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)  # 将表达数据与正常样本还是肿瘤样本合并
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)  # 我们做的是两组之间的差异，所以直接就写用这个检验做的差异分析
  conGeneMeans=mean(data[i,1:conNum])  #  在正常样本里面基因表达之的均值
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])  #  在肿瘤样本里面基因表达之的均值
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)  #log2肿瘤样本里面的均值 - log2正常样本里面的均值
  pvalue=wilcoxTest$p.value  # 得到P值
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }  # 输出数据：基因名称、正常样本里面的表达量、肿瘤样本里面的表达量、logFC、pValue
}
# 每一个基因都会经过一边上面的循环，每一个基因都可以得到一个p值
pValue=outTab[,"pValue"]  # 提取所以基因的p值
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")  # 对p值进行矫正，可以得到fdr值
outTab=cbind(outTab,fdr=fdr)  # 将fdr值加入到表格中，就可以得到差异基因的情况

setwd(path_10)
#输出所有基因的差异情况
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

#输出差异的表达量
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

#绘制火山图
pdf(file="vol.pdf",height=5,width=5)
xMax=15
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
diffSub_up <- diffSub
save(diffSub_up,file = "diffSub_up.RData")

points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
diffSub_down <- diffSub
save(diffSub_down,file = "diffSub_down.RData")


points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
geneNum=40
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
hmExp=data[hmGene,]
hmExp=log2(hmExp+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=6,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(40),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=8)
dev.off()

