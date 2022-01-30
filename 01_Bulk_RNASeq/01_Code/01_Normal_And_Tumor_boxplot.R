######
# 01 #
####################################################################################################
####################################################################################################
##################################### Normal And Tumor boxplot #####################################
####################################################################################################
####################################################################################################
# install.packages("beeswarm")
library(beeswarm)
yMin=0
yMax=1
ySeg=yMax*0.92

load("OV_mRNAsi.RData")
rt <- na.omit(OV_mRNAsi)
Type <- rt[, 2]

colnames(rt)=c("var","Type")

#正常组和肿瘤组差异比较
wilcoxTest<-wilcox.test(var ~ Type, data=rt)
wilcoxP=wilcoxTest$p.value
pvalue=signif(wilcoxP,4)
pval=0
if(pvalue<0.001){
  pval=signif(pvalue,4)
  pval=format(pval, scientific = TRUE)
}else{
  pval=round(pvalue,3)
}
outFile=paste0("mRNAsi",".pdf")
pdf(file=outFile,width=7,height=5)
par(mar = c(5,7,2,3))
boxplot(var ~ Type,data = rt,
        cex.main=1.5, cex.lab=1.3, cex.axis=1.2,ylim=c(yMin,yMax))
beeswarm(var ~ Type, data = rt, col = c("blue","red"),lwd=0.1,
         pch = 16, add = TRUE, corral="wrap")
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
text(1.5,ySeg*1.05,labels=paste("p=",pval,sep=""),cex=1.2)
dev.off()
####################################################################################################
####################################################################################################
