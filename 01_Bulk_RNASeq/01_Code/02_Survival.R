######
# 02 #
####################################################################################################
####################################################################################################
############################################# Survival #############################################
####################################################################################################
####################################################################################################
#install.packages("survival")
# install.packages("survminer")
# install.packages("gridtext")
# install.packages("ggtext")
library("survival")
library("survminer")
library("gridtext")
library("ggtext")

load("mRNAsiTime.RData")
rt <- mRNAsiTime
rt$futime=as.numeric(rt$futime)/365
var="mRNAsi"

rt$futate <- as.numeric(rt$futate)
rt$mRNAsi <- as.numeric(rt$mRNAsi)
mRNAsi_sort <- sort(rt[,var])
quantile <- c(10:90) / 100

a=ifelse(rt[,var]<=mRNAsi_sort[floor(413 * quantile[13])],"low","high")

diff=survdiff(Surv(futime, futate) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit=survfit(Surv(futime, futate) ~ a, data = rt)
if(pValue<0.001){
  pValue="<0.001"
}else{
  pValue=paste0("=",round(pValue,3))
}

surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=paste0("p",pValue),
                   pval.size=5,
                   risk.table=T,
                   legend.labs=c("high","low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table.height=.25)          
pdf(file=paste0(var,".survival.pdf"), width = 6.5, height = 5.5,onefile = FALSE)
print(surPlot)
dev.off()
####################################################################################################
####################################################################################################
