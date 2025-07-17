##数据输入及预处理
dv<-read.csv("D:/004RNAseq/platform/Spearman/AA004.csv")#因变量输入路径
RNAseq<-read.csv("D:/004RNAseq/AA/RNAseq004.csv",header=T)#自变量输入路径
row.names(RNAseq)<-RNAseq$X
RNAseq<-RNAseq[,-1]
dv<-t(dv)
colnames(dv)<-colnames(RNAseq)
data<-rbind(dv,RNAseq)
datatran<-t(data)
datatran1<-as.data.frame(datatran)

dv<-as.numeric(as.character(datatran1$dv))
round(dv,0)
datatran1$dv<-dv

pFilter=1
fd<-datatran1
fd[,2:ncol(fd)]<-lapply(fd[,2:ncol(fd)],as.character)
fd[,2:ncol(fd)]<-lapply(fd[,2:ncol(fd)],as.numeric)

##spearman相关性分析
outTable=data.frame()
for(i in colnames(fd[,2:ncol(fd)])){
  cor<-cor.test(fd[,i],fd$dv,data=fd,method="spearman",alternative="two.sided",
  exact=FALSE)
  corP=cor$p.value
  corP[is.na(corP)]<-1
  if(corP<pFilter){
  outTable=rbind(outTable,
               cbind(id=i,
                     rho=cor$estimate,
                     pvalue=cor$p.value)
  )
 }
}

head(outTable)

##数据输出
setwd("D:\\004RNAseq\\platform\\Spearman")#输出路径
write.table(outTable,file="uniCox.txt",sep="\t",row.names=F,quote=F)