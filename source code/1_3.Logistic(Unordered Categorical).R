##数据输入及预处理
dv<-read.csv("D:/004RNAseq/platform/Logistic/ORR004.csv")#因变量导入路径
RNAseq<-read.csv("D:/004RNAseq/ORR/RNAseq004.csv",header=T)#自变量导入路径
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

#logistic回归
pFilter=1
outTable=data.frame()
for(i in colnames(fd[,2:ncol(fd)])){
  logit<-glm(dv~fd[,i],data=fd,family=binomial(link='logit'))
  logitSummary=summary(logit)
  logitP=logitSummary$coefficients[-1,"Pr(>|z|)"]
  if (length(logitP)==0) {logitP=1}
  if(logitP<pFilter){
outTable=rbind(outTable,
               cbind(id=i,
                   estimate=logitSummary$coefficients[-1,"Estimate"],
                   Std=logitSummary$coefficients[-1,"Std. Error"],
                   zvalue=logitSummary$coefficients[-1,"z value"],
                   pvalue=logitSummary$coefficients[-1,"Pr(>|z|)"])
  )
 }
}

head(outTable)

##数据输出
setwd("D:\\004RNAseq\\platform\\Logistic")#输出路径
write.table(outTable,file="uniCox.txt",sep="\t",row.names=F,quote=F)