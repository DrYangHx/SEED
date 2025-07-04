# 加载必要的库
library(survival)
library(survminer)
library(glmnet)
library(survivalROC)
library(caret)

# 读取数据
OSlasso <- read.csv("D:/Desktop/医学/9-22/LUNG2/2.2-enrichment_PFI_out.csv")
OSlasso[1:5,1:5]
row.names(OSlasso)<-OSlasso$X
OSlasso<-OSlasso[,-1]
OStran<-t(OSlasso)
OStran1<-as.data.frame(OStran)

##lasso回归
library(survival)
#install.packages("survminer")
library(survminer)
#install.packages("glmnet")
library(glmnet)
#install.packages("survivalROC")
library(survivalROC)
library(caret)
x<-as.matrix(OStran1[,3:ncol(OStran1)])
y<-data.matrix(survival::Surv(OStran1$time,OStran1$status))

# 生成一个从 -10000 到 10000 的随机数
random_number <- sample(-1000:1000, 1)

# 输出随机数
print(random_number)

set.seed(200)

##保存图1
cv.fit<-cv.glmnet(x,y,type.measure="deviance",family="cox",maxit=10^7,lambda.min.ratio=0.02)
plot(cv.fit)
##保存图2
cv.fit1<-glmnet(x,y,type.measure="deviance",family="cox")
cv.fit1
plot(cv.fit1,label=TRUE)

risk.score<-predict(cv.fit$glmnet.fit,newx=x,s=cv.fit$lambda.min,type="link")
OStran1$risk.score<-risk.score
dim(risk.score)
#setwd("E:\\datasetsrandom\\NKTCL\\lasso\\OS")
#write.csv(risk.score,file="risk.score.csv")


# 提取 LASSO 回归的非零系数
lasso.coef <- predict(cv.fit, s = cv.fit$lambda.min, type = "coefficients")
lasso.coef.matrix <- as.matrix(lasso.coef)


gene.names <- rownames(lasso.coef.matrix)[lasso.coef.matrix != 0]
gene.names <- gene.names[gene.names != "(Intercept)"]  # 去掉截距项


filtered.genes <- OSlasso[gene.names, ]  # 根据基因名筛选对应的行
filtered.genes <- rbind(OSlasso[1:2, ], filtered.genes)  # 保留前两行生存信息

# 保存筛选后的基因数据到文件
write.csv(filtered.genes, file = "D:/Desktop/医学/9-22/LUNG2/3-lassoPFIout0922.csv")

