# 清空环境变量
rm(list = ls())
# 分别调用该函数处理两个不同的数据文件，并将显著基因保存为向量

# 定义一个通用函数用于处理不同的数据文件并返回结果
cox_analysis <- function(input_file) {
  # 读取数据
  data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  
  # 将第一列（基因名和生存信息）设置为行名
  row.names(data) <- data$X
  data <- data[,-1]  # 删除第一列，因为它已作为行名
  
  # 转置数据集，使患者为行
  datatran <- t(data)
  datatran1 <- as.data.frame(datatran)
  
  # 将生存数据转换为数值型
  time <- as.numeric(as.character(datatran1$time))
  status <- as.numeric(as.character(datatran1$status))
  datatran1$time <- time
  datatran1$status <- status
  
  # 将剩下的列（基因表达数据）转换为数值型
  datatran1[,3:ncol(datatran1)] <- lapply(datatran1[,3:ncol(datatran1)], as.numeric)
  
  # 进行Cox回归分析
  library(survival)
  pFilter = 1  # p值阈值
  fd <- datatran1
  outTable <- data.frame()
  sigGenes <- c("time", "status")
  
  for(i in colnames(fd[,3:ncol(fd)])){
    cox <- tryCatch({
      coxph(Surv(time, status) ~ fd[,i], data = fd)
    }, error = function(e) {
      NULL  # 处理错误，跳过此基因
    })
    
    if (!is.null(cox)) {
      coxSummary <- summary(cox)
      coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
      coxP[is.na(coxP)] <- 1
      if(coxP < pFilter){
        sigGenes <- c(sigGenes, i)
        outTable <- rbind(outTable,
                          cbind(id = i,
                                HR = coxSummary$conf.int[,"exp(coef)"],
                                HR.95L = coxSummary$conf.int[,"lower .95"],
                                HR.95H = coxSummary$conf.int[,"upper .95"],
                                pvalue = coxSummary$coefficients[,"Pr(>|z|)"]))
      }
    }
  }
  
  # 筛选显著基因
  outTable$pvalue <- as.numeric(outTable$pvalue)
  significant_genes <- outTable[outTable$pvalue < 0.05, "id"]
  
  # 返回显著基因的列表
  return(significant_genes)
}

# 分别调用该函数处理两个不同的数据文件，并将显著基因保存为向量
os_genes <- cox_analysis("D:/Desktop/医学/9-22/LUNG2/0-train_OS.csv")
pfi_genes <- cox_analysis("D:/Desktop/医学/9-22/LUNG2/0-train_PFI.csv")


# 找出两个文件中的共同基因
common_genes <- intersect(os_genes, pfi_genes)

# 保存共同基因列表为CSV文件
write.csv(data.frame(Gene = common_genes), file = "D:/Desktop/医学/9-22/LUNG2/1-unicox_common_genes.csv", row.names = FALSE)


