# 加载所需库
library(survival)
library(timeROC)
library(ggplot2)
library(survminer)
library(caret)

# 设置文件路径（确保路径正确）
file_path <- "D:/Desktop/9-22/NKTCL/3-lassoOSout0922.csv"

# 数据导入和预处理
data <- read.csv(file_path)
row.names(data) <- data$X  # 保留第一列为行名
first_column <- data$X  # 保存第一列
data <- data[,-1]  # 去掉第一列
datatran <- t(data)  # 转置数据
datatran1 <- as.data.frame(datatran)  # 转换为数据框

# 将time和status数据转换为数值型
datatran1$time <- as.numeric(as.character(datatran1$time))
datatran1$status <- as.numeric(as.character(datatran1$status))
datatran1[, 3:ncol(datatran1)] <- lapply(datatran1[, 3:ncol(datatran1)], as.numeric)

# 动态生成多因素Cox回归公式
genes <- colnames(datatran1)[3:ncol(datatran1)]
cox_formula <- as.formula(paste("Surv(time, status) ~", paste(genes, collapse = "+")))

# 多因素Cox回归
mucox <- coxph(cox_formula, data = datatran1)
summary(mucox)

# 逐步回归
stepdata <- step(mucox, direction = "both")
summary(stepdata)

# 从逐步回归模型中提取显著基因并计算风险评分
selected_genes <- names(coef(stepdata))
coefficients <- coef(stepdata)


survival_info <- data[1:2, ]

# 合并生存信息和基因信息
final_output <- rbind(survival_info, gene_rows)

# 保存为 CSV 文件
write.csv(final_output, "D:/Desktop/9-22/NKTCL/4-step_os.csv", row.names = TRUE)

cat("生存信息和基因信息已保存至 '4-step_os.csv'.\n")
