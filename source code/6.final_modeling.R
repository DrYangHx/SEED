# 加载所需库
library(survival)
library(timeROC)
library(ggplot2)
library(survminer)
library(glmnet)
library(survivalROC)
library(caret)

# 加载训练数据
train_data <- read.csv("D:/Desktop/医学/9-22/成功/BRCA/5-intersection_OS_out.csv")
row.names(train_data) <- train_data$X
train_data <- train_data[,-1]

# 转置训练数据并进行数据处理
train_trans <- t(train_data)
train_df <- as.data.frame(train_trans)

# 将time和status数据转换为数值型
train_df$time <- as.numeric(as.character(train_df$time))
train_df$status <- as.numeric(as.character(train_df$status))
train_df[, 3:ncol(train_df)] <- lapply(train_df[, 3:ncol(train_df)], as.numeric)

# 动态生成多因素Cox回归公式
genes <- colnames(train_df)[3:ncol(train_df)]
cox_formula <- as.formula(paste("Surv(time, status) ~", paste(genes, collapse = "+")))

# 多因素Cox回归训练模型
cox_model <- coxph(cox_formula, data = train_df)
summary(cox_model)

# 提取模型系数用于测试数据
selected_genes <- names(coef(cox_model))
coefficients <- coef(cox_model)
coefficients
# 加载测试数据
validation_data <- read.csv("D:/Desktop/医学/9-22/成功/BRCA/0-train_OS.csv", header = TRUE, row.names = 1)

# 将status和time分别提取为数值型数据
status <- as.numeric(validation_data[1, ])
time <- as.numeric(validation_data[2, ])

# 提取基因表达数据并转置
gene_data <- t(validation_data[-c(1, 2), ])

# 将基因名与训练模型的基因名进行匹配
selected_gene_data <- gene_data[, selected_genes]

# 计算线性预测得分（基因表达值乘以相应的系数并求和）
scores <- as.matrix(selected_gene_data) %*% coefficients

# 将结果添加到数据框中
test_df <- data.frame(Sample = rownames(gene_data), Time = time, Status = status, multicoxscore = scores)

# 确认multicoxscore已经正确计算
test_df$multicoxscore

# 使用timeROC进行AUC分析
result <- with(test_df, timeROC(T = Time, delta = Status, marker = multicoxscore, cause = 1, times = c(12, 36, 60), iid = TRUE))

dat <- data.frame(fpr = as.numeric(result$FP),
                  tpr = as.numeric(result$TP),
                  time = rep(as.factor(c(12, 36, 60)), each = nrow(result$TP)))

# 使用 ggplot2 进行绘图
ggplot() + 
  geom_line(data = dat, aes(x = fpr, y = tpr, color = time), size = 1) + 
  scale_color_manual(name = NULL, values = c("#E69F00", "#56B4E9", "#009E73"),
                     labels = paste0("AUC of ", c(1, 2, 5), "-y survival: ",
                                     format(round(result$AUC, 2), nsmall = 2))) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "#2C3E50", size = 1) + 
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.background = element_blank(),
    legend.position = c(0.7, 0.2),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 确定最佳cutoff值
cutoff <- surv_cutpoint(test_df, time = "Time", event = "Status", variables = "multicoxscore")
cutoff
# 生成分类并绘制生存曲线
groups <- surv_categorize(cutoff)
fit <- survfit(Surv(Time, Status) ~ multicoxscore, data = groups)

ggsurvplot(fit, 
           data = groups, 
           pval = TRUE, 
           pval.method = TRUE, 
           palette = c("#E69F00", "#56B4E9"), 
           risk.table = TRUE, 
           conf.int = TRUE,
           ggtheme = theme_minimal() +
             theme(
               axis.text = element_text(size = 12, face = "bold", color = "black"),
               axis.title = element_text(size = 14, face = "bold", color = "black"),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
               legend.position = "right",
               legend.title = element_text(size = 12, face = "bold"),
               legend.text = element_text(size = 12, face = "bold"),
               panel.grid = element_blank(),
               panel.border = element_blank(),
               axis.line = element_line(size = 1, color = "black")
             )
)
