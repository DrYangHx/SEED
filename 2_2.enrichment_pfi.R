# 读取数据
#富集后的文件，final_GO
data <- read.csv("D:/Desktop/医学/9-22/LUNG2/enrichment/Enrichment_GO/_FINAL_GO.csv", header = TRUE)
#源数据，带筛选的源数据
os_cleaned <- read.csv("D:/Desktop/医学/9-22/LUNG2/0-train_PFI.csv", header = TRUE)  # 读取os_cleaned文件

# 查看列名以确保列名正确
colnames(data)

# 筛选出所有 Log(q-value) < -1.3 的行，假设列名是 Log.q.value.
filtered_hits <- subset(data, Log.q.value. < -1.3)

# 确认 Hits 列是否存在基因名
if (!"Hits" %in% colnames(filtered_hits)) {
  stop("Hits 列未找到，请确认列名正确。")
}

# 分割Hits列中的基因，Hits列使用"|"作为分隔符，修正分隔符为单个符号 "|"
genes_list <- unlist(strsplit(as.character(filtered_hits$Hits), split = "\\|"))

# 去除空值
genes_list <- genes_list[genes_list != ""]
genes_list
# 保留唯一的基因名称
unique_genes <- unique(genes_list)

# 检查基因名输出是否正确
print(unique_genes)

# 将结果转换为数据框，不用写入CSV文件，直接保留在内存中
result <- data.frame(Gene = unique_genes)

# 提取os_cleaned中的前两行（生存信息）
survival_info <- os_cleaned[1:2, ]

# 根据unique_genes从os_cleaned中提取对应的行，假设基因列名是"X"或其他列名
gene_rows <- os_cleaned[os_cleaned$X %in% result$Gene, ]  # 请根据实际基因列名修改"X"

# 合并前两行的生存信息与基因数据
final_data <- rbind(survival_info, gene_rows)

# 将结果写入新的CSV文件
write.csv(final_data, file = "D:/Desktop/医学/9-22/LUNG2/2.2-enrichment_PFI_out.csv", row.names = FALSE)
