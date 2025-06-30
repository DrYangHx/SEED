# 读取两个文件
pfi_data <- read.csv("D:/Desktop/医学/9-22/LUNG2/3-lassoPFIout0922.csv", header = TRUE)
os_data <- read.csv("D:/Desktop/医学/9-22/LUNG2/3-lassoOSout0922.csv", header = TRUE)

#pfi_data <- read.csv("D:/Desktop/医学/9-22/LUNG2/4-step_pfi.csv", header = TRUE)
#os_data <- read.csv("D:/Desktop/医学/9-22/LUNG2/4-step_os.csv", header = TRUE)

# 提取基因名称（假设基因名称在第3行及以后）
pfi_genes <- pfi_data$X[3:nrow(pfi_data)]
os_genes <- os_data$X[3:nrow(os_data)]

# 找出两个文件中的共同基因
common_genes <- intersect(pfi_genes, os_genes)

# 从两个文件中筛选出共同基因的行
pfi_filtered <- pfi_data[pfi_data$X %in% common_genes, ]
os_filtered <- os_data[os_data$X %in% common_genes, ]

# 将筛选后的数据重新加上前两行的生存信息
pfi_cleaned <- rbind(pfi_data[1:2, ], pfi_filtered)
os_cleaned <- rbind(os_data[1:2, ], os_filtered)

# 将结果保存为新的CSV文件
write.csv(pfi_cleaned, file = "D:/Desktop/医学/9-22/LUNG2/5-intersection_PFI_out.csv", row.names = FALSE)
write.csv(os_cleaned, file = "D:/Desktop/医学/9-22/LUNG2/5-intersection_OS_out.csv", row.names = FALSE)


