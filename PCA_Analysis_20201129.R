# Load library
library(tidyr)
library(dplyr)
library(ggplot2)
#Loading data
data1 <- read.table("C:/Users/Yangyang/Desktop/xyq/GSE106688_genes.fpkm_table.txt", 
                    header = TRUE, quote = '\t')
head(data1)
# Change colomn names
Col_Names <- c("hESC R1","hESC R2","MES R1","MES R2","CP R1","CP R2", 
               "CM R1", "CM R2", "Fetal R1", "Fetal R2")
names(data1)[2:11] <- Col_Names
head(data1)
#去掉全为零的行
data2 <- data1[, -1]
head(data2)
data3 = data2[which(rowSums(data2==0)==0),]
head(data3)
#翻转数据
data_t <- t(data3)
head(data_t)
#Read groups
group <- read.delim("C:/Users/Yangyang/Desktop/xyq/group.txt", header = TRUE)
group
#主成分计算
pca_data <- prcomp(data_t, scale = T)
#查看合适主成分个数
screeplot(pca_data, type = "lines")
summary(pca_data)
#查看行名，确认是否为36个样品的名称
rownames(pca_data$x)
#提取PC1的百分比
x_per <- round(summary(pca_data)$importance[2, 1]*100, 0)
#提取PC2的百分比
y_per <- floor(summary(pca_data)$importance[2, 2]*100)
#按照样品名称添加组并且合并
df_sample <- data.frame(samplenames=rownames(pca_data$x), pca_data$x) %>%
  left_join(group, by = "samplenames") 
#数据尾部添加变量Group，只截取部分数据展示
df_sample
#ggplot绘图
#设置适合科研的背景色
theme_set(theme_bw())
#绘图
plot_pca <- ggplot(df_sample, mapping = aes(x = PC1, y = PC2, label = samplenames, color = Group)) +
  geom_point(aes(color=Group),size=5) + 
  
  geom_text(size = 3, ,hjust = 0.5, vjust = -1)+
  xlab(paste("PC1","(", x_per,"%)",sep=" ")) +
  ylab(paste("PC2","(", y_per,"%)",sep=" ")) + 
  scale_color_discrete(limits = c(paste("G", seq(1:5),sep = ""))+ 
                         theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

print(plot_pca)