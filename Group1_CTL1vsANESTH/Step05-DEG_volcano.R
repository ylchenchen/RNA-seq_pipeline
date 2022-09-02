rm(list = ls())
options(stringsAsFactors = F)

# 加载原始表达矩阵
load(file = "data/Step01-airwayData.Rdata")

# 读取3个软件的差异分析结果
load(file = "data//Step03-edgeR1_nrDEG.Rdata")
ls()

# 根据需要修改DEG的值
data <- DEG_edgeR1
colnames(data)


# 绘制火山图
library(ggplot2)
colnames(data)
p <- ggplot(data=data, aes(x=logFC, y=-log10(PValue),color=regulated)) + 
     geom_point(alpha=0.5, size=1.8) + 
  theme_set(theme_set(theme_bw(base_size=20))) + 
     xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c('blue','black','red'))
p
ggsave(p,filename = "result/CTLvsANESTH_EDGvolcano.pdf")






