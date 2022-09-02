rm(list = ls())
options(stringsAsFactors = F)

# 加载原始表达矩阵
load(file = "data/Step01-airwayData.Rdata")

# 读取3个软件的差异分析结果
load(file = "data/Step03-edgeR1_nrDEG.Rdata")
ls()

# 提取所有差异表达的基因名
edgeR1_sigGene <- DEG_edgeR1[DEG_edgeR1$regulated!="normal",1]
head(edgeR1_sigGene)

# 绘制热图
dat <- express_cpm[match(edgeR1_sigGene,rownames(express_cpm)),]
dat[1:4,1:4]
group <- data.frame(group=group_list)
rownames(group)=colnames(dat)
group

# 加载包
library(pheatmap)
p <- pheatmap(dat,scale = "row",show_colnames =T,show_rownames = F, 
              cluster_cols = T, 
              annotation_col=group,
              main = "edgeR's DEG")
p
ggsave(p,filename = "result/CTLvsANESTH_EDGpheatmap.pdf")









