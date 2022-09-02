rm(list = ls())
options(stringsAsFactors = F)

# 读取基因表达矩阵信息并查看分组信息和表达矩阵数据
lname <- load(file = "data/Step01-airwayData.Rdata")
lname

exprSet <- filter_count
dim(exprSet)
exprSet[1:4,1:4]
table(group_list)

# 加载包
library(edgeR)

# 假设数据符合正态分布，构建线性模型。0代表x线性模型的截距为0
design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(exprSet)
colnames(design) <- levels(factor(group_list))
design

# 构建edgeR的DGEList对象
DEG <- DGEList(counts=exprSet, 
               group=factor(group_list))
DEG$samples$lib.size

# 归一化基因表达分布
DEG <- calcNormFactors(DEG)
DEG$samples$norm.factors

# 计算线性模型的参数
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# 拟合线性模型
fit <- glmFit(DEG, design)

# 进行差异分析，1,-1意味着前比后
lrt <- glmLRT(fit, contrast=c(1,-1)) 

# 提取过滤差异分析结果
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

# 筛选上下调，设定阈值
fc_cutoff <- 2.0
pvalue <- 0.05

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),
                    which(DEG_edgeR$PValue<pvalue))
loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)


# 添加一列gene symbol
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENSEMBL", 
                  toType = "SYMBOL", 
                  OrgDb = org.Mm.eg.db)
head(id2symbol)
DEG_edgeR <- cbind(GeneID=rownames(DEG_edgeR),DEG_edgeR)
head(DEG_edgeR)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all=FALSE)
head(DEG_edgeR_symbol)

# 保存
write.table(DEG_edgeR_symbol,"result/DEG_edgeR_all.xls", row.names = F,
            sep="\t",quote = F)


## 取表达差异倍数和p值,矫正后的pvalue并保存
colnames(DEG_edgeR_symbol)
DEG_edgeR1 <- DEG_edgeR_symbol[,c(1,2,3,6,7,8)]
head(DEG_edgeR1)
save(DEG_edgeR1, file = "data/Step03-edgeR1_nrDEG.Rdata")


## 检查是否上下调设置错了
# 挑选一个差异表达基因
head(DEG_edgeR1)

exp <- c(t(express_cpm[match("ENSMUSG00000000001",rownames(express_cpm)),]))
exp
test <- data.frame(value=exp, group=group_list)

library(ggplot2)
ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()
ggsave('gene_exp_groups.pdf')



