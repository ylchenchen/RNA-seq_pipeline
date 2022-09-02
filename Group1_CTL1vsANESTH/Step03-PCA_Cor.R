# 魔幻操作，一键清空
rm(list = ls())  
options(stringsAsFactors = F)

# 加载数据并检查
lname <- load(file = 'data/Step01-airwayData.Rdata')
lname

dat <- express_cpm
dat[1:4,1:4]
dim(dat)


## 1.样本之间的相关性-层次聚类树----
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)

pdf(file = "result/sample_Treeplot.pdf",width = 6,height = 8)
plot(sampleTree)
dev.off()


## 2.样本之间的相关性-PCA----
# 第一步，数据预处理
dat <- as.data.frame(t(dat))
dat$group_list <- group_list


# 第二步，绘制PCA图
library(FactoMineR)
library(factoextra)

# 画图仅需要数值型数据，去掉最后一列的分组信息
dat_pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
class(dat_pca)

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "point", # 只显示点，不显示文字
                  col.ind = dat$group_list, # 用不同颜色表示分组
                  palette = c("#00AFBB", "#E7B800"),
                  addEllipses = T, # 是否圈起来
                  legend.title = "Groups") + theme_bw()
p
ggsave(p,file = "result/sample_PCA.pdf")


## 3.样本之间的相关性-cor----
exprSet <- express_cpm

library(corrplot)
dim(exprSet)

# 计算相关性
M <- cor(exprSet)
M
g <- corrplot(M,order = "AOE",addCoef.col = "white")

corrplot(M,order = "AOE",type="upper",tl.pos = "d",method = "circle")
corrplot(M,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,
         tl.pos="n", cl.pos="n")

# 绘制样本相关性的热图
library(pheatmap)
anno <- data.frame(sampleType=group_list)
rownames(anno) <- colnames(exprSet)
anno

p <- pheatmap::pheatmap(M,display_numbers = T,
                        annotation_col = anno,
                        fontsize = 12,cellheight = 30,
                        cellwidth = 30,cluster_rows = T,
                        cluster_cols = T)
p

pdf(file = "result/sample_cor.pdf")
print(p)
dev.off()




