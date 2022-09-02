rm(list = ls())
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Mm.eg.db)

# 读取差异分析结果
load("data/Step01-airwayData.Rdata")

load(file = "data/Step03-edgeR1_nrDEG.Rdata")

# 提取所有差异表达的基因名
edgeR1_sigGene <- DEG_edgeR1[DEG_edgeR1$regulated!="normal",1]
head(edgeR1_sigGene)

# 根据需要更改DEG的值
DEG <- edgeR1_sigGene
head(DEG)
gene_all <- rownames(filter_count)


#### 第一步，从org.Hs.eg.db提取ENSG的ID 和GI号对应关系
keytypes(org.Mm.eg.db)

# bitr in clusterProfiler
allID <- bitr(gene_all, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Mm.eg.db )
head(allID)
degID <- bitr(DEG, fromType = "ENSEMBL", 
              toType = c( "ENTREZID" ), 
              OrgDb = org.Mm.eg.db )
head(degID)


R.utils::setOption( "clusterProfiler.download.method",'auto' )
# KEGG analysis----
# 设置pvalue与qvalue为最大值，输出所有结果，
# 然后根据结果来筛选显著性通路，
# 这样就不必因为没有显著性结果重新跑一边富集过程
enrich <- enrichKEGG(gene =degID[,2],
                     organism='mmu',
                     universe=allID[,2],
                     pvalueCutoff=1,
                     qvalueCutoff=1)
##注意物种，不然下面的代码会报错
# 计算富集因子
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))
head(GeneRatio)

BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])  ))
head(BgRatio)

enrich_factor <- GeneRatio/BgRatio

out <- data.frame(enrich$ID,
                  enrich$Description,
                  enrich$GeneRatio,
                  enrich$BgRatio,
                  round(enrich_factor,2),
                  enrich$pvalue,
                  enrich$qvalue,
                  enrich$geneID)

colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
write.table(out,"result/trut_VS_untrt_enrich_KEGG.xls",row.names = F,sep="\t",quote = F)

out_sig0.05 <- out[out$qvalue<0.01,]

# barplot
bar <- barplot(enrich,showCategory=20,title="KEGG Pathway",
               colorBy="p.adjust")
bar

# 保存
pdf(file = "result/kegg_bar_plot.pdf",width = 8,height = 6)
print(bar)
dev.off()

# dotplot
dot <- dotplot(enrich,x="geneRatio",showCategory=10,font.size=12,title="KEGG Pathway")
dot

# 保存
pdf(file = "result/kegg_dot_plot.pdf",width = 8,height = 6)
print(dot)
dev.off()


# GO 
enrich <- enrichGO(gene =degID[,2],OrgDb='org.Mm.eg.db',
                   ont="BP",universe=allID[,2],pvalueCutoff=1,qvalueCutoff=1)

# 计算富集因子
GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))

BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) 
  as.numeric(x[1])/as.numeric(x[2])))

enrich_factor <- GeneRatio/BgRatio

out <- data.frame(enrich$ID,
                  enrich$Description,
                  enrich$GeneRatio,
                  enrich$BgRatio,
                  round(enrich_factor,2),
                  enrich$pvalue,
                  enrich$qvalue,
                  enrich$geneID)

colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
write.table(out,"result/trut_VS_untrt_enrich_GO.xls",row.names = F,sep="\t",quote = F)

out_sig0.05 <- out[out$qvalue<0.01,]


# barplot
bar <- barplot(enrich,showCategory=10,title="Biological Pathway",colorBy="p.adjust")
bar

# 保存
pdf(file = "result/BP_bar_plot.pdf",width = 6,height = 6)
print(bar)
dev.off()

# dotplot
dot <- dotplot(enrich,x="geneRatio",showCategory=10,font.size=12,title="Biological Pathway")
dot

# 保存
pdf(file = "result/BP_dot_plot.pdf",width = 6,height = 6)
print(dot)
dev.off()





