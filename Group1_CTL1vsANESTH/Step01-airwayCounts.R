rm(list = ls())
options(stringsAsFactors = F)
library(stringr)

dir.create("data")
dir.create("result")

## ====================1.读取数据
# 读取raw count表达矩阵
# #Gse <- read.table("GSE182793_0620-rawcounts.txt",row.names = 1, 
#                        sep = "\t", header = T)

rawcount <- read.table("./raw_counts.txt",row.names = 1, 
                  sep = "\t", header = T)

colnames(rawcount) <- c("CTL1","CTL2","CTL3","CTL4",'ANESTH1','ANESTH2',
                        'ANESTH3','ANESTH4',"CONCIOUS1","CONCIOUS2","CONCIOUS3","CONCIOUS4")
a=head(rawcount)

# 查看表达谱
rawcount[1:4,1:4]

##选取需要比较的分组和提取表达矩阵子集
# 差异分析方案为：CTL vs ISCH_ANESTH
#有12个样，表达矩阵里的列名没有给具体的列名，顺序默认和GEO上的一样，
#前四个样（1：4）为CTL，5：8为ISCH_ANESTH（简写ANESTH），9：:12为ISCH_CONCIOUS（简写CONCIOUS）
rawcount=rawcount[,1:8]
## 去除前的基因表达矩阵情况
dim(rawcount)

# 获取分组信息
group=rep(c("CTL","ANESTH"),each=4)
group_list=factor(group,levels = c("CTL","ANESTH"))
table(group_list)

## =================== 2.表达矩阵预处理
# 过滤低表达基因
# 1.过滤在至少在75%的样本中都不表达的基因
# 2.过滤平均值count<10的基因
# 3.过滤平均cpm <10 的基因
keep <- rowSums(rawcount>0) >= floor(0.5*ncol(rawcount))
table(keep)

# FALSE  TRUE 
# 42725 12689 

filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

# 加载edgeR包计算counts per millio(cpm) 表达矩阵,并对结果取log2值
library(edgeR)
express_cpm <- log2(cpm(filter_count)+ 1)
express_cpm[1:6,1:6]

# 保存表达矩阵和分组结果
save(filter_count,express_cpm,group_list,
     file = "data/Step01-airwayData.Rdata")




