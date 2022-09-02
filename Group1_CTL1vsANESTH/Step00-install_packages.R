rm(list = ls())

## Installing R packages
bioPackages <-c( 
  "corrplot","ggrepel",           # 绘制相关性图形
  "stringr",                      # 处理字符串的包
  "FactoMineR","factoextra",      # PCA分析软件
  "limma","edgeR","DESeq2",       # 差异分析的三个软件包
  "clusterProfiler", "org.Hs.eg.db",  # 安装进行GO和Kegg分析的扩展包
  "GSEABase","GSVA"               # 安装进行GSEA分析的扩展包
  )

## If you are in China, run the command below
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")) 
options(download.file.method = 'libcurl')
options(url.method='libcurl')

# 检查是否设定完毕
getOption("BioC_mirror")
getOption("CRAN")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 安装devtools管理github上的软件包
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")


## Installing missing packages
lapply( bioPackages, 
        function( bioPackage ){
          if(!bioPackage %in% rownames(installed.packages())){
              CRANpackages <- available.packages()
              if(bioPackage %in% rownames(CRANpackages)){
                 install.packages( bioPackage)
              }else{
                 BiocManager::install(bioPackage,suppressUpdates=F,ask=F)
              }
          }
        })


## 验证R包是否安装成功
library(limma)
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)

