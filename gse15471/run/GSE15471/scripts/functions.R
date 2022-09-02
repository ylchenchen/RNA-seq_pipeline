#获取探针矩阵
probeM2geneM <- function(ids,probeM){
  #01去除注释与矩阵不匹配的探针，保留共有探针，41297个探针共有
  ids=as.data.frame(ids)
  dat=probeM
  table(rownames(dat) %in% ids$probe_id)
  dat=dat[rownames(dat) %in% ids$probe_id,] #匹配ids与probeM共有探针
  ids=ids[match(rownames(dat),ids$probe_id),]#匹配ids与probeM共有探针
  dim(ids)
  dim(dat)
  
  #共有探针中含有冗余探针（冗余探针代表多个探针对应一个基因）
  #探针ID去冗余-取中位数表达量最大值（三步：新建一列，排序，去重），去冗余
  ids$median=apply(dat,1,median) #新建一列为每行的中位数
  ids=ids[order(ids$symbol,ids$median,decreasing = T),] #中位数排序
  ids=ids[!duplicated(ids$symbol),]#保留最大中位数对应的列（降序第一个不重）
  
  #去冗余后，获得基因表达矩阵
  geneM=dat[ids$probe_id,] #获得去冗余之后的dat/exp
  rownames(geneM)=ids$symbol #把ids的symbol列的每一行给dat作为行名
  geneM[1:4,1:4]  ##看是否成功将探针矩阵转换为基因矩阵
  return(geneM)
}



DEG_limma <- function(Group,probeM){
  #01输入数据准备(limma包分析获得6列数据，与差异统计相关)
  design=model.matrix(~Group) #limma三步得差异基因
  fit=lmFit(probeM,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf) 
  deg$probe_id=rownames(deg)
  deg <- inner_join(deg,ids,by="probe_id") #重获symbol列
  #为deg数据框加几列（上下调标记与富集分析的EntrezeID以及用于富集分析得EntrezeID）
  return(deg)
}

#：GO富集
GO_enrich <- function(deg){
  #(1)输入数据
  gene_up = deg$ENTREZID[deg$change == 'up'] 
  gene_down = deg$ENTREZID[deg$change == 'down'] 
  gene_diff = c(gene_up,gene_down)
  f=paste0(gse_number,"_GO.Rdata")
  #(2)运行--富集步骤耗时很长，所以用Rdata进行变量存储，方便下次运行
  if(file.exists(f)){
    load(f) #载入
  }else{
    ego <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,ont = "ALL",readable = TRUE)#GO的3个过程都画
    ego_BP <- enrichGO(gene = gene_diff,OrgDb= org.Hs.eg.db,ont = "BP",readable = TRUE) #注可只画BP
    #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    save(gene_up,gene_down,gene_diff,ego,ego_BP,file = f) #保存
  }
  return(ego)
}

# KEGG富集
KEGG_enrich <- function(deg){
  #（1）输入数据
  gene_up = deg[deg$change == 'up','ENTREZID'] 
  gene_down = deg[deg$change == 'down','ENTREZID'] 
  gene_diff = c(gene_up,gene_down)
  #（2）对上调/下调/所有差异基因进行富集分析
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  f=paste0(gse_number,"_KEGG.Rdata")
  if(!file.exists(f)){
    kk.up <- enrichKEGG(gene = gene_up,organism = 'hsa')
    kk.down <- enrichKEGG(gene =  gene_down,organism  = 'hsa')
    kk.diff <- enrichKEGG(gene = gene_diff,organism  = 'hsa')
    save(kk.diff,kk.down,kk.up,file = f)
  }
  load(f)
}

#按照q值挑选富集最显著的10条上下调KEGG通路
##小提示：选太多会导致图形很丑
top10updownKEGG <- function(kk.down,kk.up){
  down_kegg <- kk.down@result %>%
    filter(pvalue<0.05) %>% #筛选行
    mutate(group=-1) #新增列
  
  down_kegg=down_kegg[order(down_kegg$qvalue),] 
  if(nrow(down_kegg)>9){
    down_kegg=down_kegg[1:10,] 
  }else{
    down_kegg=down_kegg[1:nrow(down_kegg),] 
  }#取最显著的前十下调，少于十个全选
  
  up_kegg <- kk.up@result %>%
    filter(pvalue<0.05) %>%
    mutate(group=1)
  if(nrow(up_kegg)>9){
    up_kegg=up_kegg[1:10,] 
  }else{
    up_kegg=up_kegg[1:nrow(up_kegg),] 
  }#取最显著的前十上调，少于十个全选
  up.data=up_kegg
  down.data=down_kegg
  save(up.data,down.data,file="top10updown.Rdata")
}

# KEGG的GSEA富集
KeggGSEA_enrich <- function(deg){
  #1.输入数据准备
  geneList=deg$logFC
  names(geneList)=deg$ENTREZID
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  library(clusterProfiler)
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',#按需替换
                    #nPerm        = 1000,
                    minGSSize    = 10,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  tmp=kk_gse@result
  dim(tmp)
  kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')#按需替换
  #DOSE::setReadable():mapping geneID to gene Symbol
  tmp=kk@result
  dim(tmp)
  pro='comp1'
  write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
  save(kk,kk_gse,file = 'Rdata/gsea_kk.Rdata')
}


#挑选KEGG GSEA富集分析前6个最显著上调下调通路展示(用head与tail)
top6downupGSEA <- function(deg){
  down_k <- kk_gse[tail(order(kk_gse$enrichmentScore,decreasing = F)),];down_k$group=-1
  up_k <- kk_gse[head(order(kk_gse$enrichmentScore,decreasing = F)),];up_k$group=1
  dat=rbind(up_k,down_k)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  save(down_k,up_k,file="top6downup.Rdata")
  return(dat)
}





