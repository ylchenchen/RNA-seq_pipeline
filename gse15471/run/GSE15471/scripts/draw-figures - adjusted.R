draw_p1_boxplot <- function(probeM){
  
  # 箱线图，使用ggplot需宽数据变长数据
  data <- as.data.frame(probeM)
  data <- melt(data)
  head(data)
  title <- paste (gse_number, "/", gpl_number, sep ="")
  p1 <- ggplot( data,aes(x=variable,y=value))+
    geom_boxplot()+
    theme_ggstatsplot()+
    theme(panel.grid = element_blank(),
          axis.text=element_text(size=10,face = 'bold'),
          axis.text.x=element_text(angle=90),
          plot.title = element_text(hjust = 0.5,size =15))+
    xlab('')+
    ylab('')+
    ggtitle(title)
  p1    ##利用箱线图初看表达矩阵整体表达情况
  ggsave(p1,filename = "png/p1.png")  ##想存你可以存在当前目录下
 
}

draw_p2_PCA <- function(probeM,gse_number,Group){
  #PCA
  exp=t(probeM)#画PCA图时要求行名是样本名，列名是探针名，因此需要转换
  exp=as.data.frame(exp)#将matrix转换为data.frame 
  dat.pca <- PCA(exp, graph =FALSE)
  # 画图，PCA图p2
  this_title <- paste0(gse_number,'_PCA')
  p2 <- fviz_pca_ind(dat.pca,
                     geom.ind = "point", # show points only (not "text")
                     col.ind = Group, # color by groups
                     palette = c("#00AFBB", "#FC4E07"),
                     addEllipses = TRUE, # Concentration ellipses
                     legend.title = "Groups")+
    ggtitle(this_title)+
    theme_ggstatsplot()+
    theme(plot.title = element_text(size=15,hjust = 0.5))
  
  p2
  ggsave(p2,filename = "png/p2.png") ##ggplot包绘图的图片保存方式ggsave
  return(p2)
}


draw_p3_pheatmap <- function(geneM,Group,gse_number){
  cg=names(tail(sort(apply(geneM,1,sd)),1000))#sd排序前1000个
  n=t(scale(t(geneM[cg,]))) 
  n[n> 2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  ac=data.frame(Group)
  rownames(ac)=colnames(n)
  # 画图，高变基因的表达量热图p3
  png("png/p3.png") ##保存基础包绘制图三部曲；png();画图；dev.off()保存
  p3 <- pheatmap::pheatmap(n,
                           show_colnames =F,
                           show_rownames = F,
                           main = gse_number,
                           annotation_col=ac,
                           breaks = seq(-3,3,length.out = 100),color = colorRampPalette(colors = c("blue","white","red"))(100)))
  #因为已经手动设置了表达量最大值，所以，可以不用设置break
  p3
  dev.off()
  return(p3)
}

draw_p4_colsample <- function(Group,exp,gse_number) {
  colD=data.frame(Group)
  exp=t(probeM)#画PCA图时要求行名是样本名，列名是探针名，因此需要转换
  exp=as.data.frame(exp)#将matrix转换为data.frame 
  exprSet=t(exp)#转置样品名为行名
  rownames(colD)=colnames(exprSet)
  png("png/p4.png")
  p4 <- pheatmap::pheatmap(cor(exprSet),
                           annotation_col = colD,
                           show_rownames = F,
                           show_colnames = F,
                           main = gse_number,
                           color = colorRampPalette(colors = c("blue","white","red"))(100)
  )
  p4
  dev.off()
  return(p4)
}



draw_p5_volcano <- function(Group,deg,gse_number) {
  #火山图绘制
  p5 <- ggplot(data = deg, 
               aes(x = logFC, 
                   y = -log10(P.Value))) +
    geom_point(alpha=0.4, size=3.5, 
               aes(color=change)) +
    ylab("-log10(Pvalue)")+
    scale_color_manual(values=c("blue", "grey","red"))+
    geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) + 
    geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
    theme_bw()##有geom_代表是图层
  p5
  ggsave("png/p5.png")
  return(p5)
}



draw_p6_geneheatmap <- function(Group,n) {
  png("png/p6.png")
  annotation_col=data.frame(group=Group)
  rownames(annotation_col)=colnames(n) 
  #rownames(annotation_col)=colnames(n) 
  p6 <- heatmap_plot <- pheatmap(n,show_colnames =F,
                                 scale = "row",
                                 #cluster_cols = F, 
                                 annotation_col=annotation_col,
                                 breaks = seq(-3,3,length.out = 100),
                                 color = colorRampPalette(colors = c("blue","white","red"))(100)
  ) 
  p6
  dev.off()
  return(p6)
}

# GO图
draw_p7_GObarplot<- function(ego){
  p7<-barplot(ego, split = "ONTOLOGY", font.size = 10, 
              showCategory = 5) + 
    facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") 
  ggsave("png/p7.png")
}

# KEGG图
draw_p8_KEGGBarplot <- function(up.data,down.data){
  dat=rbind(up.data,down.data)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  gk_plot <- ggplot(dat,aes(reorder(Description, pvalue), y=pvalue)) +
    geom_bar(aes(fill=factor((pvalue>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
    coord_flip() +
    scale_fill_manual(values=c("#0072B2", "#B20072"), guide="none") +
    labs(x="", y="" ) +
    theme_pander()  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.line.x = element_line(size = 0.3, colour = "black"),#x轴连线
          axis.ticks.length.x = unit(-0.20, "cm"),#修改x轴刻度的高度，负号表示向上
          axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),##线与数字不要重叠
          axis.ticks.x = element_line(colour = "black",size = 0.3) ,#修改x轴刻度的线                         
          axis.ticks.y = element_blank(),
          axis.text.y  = element_text(hjust=0),
          panel.background = element_rect(fill=NULL, colour = 'white')
    )
  ggsave("png/p8.png")
}


# GSEA6条上下通路
draw_p9_gsea <- function(dat){
  p9 <- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="#34bfb5",high="#ff6633",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + 
    theme_ggstatsplot()+
    theme(plot.title = element_text(size = 15,hjust = 0.5),  
          axis.text = element_text(size = 12,face = 'bold'),
          panel.grid = element_blank())+
    ggtitle("Pathway Enrichment") 
  p9
  ggsave("png/p9.png")
  return(p9)
}

#
draw_p10_GSEAdown <- function(kk){
  p10 <- gseaplot2(kk, geneSetID = rownames(down_k))
  p10
  ggsave("png/p10.png")
}

draw_p11_GSEAup <- function(kk){
  p11=gseaplot2(kk, geneSetID = rownames(up_k))
  p11
  ggsave("p11.png")
}



