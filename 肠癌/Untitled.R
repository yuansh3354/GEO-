### ---
### Title: "Untitled"
### Create: "Yuansh"
### Date: "5/01/2020"
### Email: yuansh3354@163.com
### output: html_document
### ---


### step0 准备

##### 1. 设置镜像
if(T){
  options()$repos 
  options()$BioC_mirror
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  #setwd('/Volumes/Lexar/ZG/')
  #BiocManager::install('randomForestSRC')
  #install.packages('包')
}

##### 2. 导入包
if(T){
  library(limma)
  library(GEOquery)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(affy)
  library(gcrma)
  library(stringr)
  library(hgu133plus2.db )
  library(org.Hs.eg.db)
  library(GenomicFeatures)
  library(rtracklayer)
  library(biomaRt)
  library(glmnet)
  library(survival)
  library(Hmisc)
  library(clusterProfiler)
}

# 训练
file = 'GSE8671'
setwd('/Volumes/Lexar/ZG/肠癌')
if(T){
  gset <- getGEO(file, destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)
  a=gset[[1]]
  tab <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID"))
  e <- exprs(a)
  geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
    colMeans(e[ids,,drop=FALSE])
  }))
  geneExpr=geneExpr[apply(geneExpr,1,sd)>0,]
  geneExpr=log2(geneExpr+1)
  geneExpr=normalizeBetweenArrays(geneExpr)
  pd=pData(a) 
  # 1. 将样本信息导出
  write.csv(pd,file = paste(file,'_clinic.csv',sep = ''))
  write.csv(geneExpr,file = paste(file,'.csv',sep = ''))
  rm(a,gset)
}
paste(file,'_clinic.csv',sep = '')
pd = read.csv(paste(file,'_clinic.csv',sep = ''),header = T, row.names = 1)
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)

### step1 差异基因

##### 1.提取grouplist 并构建分组矩阵
if(T){
  group_list = as.character(pd[,45])
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(df_expr)
  print(paste(colnames(design)[1],colnames(design)[2],sep = '-'))
  contrast.matrix<-makeContrasts(paste(colnames(design)[1],colnames(design)[2],sep = '-'),
                                 levels = design)
}

##### 2. 寻找差异基因
deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}
deg = deg(df_expr,design,contrast.matrix)
deg = deg[order(deg$adj.P.Val),]
# 输出文件夹中的GSE13911-deg.csv就是差异基因
top = deg
top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),]
write.csv(top,paste(file,'-deg.csv',sep = ''))
##### 3. 提取满足条件的差异基因
top = deg
int_gene = rownames(top[which(top$adj.P.Val<0.05 & abs(top$logFC)>1),])[1:100]
id = int_gene
write.csv(top,'GSE13911-top-deg.csv')

##### 4.聚类
if(T){
  library(pheatmap)
  n = t(scale(t(df_expr[int_gene,])))
  n[n>2] = 2
  n[n< -2] = -2
  ac = data.frame(g=group_list)
  ac$names = colnames(n) 
  ac = ac[order(ac[,1]),]
  rownames(ac) = ac$names
  a = as.data.frame(ac[,1])
  colnames(a) = 'Type'
  rownames(a) = rownames(ac)
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           annotation_col=a,filename = '符合条件的基因聚类.png')
}

### 5.通路富集
##### 其中go是功能富集,kegg是通路富集
if(T){
  logFC_t=1
  top$g= 0
  top[which(top$logFC<= -1),]$g = 'DOWN'
  top[which(top$logFC>= 1),]$g = 'UP'
  
  top$ENTREZID=rownames(top)
  
  gene_up= top[top$g == 'UP','ENTREZID'] 
  gene_down=top[top$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  source('kegg_and_go_up_and_down.R')
  run_go(gene_up,gene_down,pro='NORMAL-TUMOR')
}
if(T){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_diff,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  barplot(kk.up )
  ggsave('kk.up.barplot.png')

  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_diff,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  barplot(kk.down )
  ggsave('kk.down.barplot.png')
  
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  barplot(kk.diff )
  ggsave('kk.diff.barplot.png')
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'kegg_up_down.png')
}#kegg
