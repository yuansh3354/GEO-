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
  library(clusterProfiler)
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
}

##### 3. 数据下载
file = 'GSE66229'
setwd('/Volumes/Lexar/ZG/胃癌')
df_expr = read.csv(paste(file,'.csv',sep = ''),header = T, row.names = 1)
pd = read.csv('survival.csv',header = T)
pd = na.omit(pd)
geoid = pd[,1]
sit = which(colnames(df_expr) %in% geoid)
GEO_ID = c(colnames(df_expr)[sit],colnames(df_expr)[-sit])
label = c(rep('tumor',300),rep('normal',100))
clinic = cbind(GEO_ID,label)
pd = merge(clinic,pd,all = T)
rownames(pd) = pd$GEO_ID
pd$GEO_ID == colnames(df_expr)
sit = which(pd$class!=3)
pd = pd[-sit,]
df_expr = df_expr[,-sit]
pd$Death = ifelse(pd$Death == 1 ,0,1)
if(T){
  group_list = as.character(pd[,2])
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
top = top[which(top$adj.P.Val<0.05 & abs(top$logFC)>0.3),]
write.csv(top,paste(file,'-deg.csv',sep = ''))

##### 3. 提取满足条件的差异基因
if(T){
  top = deg
  top = top[which(top$adj.P.Val<0.05),]
  top = top[order(top$adj.P.Val),]
  top = top[1:100,]
  int_gene = rownames(top[which(top$adj.P.Val<0.05& abs(top$logFC)>0.5),])
  id = rownames(int_gene)
  write.csv(top,paste(file,'-top-deg.csv',sep = ''))
}
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
  top = deg
  top = top[which(top$adj.P.Val<0.05& abs(top$logFC)>0.3),]
  logFC_t=1
  top$g= 0
  top[which(top$logFC<= 0),]$g = 'DOWN'
  top[which(top$logFC>= 0),]$g = 'UP'
  
  top$ENTREZID=rownames(top)
  
  gene_up= top[top$g == 'UP','ENTREZID'] 
  gene_down=top[top$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  source('kegg_and_go_up_and_down.R')
  #run_go(gene_up,gene_down,pro='NORMAL-TUMOR')
}
if(T){
if(T){
    kk.up <- enrichKEGG(gene         = gene_up,
                        organism     = 'hsa',
                        universe     = gene_diff,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
    head(kk.up)[,1:6]
    barplot(kk.up )
    ggsave('kk.up.barplot.png')
    write.csv(kk.up,paste(file,'-kk_up.csv',sep = ''),row.names = F)
    
    kk.down <- enrichKEGG(gene         =  gene_down,
                          organism     = 'hsa',
                          universe     = gene_diff,
                          pvalueCutoff = 0.9,
                          qvalueCutoff =0.9)
    head(kk.down)[,1:6]
    barplot(kk.down )
    ggsave('kk.down.barplot.png')
    write.csv(kk.down,paste(file,'-kk_down.csv',sep = ''),row.names = F)
    
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
}
### step2 cox分析
##### 1.构建数据:其中sample是所需表达谱,survival是所需生存信息
if(T){
  top = deg
  top = top[which(top$adj.P.Val<0.05& abs(top$logFC)>0.3),]
  top = top[order(top$adj.P.Val),]
  IDs = rownames(top)
  gene = paste0('ID_',IDs)
  sample = df_expr[IDs,]
  survival = pd[,c(5,4)] 
  survival = na.omit(survival)
  sample = sample[,colnames(sample) %in% rownames(survival)]
  gene = gsub(gene,pattern = '-', replacement = '_')
  rownames(sample) = gene
  colnames(survival) = c('time', 'status')
  cox_dat = as.data.frame(cbind(survival,t(sample))) 
  cox_dat[,1] = as.numeric(cox_dat[,1])
  cox_dat[,2] = as.numeric(cox_dat[,2])
  
} 
##### 2.构建cox模型 识别和预后相关的基因
if(T){
  library("survival")
  library("survminer")
  library(clusterProfiler)
  library(stringr)
  cox_analy = function(gene,survival_info){
    uni_cox = function(single_gene){
      formula = as.formula(paste0('Surv(time,status)~',single_gene))
      surv_uni_cox = summary(coxph(formula, data = cox_dat)) 
      ph_hypothesis_p = cox.zph(coxph(formula,data = cox_dat))$table[1:3]
      if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p > 0.05){
        single_cox_report = data.frame(
          'ID'=single_gene,
          'beta' = surv_uni_cox$coefficients[,1],
          'Hazard_Ratio' = exp(surv_uni_cox$coefficients[,1]),
          'z_pvalue'=surv_uni_cox$coefficients[,5],
          'Wald_pvalue'= as.numeric(surv_uni_cox$waldtest[3]),
          'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
        single_cox_report
      }
    }
    uni_cox_list = lapply(gene,uni_cox)
    do.call(rbind,uni_cox_list)
  }
  a = gene
  uni_cox_df = cox_analy(a,cox_dat)
  cox_IDs = str_split(uni_cox_df[,1],'_',simplify = T) 
  uni_cox_df[,1] = cox_IDs[,2]
  cox_IDs = cox_IDs[,2]
}
uni_cox_df$z_pvalue_adjust = p.adjust(uni_cox_df$z_pvalue ,method = "BH")
write.csv(uni_cox_df,paste(file,'-cox-gene.csv',sep = ''),row.names = F)





# 生存分析
if(T){
  ID = uni_cox_df[,1]
  geneIDselect <-select(org.Hs.eg.db, #.db是这个芯片数据对应的注释包
                        keys=ID,
                        columns=c("SYMBOL"), #clolumns参数是你要转换的ID类型是什么，这里选择三个。
                        keytype="ENTREZID" )#函数里面的keytype与keys参数是对应的，keys是你输入的那些数据，keytype是指这些数据是属于什么类型的数据。 
  
  survival_exp = df_expr[ID,]
  rownames(survival_exp) = geneIDselect[,2]
  survival_table = survival
  names = rownames(survival_table)
  sit = which(colnames(survival_exp) %in% names)
  survival_exp = survival_exp[,sit]
  n= 2#  2
  rownames(survival_exp[n,])
  survival_table$ADHFE1 = ifelse(t(survival_exp[n,]) > mean(as.numeric(survival_exp[n,])),'hig','low')
  meta = na.omit(survival_table)
  meta$time = as.numeric(meta$time)
  sfit1 <- survfit(Surv(time, status)~ADHFE1, data=meta) #primary_ER/CTC_ER
  ## more complicate figures.
  ggsurvplot(sfit1,palette =c("#1EB2A6","#F67575"),
             risk.table =TRUE,
             pval =TRUE,
             alpha=0.75,
             conf.int = FALSE,
             xlab ="Time in months", 
             ylab ="Duration of the treatment",
             ncensor.plot = F)
}
