### Title: "生信人"
### 
### Date: "06/22/2020"
### 
### Tips : None
### 
### Best Regards,
### 
### Yuan.SH
### 
### Please contact with me via follow way:
###   
### (a) E-mail :yuansh3354@163.com
### 
### (b) QQ :1044532817
### 
### (c) WeChat :YuanSh181014

setwd('/Users/yuansh/Desktop/shengxinren')

######################################## Step.1_load_packages & functions ########################################

### packages ###
library(TCGAbiolinks)
library(VennDiagram)
library(magrittr)
library("FactoMineR")
library("factoextra")
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(genefu)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
library(stringr)
library(devtools) 
library(DESeq2)
library(survival)
library(survminer)

### functions  ###
# heatmap
heat_map = function(exprSet,group_list,ids = NULL ,names,max_expr = 4 ,show_rownames=F,show_colnames=F){
  library(pheatmap)
  colD=data.frame(group_list=group_list)
  
  if(is.null(ids)){
    ids=names(tail(sort(apply(exprSet,1,sd)),100))
    exprSet = t(scale(t(exprSet[ids,])))
    pdf = '-heatmap.pdf'
  }else{
    exprSet = t(scale(t(exprSet[ids,])))
    pdf = '-ids-heatmap.pdf'
  }
  
  rownames(colD)=colnames(exprSet)
  exprSet[exprSet>max_expr]=max_expr
  exprSet[exprSet< -max_expr]= -max_expr
  pheatmap(exprSet,
           show_rownames=show_rownames,
           show_colnames=show_colnames,
           annotation_col = colD,
           filename = paste(names,pdf,sep = ''))
}

### limma ### 
deg_limma = function(exprSet,design,contrast.matrix){
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

######################################## Step.2_load_datasets ########################################

### 表达谱 ###
expdf = read.table('AgilentG4502A_07_3.txt', header = T, row.names = 1)

### 直接导入列名会出现. 因此需要替换确保列名和clinic的行名一致 ###
t=str_split(colnames(expdf),'\\.',simplify = T) 
colnames(expdf)=paste(t[,1],t[,2],t[,3],t[,4],sep = '-')
expdf = na.omit(expdf)

### 临床信息 ###
clinic = read.csv("COAD_survival.csv",header = T, row.names = 1)
sit = rownames(clinic) %in% colnames(expdf) %>% which
clinic = clinic[sit,]
clinic = clinic[colnames(expdf),]
### grouplist ###
### TCGA命名规则-01 ～-09为疾病样本， -11 为对照样本 ###
group_list = ifelse(as.numeric(substr(colnames(expdf),14,15)) < 10,'tumor','normal')

######################################## Step.3_limma ########################################

design = model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
exprSet=expdf
rownames(design)=colnames(exprSet)
colnames(design)=c('normal','tumor')
contrast.matrix<-makeContrasts("tumor-normal",
                               levels = design)

deg = deg_limma(exprSet,design,contrast.matrix)

######################################## Step.4_volcano ########################################

df=deg
df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
df$g=ifelse(df$P.Value>0.01,'stable',ifelse( df$logFC >2,'up',ifelse( df$logFC < -2,'down','stable') ))
df$abs = abs(df$logFC)
df$name=rownames(df)
top_genes  = rownames(df[order(df$abs,decreasing = T),])[1:10]
ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.8,
          label = "name", repel = T,
          label.select = top_genes, 
          palette = c("#00AFBB", "grey", "#FC4E07") )
ggsave('volcano.png')
sit  = which(df$abs>2 & df$adj.P.Val<0.05)
diff_gene = df[sit,] %>% rownames()
######################################## Step.5_immport ########################################

immport = read.csv('Geneappend3.csv',header = T, stringsAsFactors = F)

### 提取免疫相关差异基因 ###
immune_gene = intersect(immport$Symbol,rownames(deg)) 
ids = intersect(immune_gene,diff_gene)
lt=list(immune_gene = immune_gene,diff_gene = diff_gene)
par(bty="n", mgp = c(2,.33,0), mar=c(3,3,1,0)+.1, las=1, tcl=-.25)
grid.newpage()
pushViewport(viewport(w = .9, h = .9))
grid.draw(venn.diagram(lt, filename = 'Veen.png',
                       fill=c("cornflowerblue","yellow"),
                       alpha=0.85,
                       euler.d=TRUE, 
                       cex = 1.6, 
                       cat.cex = 1.6,
                       cat.fontface="italic", 
                       euler.diagram=TRUE))

######################################## Step.6_heatmap ########################################

### 原始数据聚类 ###
heat_map(exprSet = expdf,group_list = group_list,ids = ids,max_expr = 4,show_rownames = F,names = 'CTL-TUM')
heat_map(exprSet = expdf,group_list = group_list,ids = NULL,max_expr = 4,show_rownames = F,names = 'CTL-TUM')

### 对肿瘤组织聚类 ###
sit = which(group_list == 'tumor')
exp_tumor = expdf[,sit]
tum_list = group_list[sit]
heat_map(exprSet = exp_tumor,group_list = tum_list,ids = ids,max_expr = 3,show_rownames = F,names = 'TUM')
heat_map(exprSet = exp_tumor,group_list = tum_list,ids = NULL,max_expr = 3,show_rownames = F,names = 'TUM')

### 定义肿瘤样本label ###
setwd('/Users/yuansh/Desktop/shengxinren/geneheatmap')

for (i in 1:length(ids)) {
  print(ids[i])
  mean_exp = exp_tumor[ids[i],] %>% as.numeric %>% mean()
  label = ifelse(exp_tumor[ids[i],] < mean_exp,'low_exp','hig_exp') %>% as.character()
  heat_map(exprSet = exp_tumor,group_list = label,ids = NULL,max_expr = 2,show_rownames = F,names = ids[i])
}

setwd('/Users/yuansh/Desktop/shengxinren/')

######################################## Step.7_survival ########################################

### 原始标签生存分析 ###
clinic_survival = data.frame(clinic[,c(2,3)] )  %>% set_colnames(.,c('event','time')) 
clinic_survival[,2] = clinic_survival[,2]/12
clinic_survival$label = group_list
sfit <- survfit(Surv(time, event)~label, data=clinic_survival)
gg = ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =F,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
ggsave('CTL-VS-TUM.tiff',plot = print(gg), width = 18, height = 8.5, units = "in")

### 定义label ###
setwd('/Users/yuansh/Desktop/shengxinren/survival')
sit = which(group_list == 'tumor')

for (i in 1:length(ids)) {
  print(ids[i])
  mean_exp = exp_tumor[ids[i],] %>% as.numeric %>% mean()
  label = ifelse(exp_tumor[ids[i],] < mean_exp,'low_exp','hig_exp') %>% as.character()
  clinic_survival = data.frame(clinic[sit,c(2,3)] )  %>% set_colnames(.,c('event','time')) 
  clinic_survival[,2] = clinic_survival[,2]/12
  clinic_survival$label = label
  sfit <- survfit(Surv(time, event)~label, data=clinic_survival)
  gg = ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
                  risk.table =F,pval =TRUE,
                  conf.int =TRUE,xlab ="Time in months", 
                  ggtheme =theme_light(), 
                  ncensor.plot = TRUE)
  ggsave(paste(ids[i],'survival.tiff',sep = ''),plot = print(gg), width = 4, height = 3, units = "in")
  
  
}

setwd('/Users/yuansh/Desktop/shengxinren/')





