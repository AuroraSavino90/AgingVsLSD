rm(list=ls())
load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Psychedelics metanalysis/Psychedelics_PFC.RData")
load("Results/RData/Alldata_20Sep.RData")
load("Results/RData/ExerciseAndAlcohol.RData")
load("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Dementia/Dementia_alldata.RData")
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
library(pheatmap)
source("Scripts/Dementia_meta.R")

DExpr<-function(cond1, cond2){
  cond1<-data.frame(cond1)
  cond2<-data.frame(cond2)
  gs <- factor(c(rep("Cond1",ncol(cond1)),rep("Cond2",ncol(cond2))))
  #}
  design <- model.matrix(~gs + 0, cbind.data.frame(cond1, cond2))
  #colnames(design) <- levels(gs)
  
  fit <- lmFit(cbind.data.frame(cond1, cond2), design)  # fit linear model
  
  # set up contrasts of interest and recalculate model coefficients
  #cts <- paste(levels(gs)[1], levels(gs)[2], sep="-")
  cts<-"gsCond1-gsCond2"
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  return(tT)
}

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

##########################
###### MSIGDB
###########################
library(msigdbr)
library(clusterProfiler)
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy" ) %>% 
  dplyr::select(gs_name, gene_symbol)

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

DEG_down_conv<-unique(homologs[which(homologs[,2] %in% DEG_down),3])
DEG_up_conv<-unique(homologs[which(homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down_LSD<- enricher(gene = DEG_down_conv,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe=universe, 
                   TERM2GENE = m_df)

TF_up_LSD<- enricher(gene = DEG_up_conv,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 universe=universe, TERM2GENE = m_df)



region<-"DLPFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))

TF_down<-list()
TF_up<-list()
i<-0
dat_names<-c()
for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))

  cond2<-data[,disease=="Healthy"]#ctrl
  cond1<-data[,disease=="Alzheimer"]#passive administration
  
  DE<-DExpr(cond1, cond2)
  
  genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
  genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
   
  TF_up[[i]]<- enricher(gene = genes_up,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     universe=rownames(DE), 
                     TERM2GENE = m_df)
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=rownames(DE), 
                   TERM2GENE = m_df)
  
  dat_names<-c(dat_names, dd)
  
  }

region<-"PFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))


for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))

  cond2<-data[,disease=="Healthy"]#ctrl
  cond1<-data[,disease=="Alzheimer"]#passive administration
  
  DE<-DExpr(cond1, cond2)
  genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
  genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
  
  TF_up[[i]]<- enricher(gene = genes_up,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     universe=rownames(DE), 
                     TERM2GENE = m_df)
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=rownames(DE), 
                   TERM2GENE = m_df)
  
  dat_names<-c(dat_names, dd)
}


region<-"PFC"
diagnosis<-c("Huntington")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","Huntington")



for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Huntington"))

  cond2<-data[,disease=="Healthy"]#ctrl
  cond1<-data[,disease=="Huntington"]#passive administration
  
  DE<-DExpr(cond1, cond2)
  
  genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
  genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
  
  TF_up[[i]]<- enricher(gene = genes_up,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     universe=rownames(DE), 
                     TERM2GENE = m_df)
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=rownames(DE), 
                   TERM2GENE = m_df)
  
  dat_names<-c(dat_names, dd)
}


region<-"DLPFC"
diagnosis<-c("PDD")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","PDD")


for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "PDD"))

  cond2<-data[,disease=="Healthy"]#ctrl
  cond1<-data[,disease=="PDD"]#passive administration
  
  DE<-DExpr(cond1, cond2)
  
  genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
  genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
  
  TF_up[[i]]<- enricher(gene = genes_up,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     universe=rownames(DE), 
                     TERM2GENE = m_df)
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=rownames(DE), 
                   TERM2GENE = m_df)
  
  dat_names<-c(dat_names, dd)
}


region<-"DLPFC"
diagnosis<-c("DLB")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","DLB")


for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "DLB"))

  cond2<-data[,disease=="Healthy"]#ctrl
  cond1<-data[,disease=="DLB"]#passive administration
  
  DE<-DExpr(cond1, cond2)
  
  genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
  genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
  
  TF_up[[i]]<- enricher(gene = genes_up,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1,
                     universe=rownames(DE), 
                     TERM2GENE = m_df)
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=rownames(DE), 
                   TERM2GENE = m_df)
  
  dat_names<-c(dat_names, dd)
}


######evaluation collapsing the pvalue

library(metap)
allpaths<-c()
for(i in 1:length(TF_down)){
  allpaths<-union(allpaths, TF_down[[i]]$ID)
}

allp_mat<-matrix(NA,nrow=length(allpaths), ncol=length(TF_down))
rownames(allp_mat)<-allpaths
for(i in 1:length(TF_down)){
  allp_mat[TF_down[[i]]$ID,i]<-TF_down[[i]]$pvalue
}


allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% TF_up_LSD$ID),]
pup<-list()
ind<-0
for(up in TF_up_LSD$ID[which(TF_up_LSD$ID %in% rownames(allp_mat_sel))]){
  ind<-ind+1
   istwo <- rep(F, length(genes_cor))
  toinvert <- rep(F,  length(genes_cor))
  
    pup[[ind]]<-sumlog(two2one(na.omit(allp_mat_sel[up, ]), two = istwo, invert = toinvert))[[3]]
 }
names(pup)<-TF_up_LSD$ID[which(TF_up_LSD$ID %in% rownames(allp_mat_sel))]
pup_dementia<-pup
allp_mat_sel_dementia<-allp_mat_sel
save(pup_dementia, file="Results/RData/pup_dementia_TF.RData")
save(allp_mat_sel_dementia, file="Results/RData/allp_dementia_TF.RData")




