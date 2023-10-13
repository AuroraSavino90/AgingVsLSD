rm(list=ls())
setwd("workdir")#workdir = working directory

load("Results/RData/DE_GSE179379.RData")
load("Data/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

##########################
#### GSEA GO
##########################

library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(org.Hs.eg.db)

region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier
dat_names<-c()
genes_cor<-list()
n<-0
for(dd in dat){
  n<-n+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]
    
    genes_cor[[n]]<-cor(t(data), as.numeric(age))
    dat_names<-c(dat_names, dd)
  }
}

region<-"PFC"
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier

for(dd in dat){
  n<-n+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]
    
    genes_cor[[n]]<-cor(t(data), as.numeric(age))
    dat_names<-c(dat_names, dd)
  }
}

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]

save(genes_cor, file="Results/genes_cor.RData")

gsea_GO<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  
  gsea_GO[[n]]<-gseGO(forgesea,
                      ont = "BP",
                      "org.Hs.eg.db",
                      keyType = "SYMBOL",
                      exponent = 1,
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 1e-10,
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      verbose = TRUE,
                      seed = FALSE,
                      by = "fgsea"
  )
  
}
names(gsea_GO)<-dat_names

save(gsea_GO, file="Results/RData/gsea_GO.RData")
