rm(list=ls())
setwd("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")
try(load("RData/DE_GSE179379.RData"))
try(load("RData/Alldata_7Jul.RData"))
homologs<-read.csv("Data/Human rat homologs.txt")
library(clusterProfiler)
library(fgsea)
library(ggplot2)

#######################
### random drugs
######################

load("Results/RData/CMAP_NPC_mean.RData")
load("Results/RData/compounds_NPC.RData")

fgsea_res_drugs<-list()
for(comp in 1:ncol(data_mean_all)){
  print(comp)
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE5388"))
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
  genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  genes_up<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=T)[1:length(genes_up)]]
  genes_dn<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=F)[1:length(genes_dn)]]
  
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
    if(length(which(age>6570))>=10){
      sample<-sample[which(age>6570)]
      data<-data[, sample]
      age<-age[which(age>6570)]
      
      genes_cor[[n]]<-cor(t(data), as.numeric(age))
      dat_names<-c(dat_names, dd)
    }
  }
  
  region<-"PFC"
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  for(dd in dat){
    n<-n+1
    data<-get(dd)
    sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    age<-as.numeric(age)
    if(length(which(age>6570))>=10){
      sample<-sample[which(age>6570)]
      data<-data[, sample]
      age<-age[which(age>6570)]
      
      genes_cor[[n]]<-cor(t(data), as.numeric(age))
      dat_names<-c(dat_names, dd)
    }
  }
  #LSD
  #rat_homologs<-unique(homologs[homologs[,2] %in% rownames(GSE179379)[which(rowSums(GSE179379>10)>1)],3])
  
  genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
  fgsea_res_drugs[[comp]]<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){
    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    if(dat_names[n]=="GTEx_PFC"){
      fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea, nPermSimple=100000)
    } else {
      fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
    }
  }
  
}
save(fgsea_res_drugs, file="RData/fgsea_res_drugs_RAW.RData")
