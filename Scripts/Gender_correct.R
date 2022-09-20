rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")
load("RData/Alldata_27May.RData")

##################################
load("RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

library(fgsea)
library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE30272"))#tmp
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
dat<-setdiff(dat, c("GSE5388"))#outlier

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

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
  gender<-metadata$Gender[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    gender<-gender[which(age>6570)]
    age<-age[which(age>6570)]
    
    if(length(table(gender))>1){
      data_resid<-apply(data,1, function(x){residuals(lm(unlist(x)~gender))})
      data_resid<-t(data_resid)
      
    genes_cor[[n]]<-cor(t(data_resid), as.numeric(age))
    dat_names<-c(dat_names, dd)
    }
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
  gender<-metadata$Gender[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    gender<-gender[which(age>6570)]
    age<-age[which(age>6570)]
    
    if(length(table(gender))>1){
      data_resid<-apply(data,1, function(x){residuals(lm(unlist(x)~gender))})
      data_resid<-t(data_resid)
      
      genes_cor[[n]]<-cor(t(data_resid), as.numeric(age))
      dat_names<-c(dat_names, dd)
    }
  }
}
#LSD
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_res<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
  
  p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
  
}


df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}


fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/GSEA_all_datasets_gender.pdf", 8, 8)
ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x")
dev.off()


library(metap)

meta_dn<-function(x){
  istwo <- rep(T, 10)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up<-function(x){
  istwo <- rep(T, 10)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

pdn<-meta_dn(fgsea_res)
pup<-meta_up(fgsea_res)

