rm(list=ls())
setwd("workdir")#workdir = working directory

load("Data/Alldata_20Sep.RData")
signature<-read.csv("Data/ACEL-17-e12819-s002.csv")
signature_up<-signature$Gene.Symbol[signature$Microarray=="up"]
signature_dn<-signature$Gene.Symbol[signature$Microarray=="down"]


##################################
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
load("Results/RData/fgsea_res.RData")

fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

###################

library(clusterProfiler)
library(fgsea)
library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741", "GSE5388"))
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
dat<-setdiff(dat, c("GSE102741", "GSE5388"))

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


range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p1<-list()
p2<-list()
df<-list()
fgsea_res<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  fgsea_res[[n]]<-fgsea(list(UP=signature_up, DN=signature_dn), forgesea)
  p1[[n]]<-plotEnrichment(unique(signature_up), forgesea)
  p2[[n]]<-plotEnrichment(unique(signature_dn), forgesea)
  df[[n]]<-rbind.data.frame(data.frame(scaled=range01(p1[[n]]$data[,1]), p1[[n]]$data, path=rep("signature_up", nrow(p1[[n]]$data)), dataset=rep(dat_names[n],nrow(p1[[n]]$data))),
                            data.frame(scaled=range01(p2[[n]]$data[,1]), p2[[n]]$data, path=rep("signature_dn", nrow(p2[[n]]$data)), dataset=rep(dat_names[n],nrow(p2[[n]]$data))))
}
df_signature<-df[[1]]
for(n in 2:length(genes_cor)){
  df_signature<-rbind.data.frame(df_signature, df[[n]])
}

pdf("Results/Figures/GSEA_siganture.pdf",8,8)
ggplot(df_signature, aes(x=rank, y=ES, colour=path))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  facet_wrap(.~dataset)+ scale_colour_manual(values=c("signature_dn"="blue", "signature_up"="red"))
dev.off()



age_up<-read.csv("Results/age_up.csv")
age_dn<-read.csv("Results/age_dn.csv")

meta_dn<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}


p_dn<-c(meta_dn(fgsea_res))
p_up<-c(meta_up(fgsea_res))

