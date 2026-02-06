rm(list=ls())
setwd("workdir")#workdir = working directory

library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
load("Data/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier

dat_names<-c()
age_all<-c()
dataset_all<-c()
for(dd in dat){

  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]
   dat_names<-c(dat_names, dd)
   age_all<-c(age_all, age)
   dataset_all<-c(dataset_all, rep(dd, ncol(data)))
  }
}

save(age_all, file="Results/RData/age_all.RData")


all<-c()
for(dd in dat_names){
  data<-get(dd)
  all<-c(all, rownames(data))
}

sum(table(all)>=11)
selgenes<-names(which(table(all)>=11))

############################ALLDATA
data<-get(dat_names[1])
sample<-metadata$Sample[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
data<-data[, sample]
age<-metadata$Age[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
age<-as.numeric(age)
if(length(which(age>=7300))>=10){
  sample<-sample[which(age>=7300)]
  data<-data[, sample]
alldata<-data[selgenes,]
}

for(i in 2:length(dat_names)){
  data<-get(dat_names[i])
  sample<-metadata$Sample[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>=7300))>=10){
    sample<-sample[which(age>=7300)]
    data<-data[, sample]
    age<-age[which(age>=7300)]


  alldata<-cbind(alldata, data[selgenes,])
  }
}

library(FactoMineR)
pca<-PCA(t(alldata))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], age=age_all/365, dataset=dataset_all)

pdf("Results/Figures/PCA_DLPF_dataset.pdf",7,5)
ggplot(df, aes(PC2, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_DLPF_age.pdf",6,5)
ggplot(df, aes(PC2, PC1, colour=age))+geom_point(size=1)+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()
pdf("Results/Figures/PC1_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC1, PC4, colour=dataset))+geom_point()+theme_classic()
ggplot(df, aes(PC1, PC4, colour=age))+geom_point()+theme_classic()


batch = dataset_all

# parametric adjustment
library(sva)
mod = model.matrix(~age_all, data=data.frame(alldata))
mod0 = model.matrix(~1,data=data.frame(t(alldata)))
n.sv = num.sv(data.frame(alldata),mod,method="leek")
svobj = sva(as.matrix(alldata),mod,mod0,n.sv=n.sv)
plot(age_all, svobj[[1]][,1], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,2], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,3], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,4], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,5], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,6], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,7], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,8], col=as.factor(dataset_all))

alldata_resid<-apply(alldata,1, function(x){residuals(lm(unlist(x)~svobj[[1]][,1:n.sv]))})
#alldata_resid<-apply(alldata,1, function(x){residuals(lm(unlist(x)~svobj[[1]][,1]+svobj[[1]][,2]+svobj[[1]][,3]+svobj[[1]][,4]+svobj[[1]][,5]+svobj[[1]][,6]+svobj[[1]][,7]))})
alldata_resid<-t(alldata_resid)

pca<-PCA(t(alldata_resid))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
pdf("Results/Figures/PCA_SVA_dataset.pdf",7,5)
ggplot(df, aes(PC2, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_SVA_age.pdf",6,5)
ggplot(df, aes(PC2, PC1, colour=age))+geom_point()+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()

ggplot(df, aes(PC3, PC4, colour=age))+geom_point()

pdf("Results/Figures/PC1_SVA_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()
cor.test(df$age, df$PC1)

save(alldata_resid, file="Results/RData/alldata_resid.RData")

#################################################
#########GSEA with the batch corrected dataset
##################################################

library(fgsea)
library(ggplot2)

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

genes_cor<-cor(t(alldata_resid), as.numeric(age_all))

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

  forgesea<-unlist(genes_cor)
  names(forgesea)<-rownames(genes_cor)
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), 
                       data.frame(p2$data, dir=rep("dn", nrow(p2$data))))
  
  p<-ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf("Results/Figures/GSEA_all_batchcorrect.pdf", 5,4)
ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))
dev.off()

