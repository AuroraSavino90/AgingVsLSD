library(FactoMineR)
load("RData/Alldata_7Jul.RData")

region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))

pc1<-list()
pc2<-list()
cor1<-list()
cor2<-list()
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
    pca<-PCA(t(data), graph = FALSE)
    pc1[[n]]<-pca$eig[1,2]
    pc2[[n]]<-pca$eig[2,2]
    cor1[[n]]<-cor.test(pca$ind$coord[,1], age)[[3]]
    cor2[[n]]<-cor.test(pca$ind$coord[,2], age)[[3]]
    
    dat_names<-c(dat_names, dd)
  }
}

region<-"PFC"
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))

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
    pca<-PCA(t(data), graph = FALSE)
    pc1[[n]]<-pca$eig[1,2]
    pc2[[n]]<-pca$eig[2,2]
    cor1[[n]]<-cor.test(pca$ind$coord[,1], age)[[3]]
    cor2[[n]]<-cor.test(pca$ind$coord[,2], age)[[3]]
    dat_names<-c(dat_names, dd)
  }
}

library(ggplot2)
library(ggrepel)
df<-data.frame(pc1=unlist(pc1), cor1=-log10(unlist(cor1)), dataset=dat_names)
ggplot(df, aes(x=pc1, y=cor1, label=dataset))+geom_point()+geom_label_repel()
plot(unlist(pc1), unlist(cor1))
