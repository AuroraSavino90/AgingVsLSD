rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")

library(ggplot2)
load("RData/DE_GSE179379.RData")
load("RData/Alldata_7Jul.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

library(pheatmap)

########################
###annotazioni
#######################
df<-data.frame(Age=as.numeric(metadata$Age)/365)
pdf("Results/Figures/Age_frequency.pdf",5,4)
ggplot(df, aes(x=Age))+geom_histogram()+theme_classic()
dev.off()

metadata$Diagnosis<-factor(metadata$Diagnosis, levels=names(sort(table(metadata$Diagnosis), decreasing = T)))
metadata$Region_simpl<-factor(metadata$Region_simpl, levels=names(sort(table(metadata$Region_simpl), decreasing = T)))
metadata$Race<-factor(metadata$Race, levels=names(sort(table(metadata$Race), decreasing = T)))

pdf("Results/Figures/Gender_frequency.pdf",2,5)
ggplot(metadata, aes(x=Gender))+geom_bar()+theme_classic()
dev.off()

pdf("Results/Figures/Diagnosis_frequency.pdf",6,5)
ggplot(metadata, aes(x=Diagnosis))+geom_bar()+theme_classic()
dev.off()

pdf("Results/Figures/Race_frequency.pdf",4,5)
ggplot(metadata, aes(x=Race))+geom_bar()+theme_classic()+ theme (axis.text.x = element_text (angle = 90, vjust = 1, hjust=1))
dev.off()

pdf("Results/Figures/Region_frequency.pdf",5,5)
ggplot(metadata, aes(x=Region_simpl))+geom_bar()+theme_classic()
dev.off()

####Rappresentare il range di et? di ogni dataset
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens"]))
binned<-matrix(nrow=length(dat), ncol=100)
rownames(binned)<-dat
age_median<-c()
for(dd in dat){
    age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens")]))
    age_median<-c(age_median, median(age, na.rm = T))
    for(bin in 1:100){
    age_bin<-c(((bin-1)*365+1):(bin*365))
    if(length(which(age %in% intersect(age_bin, age)))>0){
      binned[dd,bin]<-length(which(age %in% intersect(age_bin, age)))
    } else {
      binned[dd,bin]<-0
    }

  }
}
colnames(binned)<-c(1:100)
binned<-binned[order(age_median),]
pdf("Results/Figures/age_datasets.pdf", 16, 6)
pheatmap(log2(binned+1), cluster_cols = F, cluster_rows = F, cellwidth = 10, cellheight = 15)
dev.off()


binned2<-table(metadata[,c("Diagnosis", "Region_simpl")])
pdf("Results/Figures/samples_combn.pdf", 5, 5)
pheatmap(log2(binned2+1), cluster_cols = F, cluster_rows = F, cellwidth = 15, cellheight = 15)
dev.off()

min(as.numeric(metadata$Age), na.rm = T)
min(as.numeric(metadata$Age)[as.numeric(metadata$Age)>0], na.rm = T)
max(as.numeric(metadata$Age), na.rm = T)/365


#############
### For metadata table
############
table(metadata$Gender)
table(metadata$Diagnosis)

table(metadata$Region_simpl)
table(metadata$Race)
table(metadata$Alcohol)
table(metadata$Smoke)

#improve medication, drugs or remove

#########################################
#########################################
#########################################

age_range<-c()
for(dd in dat){
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd)]))
  age_range<- c(age_range, max(age, na.rm=T)-min(age, na.rm=T))
}
names(age_range)<-dat


cor_age<-function(organism, diagnosis, region){
dat<-na.omit(unique(metadata$Dataset[metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE30272"))#no raw
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance

#compare aging datasets
data<-get(dat[1])
unione<-rownames(data)
for(dd in dat){
  data<-get(dd)
  unione<-union(unione, rownames(data))
  print(length(unione))
}

cor_age<-matrix(nrow=length(unione), ncol=length(dat))
colnames(cor_age)<-dat
rownames(cor_age)<-unione
for(dd in dat){
  data<-get(dd)
  genes<-intersect(rownames(data), unione)
  sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl%in% region)]
  data<-data[, sample]
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd &metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl%in% region)]))
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]
    cor_age[genes,dd]<-cor(t(data[genes,]), age, method="p", use="pairwise.complete.obs")
  }
}
cor_age<-cor_age[,colSums(is.na(cor_age))!=nrow(cor_age)]
return(cor_age)

}

cor_DLPFC<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="Healthy")
cor_mPFC<-cor_age(organism="Homo sapiens", region="PFC", diagnosis="Healthy")

colnames(cor_DLPFC)<-paste(colnames(cor_DLPFC), "DLPFC")
colnames(cor_mPFC)<-paste(colnames(cor_mPFC), "mPFC")

allgenes<-unique(c(rownames(cor_DLPFC), rownames(cor_mPFC)))
nsets<-ncol(cor_DLPFC)+ncol(cor_mPFC)
cor_age_tot<-matrix(nrow=length(allgenes), ncol=nsets)
colnames(cor_age_tot)<-c(colnames(cor_DLPFC), colnames(cor_mPFC))
rownames(cor_age_tot)<-allgenes

cor_age_tot[rownames(cor_DLPFC), colnames(cor_DLPFC)]<-cor_DLPFC
cor_age_tot[rownames(cor_mPFC), colnames(cor_mPFC)]<-cor_mPFC

hist(rowSums(cor_age_tot>0, na.rm=T))
hist(rowSums(cor_age_tot<0, na.rm=T))

age_up<-rownames(cor_age_tot)[rowSums(cor_age_tot>0, na.rm=T)>=12]
age_dn<-rownames(cor_age_tot)[rowSums(cor_age_tot<0, na.rm=T)>=12]

write.csv(age_up, file="age_up.csv")
write.csv(age_dn, file="age_dn.csv")

toplot<-cor_age_tot[age_up,]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(t(toplot),   cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T)

toplot<-cor_age_tot[age_dn,]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(t(toplot),   cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T)

cor_datasets<-cor(cor_age_tot, method="s", use="pairwise.complete.obs")
hist(cor_datasets[upper.tri(cor_datasets)], breaks=20)


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_datasets), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(cor_datasets), na.rm=T)/paletteLength, max(unlist(cor_datasets), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

p<-pheatmap(cor_datasets, cellwidth=15, cellheight=15, breaks=myBreaks,
            color = myColor)

pdf("Results/Figures/Aging_similarity.pdf",8,8)
print(p)
dev.off()



cor_DLPFC<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="Healthy")
cor_mPFC<-cor_age(organism="Homo sapiens", region="PFC", diagnosis="Healthy")
cor_TC<-cor_age(organism="Homo sapiens", region="TC", diagnosis="Healthy")
cor_OPFC<-cor_age(organism="Homo sapiens", region="OPFC", diagnosis="Healthy")

colnames(cor_DLPFC)<-paste(colnames(cor_DLPFC), "DLPFC")
colnames(cor_mPFC)<-paste(colnames(cor_mPFC), "mPFC")
colnames(cor_OPFC)<-paste(colnames(cor_OPFC), "OPFC")

allgenes<-unique(c(rownames(cor_DLPFC), rownames(cor_mPFC), names(cor_TC), rownames(cor_OPFC), names(cor_SFG)))
nsets<-ncol(cor_DLPFC)+ncol(cor_mPFC)+1+ncol(cor_OPFC)
cor_age_tot<-matrix(nrow=length(allgenes), ncol=nsets)
colnames(cor_age_tot)<-c(colnames(cor_DLPFC), colnames(cor_mPFC), "TC", colnames(cor_OPFC))
rownames(cor_age_tot)<-allgenes

cor_age_tot[rownames(cor_DLPFC), colnames(cor_DLPFC)]<-cor_DLPFC
cor_age_tot[rownames(cor_mPFC), colnames(cor_mPFC)]<-cor_mPFC
cor_age_tot[names(cor_TC), "TC"]<-cor_TC
cor_age_tot[rownames(cor_OPFC), colnames(cor_OPFC)]<-cor_OPFC



cor_datasets<-cor(cor_age_tot, method="s", use="pairwise.complete.obs")
hist(cor_datasets[upper.tri(cor_datasets)], breaks=20)


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_datasets), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(cor_datasets), na.rm=T)/paletteLength, max(unlist(cor_datasets), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(cor_datasets, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)

#########################
##### Disease
#######################

cor_DLPFC_s<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="Schizophrenia")
cor_mPFC_s<-cor_age(organism="Homo sapiens", region="PFC", diagnosis="Schizophrenia")
cor_DLPFC_m<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="MDD")
cor_DLPFC_b<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="Bipolar")

colnames(cor_DLPFC_m)<-paste(colnames(cor_DLPFC_m), "m")
colnames(cor_DLPFC_b)<-paste(colnames(cor_DLPFC_b), "b")

allgenes<-unique(c(rownames(cor_DLPFC_b), rownames(cor_DLPFC_m), names(cor_DLPFC_s), names(cor_mPFC_s)))
nsets<-ncol(cor_DLPFC_b)+ncol(cor_DLPFC_m)+2
cor_age_tot_dis<-matrix(nrow=length(allgenes), ncol=nsets)
colnames(cor_age_tot_dis)<-c(colnames(cor_DLPFC_b), colnames(cor_DLPFC_m), "DLPFC_s", "mPFC_s")
rownames(cor_age_tot_dis)<-allgenes

cor_age_tot_dis[rownames(cor_DLPFC_b), colnames(cor_DLPFC_b)]<-cor_DLPFC_b
cor_age_tot_dis[rownames(cor_DLPFC_m), colnames(cor_DLPFC_m)]<-cor_DLPFC_m
cor_age_tot_dis[names(cor_DLPFC_s), "DLPFC_s"]<-cor_DLPFC_s
cor_age_tot_dis[names(cor_mPFC_s), "mPFC_s"]<-cor_mPFC_s


cor_datasets<-cor(cor_age_tot_dis, method="s", use="pairwise.complete.obs")
hist(cor_datasets[upper.tri(cor_datasets)], breaks=20)


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_datasets), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(cor_datasets), na.rm=T)/paletteLength, max(unlist(cor_datasets), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(cor_datasets, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, clustering_distance_rows = "correlation",  clustering_distance_cols = "correlation")

incommon<-intersect(rownames(cor_age_tot), rownames(cor_age_tot_dis))
cor_age_tot_tot<-cbind(cor_age_tot[incommon, 1:13], cor_age_tot_dis[incommon,-8])


cor_datasets<-cor(cor_age_tot_tot, method="s", use="pairwise.complete.obs")


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_datasets), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(cor_datasets), na.rm=T)/paletteLength, max(unlist(cor_datasets), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(cor_datasets, cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor,
         clustering_distance_rows = "correlation",  clustering_distance_cols = "correlation")

pheatmap(cor_datasets, cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor)
