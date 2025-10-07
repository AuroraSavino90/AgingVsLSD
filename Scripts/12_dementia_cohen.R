rm(list=ls())
setwd("workdir")#workdir = working directory

load("Data/Psychedelics_PFC.RData")
load("Data/Alldata_20Sep.RData")
load("Data/ExerciseAndAlcohol.RData")
load("Data/Dementia_alldata_meta.RData")
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
library(pheatmap)
library(ggplot2)

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

region<-"DLPFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))

library(GSVA)
cohen_up<-c()
cohen_dn<-c()
p_up<-c()
p_dn<-c()
dat_names<-c()
for(dd in dat){
  data<-get(dd)
  if(sum(is.na(data))>0){
  data<- data[-which(is.na(data), arr.ind = T)[,1],]
  }
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
  ## build GSVA parameter object
  gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
  
  ## estimate GSVA enrichment scores for the three sets
  ssgsea <- gsva(gsvapar)

    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
      dat_names<-c(dat_names, dd)
      p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
      p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])

  }

region<-"PFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))


for(dd in dat){
  data<-get(dd)
  if(sum(is.na(data))>0){
    data<- data[-which(is.na(data), arr.ind = T)[,1],]
  }
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))

  ## build GSVA parameter object
  gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
  
  ## estimate GSVA enrichment scores for the three sets
  ssgsea <- gsva(gsvapar)

  cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
  cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
  p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
  p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])
  dat_names<-c(dat_names, dd)
}


region<-"PFC"
diagnosis<-c("Huntington")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","Huntington")



for(dd in dat){
  data<-get(dd)
  if(sum(is.na(data))>0){
    data<- data[-which(is.na(data), arr.ind = T)[,1],]
  }
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Huntington"))

  ## build GSVA parameter object
  gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
  
  ## estimate GSVA enrichment scores for the three sets
  ssgsea <- gsva(gsvapar)

  cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Huntington"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Huntington"])^2)/2))
  cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Huntington"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Huntington"])^2)/2))
  p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Huntington"])[[3]])
  p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Huntington"])[[3]])
  dat_names<-c(dat_names, dd)
}


region<-"DLPFC"
diagnosis<-c("PDD")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","PDD")


for(dd in dat){
  data<-get(dd)
  if(sum(is.na(data))>0){
    data<- data[-which(is.na(data), arr.ind = T)[,1],]
  }
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "PDD"))

  ## build GSVA parameter object
  gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
  
  ## estimate GSVA enrichment scores for the three sets
  ssgsea <- gsva(gsvapar)

  cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="PDD"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="PDD"])^2)/2))
  cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="PDD"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="PDD"])^2)/2))
  p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="PDD"])[[3]])
  p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="PDD"])[[3]])
  dat_names<-c(dat_names, dd)
}


region<-"DLPFC"
diagnosis<-c("DLB")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
diagnosis<-c("Healthy","DLB")


for(dd in dat){
  data<-get(dd)
  if(sum(is.na(data))>0){
    data<- data[-which(is.na(data), arr.ind = T)[,1],]
  }
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  age<-as.numeric(age)
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "DLB"))

  ## build GSVA parameter object
  gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
  
  ## estimate GSVA enrichment scores for the three sets
  ssgsea <- gsva(gsvapar)

  cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="DLB"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="DLB"])^2)/2))
  cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="DLB"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="DLB"])^2)/2))
  p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="DLB"])[[3]])
  p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="DLB"])[[3]])
  dat_names<-c(dat_names, dd)
}


forheat<-rbind(cohen_up, cohen_dn)
colnames(forheat)<-dat_names
rownames(forheat)<-c("up with LSD", "down with LSD")
#forheat[ptot>0.1]<-NA
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(forheat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(forheat), na.rm=T)/paletteLength, max(unlist(forheat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(forheat[,order(forheat[1,], decreasing = T)], cluster_cols = F,  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)


df<-data.frame(cohen=c(cohen_up, cohen_dn), direction=rep(c("up with drug", "dn with drug"), each=length(cohen_up)),
               pval= -log10(c(p_up, p_dn)),
               disease=rep(c("AD", "AD", "AD", "AD", "AD", "HD", "PDD", "DLB"),2))

gg<-ggplot(df, aes(x=cohen, y=pval, colour=direction))+geom_point()+theme_classic()+ labs(x = "Cohen's d", y = "-log10(pvalue)")+
  geom_vline(xintercept=0, linetype='dashed', color='grey', size=0.5)+theme_classic()+scale_colour_manual(values=c("dn with drug"="blue", "up with drug"="red"))

pdf("Results/Figures/dementia_LSD.pdf",5,5)
print(gg)
dev.off()


library(metap)

istwo <- rep(T, length(cohen_up))
toinvert <- ifelse(cohen_up<0, T, F)

p_up_LSD<-sumlog(two2one(p_up, two = istwo, invert = toinvert))


istwo <- rep(T, length(cohen_dn))
toinvert <- ifelse(cohen_dn>0, T, F)

p_dn_LSD<-sumlog(two2one(p_dn, two = istwo, invert = toinvert))

istwo <- rep(T, length(cohen_up))
toinvert <- ifelse(cohen_up>0, T, F)

p_up_LSD_inv<-sumlog(two2one(p_up, two = istwo, invert = toinvert))


istwo <- rep(T, length(cohen_dn))
toinvert <- ifelse(cohen_dn<0, T, F)

p_dn_LSD_inv<-sumlog(two2one(p_dn, two = istwo, invert = toinvert))

#########
##Other psychedelics
###########
rat_homologs<-read.csv("Data/Human rat homologs.txt")
mouse_homologs<-read.csv("Data/Human mouse homologs.txt", sep="\t")

mouse_data<-c("DE_GSE64607","DE_GSE164798_chronic", "DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57", "DE_MDMA", "DE_GSE161626",
              "DE_GSE161626_48h","DE_GSE161626_7d","DE_GSE81672","DE_GSE26364","DE_GSE209859",
              "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
              "DE_GSE105453_Y", "DE_GSE105453_O", "DE_GSE111273_Y", "DE_GSE111273_O")
rat_data<-c("DE_GSE14720", "DE_GSE23728","DE_GSE179380","DE_DMT", "DE_pharm", "DE_harm")
##MDMA
DEGs_array<-c("DE_MDMA",
              "DE_GSE26364",
              #DOI
              "DE_GSE23728",
              #LSD
              "DE_GSE179380", 
              "DE_GSE64607","DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57"
              
)

DEGs_seq<-c("DE_DMT", "DE_GSE161626","DE_GSE161626_48h","DE_GSE161626_7d", "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
            "DE_GSE81672", "DE_harm", "DE_pharm", "DE_GSE164798_chronic", "DE_GSE111273_Y", "DE_GSE111273_O")

p_up_tot<-c()
p_dn_tot<-c()
p_up_tot_inv<-c()
p_dn_tot_inv<-c()

for(DE in DEGs_array){
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }

  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  
  region<-"DLPFC"
  diagnosis<-c("Healthy", "Alzheimer")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))
  
  library(GSVA)
  cohen_up<-c()
  cohen_dn<-c()
  p_up<-c()
  p_dn<-c()
  dat_names<-c()
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
    dat_names<-c(dat_names, dd)
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])
    
  }
  
  region<-"PFC"
  diagnosis<-c("Healthy", "Alzheimer")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"PFC"
  diagnosis<-c("Huntington")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","Huntington")
  
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Huntington"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Huntington"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Huntington"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Huntington"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Huntington"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Huntington"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Huntington"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"DLPFC"
  diagnosis<-c("PDD")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","PDD")
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "PDD"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="PDD"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="PDD"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="PDD"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="PDD"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="PDD"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="PDD"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"DLPFC"
  diagnosis<-c("DLB")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","DLB")
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "DLB"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="DLB"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="DLB"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="DLB"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="DLB"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="DLB"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="DLB"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  istwo <- rep(T, length(cohen_up))
  toinvert <- ifelse(cohen_up<0, T, F)
  
  p_up_tot<-c(p_up_tot, sumlog(two2one(p_up, two = istwo, invert = toinvert))$p)
  
  
  istwo <- rep(T, length(cohen_dn))
  toinvert <- ifelse(cohen_dn>0, T, F)
  
  p_dn_tot<-c(p_dn_tot, sumlog(two2one(p_dn, two = istwo, invert = toinvert))$p)
  
  istwo <- rep(T, length(cohen_up))
  toinvert <- ifelse(cohen_up>0, T, F)
  
  p_up_tot_inv<-c(p_up_tot_inv, sumlog(two2one(p_up, two = istwo, invert = toinvert))$p)
  
  
  istwo <- rep(T, length(cohen_dn))
  toinvert <- ifelse(cohen_dn<0, T, F)
  
  p_dn_tot_inv<-c(p_dn_tot_inv, sumlog(two2one(p_dn, two = istwo, invert = toinvert))$p)
  
}



for(DE in DEGs_seq){
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  
  region<-"DLPFC"
  diagnosis<-c("Healthy", "Alzheimer")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))
  
  library(GSVA)
  cohen_up<-c()
  cohen_dn<-c()
  p_up<-c()
  p_dn<-c()
  dat_names<-c()
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
    dat_names<-c(dat_names, dd)
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])
    
  }
  
  region<-"PFC"
  diagnosis<-c("Healthy", "Alzheimer")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Alzheimer"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Alzheimer"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Alzheimer"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Alzheimer"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Alzheimer"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Alzheimer"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"PFC"
  diagnosis<-c("Huntington")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","Huntington")
  
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "Huntington"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="Huntington"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="Huntington"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="Huntington"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="Huntington"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="Huntington"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="Huntington"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"DLPFC"
  diagnosis<-c("PDD")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","PDD")
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "PDD"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="PDD"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="PDD"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="PDD"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="PDD"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="PDD"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="PDD"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  
  region<-"DLPFC"
  diagnosis<-c("DLB")
  dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis==diagnosis & metadata$Region %in% region]))
  diagnosis<-c("Healthy","DLB")
  
  
  for(dd in dat){
    data<-get(dd)
    if(sum(is.na(data))>0){
      data<- data[-which(is.na(data), arr.ind = T)[,1],]
    }
    sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    age<-as.numeric(age)
    disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
    disease<-factor(disease, levels=c("Healthy", "DLB"))
    
    ## build GSVA parameter object
    gsvapar <- ssgseaParam(as.matrix(data), geneSets=list(up=intersect(rownames(data), genes_up),dn=intersect(rownames(data), genes_dn)))
    
    ## estimate GSVA enrichment scores for the three sets
    ssgsea <- gsva(gsvapar)
    
    cohen_up<-c(cohen_up, (mean(ssgsea[1,disease=="Healthy"])-mean(ssgsea[1,disease=="DLB"]))/sqrt((sd(ssgsea[1,disease=="Healthy"])^2+sd(ssgsea[1,disease=="DLB"])^2)/2))
    cohen_dn<-c(cohen_dn, (mean(ssgsea[2,disease=="Healthy"])-mean(ssgsea[2,disease=="DLB"]))/sqrt((sd(ssgsea[2,disease=="Healthy"])^2+sd(ssgsea[2,disease=="DLB"])^2)/2))
    p_up<-c(p_up, t.test(ssgsea[1,disease=="Healthy"],ssgsea[1,disease=="DLB"])[[3]])
    p_dn<-c(p_dn, t.test(ssgsea[2,disease=="Healthy"],ssgsea[2,disease=="DLB"])[[3]])
    dat_names<-c(dat_names, dd)
  }
  
  istwo <- rep(T, length(cohen_up))
  toinvert <- ifelse(cohen_up<0, T, F)
  
  p_up_tot<-c(p_up_tot, sumlog(two2one(p_up, two = istwo, invert = toinvert))$p)
  
  
  istwo <- rep(T, length(cohen_dn))
  toinvert <- ifelse(cohen_dn>0, T, F)
  
  p_dn_tot<-c(p_dn_tot, sumlog(two2one(p_dn, two = istwo, invert = toinvert))$p)
  
  istwo <- rep(T, length(cohen_up))
  toinvert <- ifelse(cohen_up>0, T, F)
  
  p_up_tot_inv<-c(p_up_tot_inv, sumlog(two2one(p_up, two = istwo, invert = toinvert))$p)
  
  
  istwo <- rep(T, length(cohen_dn))
  toinvert <- ifelse(cohen_dn<0, T, F)
  
  p_dn_tot_inv<-c(p_dn_tot_inv, sumlog(two2one(p_dn, two = istwo, invert = toinvert))$p)
  
}

library(ggrepel)
df<-data.frame(pup=-log10(c(p_up_LSD$p, unlist(p_up_tot))), pdn=-log10(c(p_dn_LSD$p,unlist(p_dn_tot))), 
               dataset=c("LSD", DEGs_array, DEGs_seq))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")



dataset<-c("MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise (28d)", "Alcohol (vapour chamber)", "Alcohol (vapour chamber 2)", "Alcohol (IP injection)",
           "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "Exercise (30d)","EE (Young)", "EE (Old)", 
           "LSD (chronic)")
type<-rep("Psychoplastogen", length(dataset))
type[c(5, 18:20)]<-"Positive Control"
type[c(6,7,8)]<-"Negative Control"

df<-data.frame(pup=c(-log10(c(p_up_LSD$p, unlist(p_up_tot)))+log10(c(p_up_LSD_inv$p, unlist(p_up_tot_inv)))), 
               pdn=c(-log10(c(p_dn_LSD$p,unlist(p_dn_tot)))+log10(c(p_dn_LSD_inv$p, unlist(p_dn_tot_inv)))), dataset=dataset, direction="increase aging",
               type=type)

pdf("Results/Figures/All_psychedelicsAndCTRL_dementia.pdf",6.5,5.5)
ggplot(df, aes(x=pup, y=pdn, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()




