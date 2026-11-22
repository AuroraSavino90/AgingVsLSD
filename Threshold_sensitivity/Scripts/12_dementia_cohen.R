# ============================================================
# THRESHOLD-SENSITIVITY COPY of Scripts/12_dementia_cohen.R
# Changes vs the original: (1) gene selection threshold set by
# SELECTION_MODE for ALL datasets (P_ARRAY / P_SEQ columns);
# (2) outputs written to OUT (Threshold_sensitivity/Results_<mode>);
# (3) see sensitivity_setup.R. Run via Threshold_sensitivity/run_all.R
# ============================================================
rm(list=ls())
source("Threshold_sensitivity/Scripts/sensitivity_setup.R")  # SELECTION_MODE, OUT, P_ARRAY/P_SEQ, shims

load("Data/original_RData/Psychedelics_PFC.RData")
load("Data/original_RData/Alldata_20Sep.RData")
load("Data/original_RData/ExerciseAndAlcohol.RData")
load("Data/original_RData/Dementia_alldata_meta.RData")
load("Data/intermediate_RData/DE_GSE179379.RData")
homologs<-read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022
library(pheatmap)
library(ggplot2)

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379[[P_SEQ]]<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379[[P_SEQ]]<0.05)]
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
p_up_list<-list()
p_dn_list<-list()
cohen_up_list<-list()
cohen_dn_list<-list()

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

p_up_list[["DE_GSE179379"]]<-p_up
p_dn_list[["DE_GSE179379"]]<-p_dn
cohen_up_list[["DE_GSE179379"]]<-cohen_up
cohen_dn_list[["DE_GSE179379"]]<-cohen_dn


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

pdf(file.path(OUT, "Figures/dementia_LSD.pdf"),5,5)
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
rat_homologs<-read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022
mouse_homologs<-read.csv("Data/homologs/Human mouse homologs.txt", sep="\t")  # Ensembl 107: BioMart export, human genes (GRCh38.p13) with mouse homologs, August 2022

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
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs[[P_ARRAY]]<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs[[P_ARRAY]]<0.05)]
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
  
  p_up_list[[DE]]<-p_up
  p_dn_list[[DE]]<-p_dn
  cohen_up_list[[DE]]<-cohen_up
  cohen_dn_list[[DE]]<-cohen_dn
  
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
  condition<-c()
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs[[P_SEQ]]<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs[[P_SEQ]]<0.05)]
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
    condition<-c(condition, "Alzheimer DLPFC")
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
    condition<-c(condition, "Alzheimer PFC")
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
    condition<-c(condition, "Huntington PFC")
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
    condition<-c(condition, "PDD DLPFC")
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
    condition<-c(condition, "DLB DLPFC")
  }
  
  p_up_list[[DE]]<-p_up
  p_dn_list[[DE]]<-p_dn
  cohen_up_list[[DE]]<-cohen_up
  cohen_dn_list[[DE]]<-cohen_dn
  
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



dataset<-c("LSD (chronic)", "MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise (28d)", "Alcohol (vapour chamber)", "Alcohol (vapour chamber 2)", "Alcohol (IP injection)",
           "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "Exercise (30d)","EE (Young)", "EE (Old)"
           )
type<-rep("Psychoplastogen", length(dataset))
type[c(6, 19:21)]<-"Positive Control"
type[c(7,8,9)]<-"Negative Control"

df_all<-data.frame(pup_rev= -log10(c(p_up_LSD$p,unlist(p_up_tot))),
                   pup_mimick= -log10(c(p_up_LSD_inv$p, unlist(p_up_tot_inv))),
pdn_rev= -log10(c(p_dn_LSD$p,unlist(p_dn_tot))),
pdn_mimick= -log10(c(p_dn_LSD_inv$p, unlist(p_dn_tot_inv))),
pup_diff=c(-log10(c(p_up_LSD$p, unlist(p_up_tot)))+log10(c(p_up_LSD_inv$p, unlist(p_up_tot_inv)))), 
pdn_diff=c(-log10(c(p_dn_LSD$p,unlist(p_dn_tot)))+log10(c(p_dn_LSD_inv$p, unlist(p_dn_tot_inv)))), dataset=dataset,
type=type)
                   
library(openxlsx)
write.xlsx(df_all, file=file.path(OUT, "dementia_p_all.xlsx"))

###dementia-opposition scores (APR) per intervention (Table S13)
# APR = -log10(opposition p) - (-log10(concordance p)); infinite APRs are capped at +/-300
# relative APR = APR / mean APR of the aging-reverting comparators (positive controls)
cap_apr <- function(apr, p_opp, p_conc) ifelse(is.finite(apr), apr, 300 * sign(p_opp - p_conc))
apr_up  <- cap_apr(df_all$pup_diff, df_all$pup_rev, df_all$pup_mimick)
apr_dn  <- cap_apr(df_all$pdn_diff, df_all$pdn_rev, df_all$pdn_mimick)
apr_avg <- (apr_up + apr_dn)/2
is_pc   <- df_all$type == "Positive Control"
type_label <- c("Psychoplastogen"="Psychoplastogen", "Positive Control"="Aging-Reverting Comparator",
                "Negative Control"="Aging-Concordant Comparator")

apr_table <- data.frame(
  dataset = df_all$dataset,
  "opposition p (-log10) -up-regulated genes"    = df_all$pup_rev,
  "concordance p (-log10) -up-regulated genes"   = df_all$pup_mimick,
  "opposition p (-log10) -down-regulated genes"  = df_all$pdn_rev,
  "concordance p (-log10) -down-regulated genes" = df_all$pdn_mimick,
  "APR (-log10) -up-regulated genes"   = apr_up,
  "APR (-log10) -down-regulated genes" = apr_dn,
  "average APR" = apr_avg,
  type = unname(type_label[df_all$type]),
  "APR / mean Aging-Reverting Comparators APR - up-regulated genes"   = apr_up/mean(apr_up[is_pc]),
  "APR / mean Aging-Reverting Comparators APR - down-regulated genes" = apr_dn/mean(apr_dn[is_pc]),
  "average APR fold" = apr_avg/mean(apr_avg[is_pc]),
  check.names = FALSE
)
write.xlsx(apr_table, file=file.path(OUT, "dementia_opposition_scores.xlsx"))


df<-data.frame(pup=c(-log10(c(p_up_LSD$p, unlist(p_up_tot)))+log10(c(p_up_LSD_inv$p, unlist(p_up_tot_inv)))), 
               pdn=c(-log10(c(p_dn_LSD$p,unlist(p_dn_tot)))+log10(c(p_dn_LSD_inv$p, unlist(p_dn_tot_inv)))), dataset=dataset,
               type=type)

pdf(file.path(OUT, "Figures/All_psychedelicsAndCTRL_dementia.pdf"),6.5,5.5)
ggplot(df, aes(x=pup, y=pdn, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()



###create a table with all data

DE_annot <- data.table(
  DE        = names(p_up_list),
  DE_label  = dataset,
  DE_type   = type
)
library(data.table)

final_table <- rbindlist(
  lapply(names(p_up_list), function(DE) {
    
    dt_up <- data.table(
      DE        = DE,
      dataset_id = dat_names,
      condition  = condition,   
      direction = "UP",
      pvalue    = p_up_list[[DE]],
      cohen     = cohen_up_list[[DE]]
    )
    
    dt_dn <- data.table(
      DE        = DE,
      dataset_id = dat_names,
      condition  = condition,
      direction = "DN",
      pvalue    = p_dn_list[[DE]],
      cohen     = cohen_dn_list[[DE]]
    )
    
    rbind(dt_up, dt_dn)
  })
)

final_table <- merge(
  final_table,
  DE_annot,
  by = "DE",
  all.x = TRUE
)


library(data.table)
library(metap)

setDT(final_table)

p <- as.numeric(final_table$pvalue)   

inv_std <- fifelse(
  final_table$direction == "UP", final_table$cohen < 0,
  fifelse(final_table$direction == "DN", final_table$cohen > 0, NA)
)

inv_inv <- fifelse(
  final_table$direction == "UP", final_table$cohen > 0,
  fifelse(final_table$direction == "DN", final_table$cohen < 0, NA)
)

final_table[, `:=`(p_1t = NA_real_, p_1t_inv = NA_real_)]

ok_std <- !is.na(p) & !is.na(inv_std)
ok_inv <- !is.na(p) & !is.na(inv_inv)

final_table[, `:=`(
  inverted_1t     = as.logical(inv_std),  # TRUE/FALSE/NA
  inverted_1t_inv = as.logical(inv_inv)
)]

final_table$p_1t[ok_std] <- metap::two2one(
  p[ok_std],
  two    = rep(TRUE, sum(ok_std)),
  invert = inv_std[ok_std]
)

final_table$p_1t_inv[ok_inv] <- metap::two2one(
  p[ok_inv],
  two    = rep(TRUE, sum(ok_inv)),
  invert = inv_inv[ok_inv]
)

final_table$DE<-gsub("DE_", "", final_table$DE)

write.xlsx(final_table, file=file.path(OUT, "all_GSEA_dementia.xlsx"))


###table of all dementia tests (Table S12), input of 21_cohen_figures.R
# Cohen's d above is computed as (Healthy - Dementia); here it is reported as
# (Dementia - Healthy), so that d > 0 means higher gene-set activity in dementia.
# p-values are unchanged: inverted_1t / p_1t are the dementia-reversal tests,
# inverted_1t_inv / p_1t_inv the dementia-mimicking tests.

dementia_tests <- data.frame(
  intervention_dataset   = final_table$DE,
  intervention_condition = final_table$DE_label,
  intervention_type      = final_table$DE_type,
  dementia_dataset       = final_table$dataset_id,
  dementia_condition     = final_table$condition,
  direction              = final_table$direction,
  pvalue                 = final_table$pvalue,
  cohen                  = -final_table$cohen,
  invert_rev             = final_table$inverted_1t,
  invert_mimick          = final_table$inverted_1t_inv,
  p_onetail_rev          = final_table$p_1t,
  p_onetail_mimick       = final_table$p_1t_inv
)

dementia_tests_desc <- data.frame(
  Column = colnames(dementia_tests),
  Description = c(
    "GEO ID (where available) of the intervention dataset used to define differentially expressed genes (DEGs)",
    "description of the intervention",
    "type of intervention: Psychoplastogen, Positive Control, Negative Control",
    "GEO ID of the dementia dataset used to test differences in intervention DEGs signatures",
    "Condition and brain area analysed",
    "direction of differentially expressed genes. UP= up with intervention; DN= down with intervention",
    "t-test p-value",
    "Cohen's d (dementia vs healthy) of the ssGSEA scores of the intervention DEGs",
    "whether the p-value was inverted or not to compute the one-tailed dementia reversal p-value",
    "whether the p-value was inverted or not to compute the one-tailed dementia mimicking p-value",
    "one-tail p-value for dementia reversal",
    "one-tail p-value for dementia mimicking"
  )
)

wb_dem <- createWorkbook()
addWorksheet(wb_dem, "Column description")
writeData(wb_dem, "Column description", dementia_tests_desc)
addWorksheet(wb_dem, "Full test results")
writeData(wb_dem, "Full test results", dementia_tests)
saveWorkbook(wb_dem, file.path(OUT, "dementia_cohen_all_tests.xlsx"), overwrite = TRUE)

