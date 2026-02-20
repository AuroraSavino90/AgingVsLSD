rm(list=ls())
setwd("workdir")#workdir = working directory
################################
####### gsea lsd vs aging
###############################
library(fgsea)
library(ggplot2)
library(limma)
library(ggrepel)
load("Data/ExerciseAndAlcohol.RData")
load("Data/Psychedelics_PFC.RData")
load("Data/Alldata_20Sep.RData")
source("Scripts/utilities.R")
rat_homologs<-read.csv("Data/Human rat homologs.txt")
mouse_homologs<-read.csv("Data/Human mouse homologs.txt", sep="\t")

mouse_data<-c("DE_GSE64607","DE_GSE164798_chronic", "DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57", "DE_MDMA", "DE_GSE161626",
              "DE_GSE161626_48h","DE_GSE161626_7d","DE_GSE81672","DE_GSE26364","DE_GSE209859",
              "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
               "DE_GSE111273_Y", "DE_GSE111273_O")
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
            "DE_GSE81672", "DE_harm", "DE_pharm", "DE_GSE164798_chronic",  "DE_GSE111273_Y", "DE_GSE111273_O")
library(metap)

meta_dn_rev<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up_rev<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

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


pup<-list()
pdn<-list()
cor_up_tot<-list()
cor_dn_tot<-list()
n<-0
i<-0
for(DE in DEGs_array){
  i<-i+1
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
  n<-n+1
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance

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
dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance

for(dd in dat){
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
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
fgsea_res<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
  
  p[[n]]<-ggplot(df[[n]], aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
  
}

pdn[[i]]<-meta_dn(fgsea_res)
pup[[i]]<-meta_up(fgsea_res)


df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}


fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf(paste("Results/Figures/GSEA_all_datasets_",DE,".pdf", sep=""), 8, 8)
print(ggplot(df_tot, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
dev.off()
save(fgsea_res, file=paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))

}

######################### RNA-seq

cor_up_tot<-list()
cor_dn_tot<-list()
n<-0

for(DE in DEGs_seq){
  i<-i+1
  n<-n+1
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
  
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance
  
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
  dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance
  
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
  fgsea_res<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){
    
    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
    
    p1<-plotEnrichment(genes_up, forgesea)
    p2<-plotEnrichment(genes_dn, forgesea)
    
    df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
    
    p[[n]]<-ggplot(df[[n]], aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
      scale_colour_manual(values=c("dn"="blue", "up"="red"))
    
  }
  
  pdn[[i]]<-meta_dn_rev(fgsea_res)
  pup[[i]]<-meta_up_rev(fgsea_res)
  
  df_tot<-df[[1]]
  for(n in 2:length(genes_cor)){
    df_tot<-rbind.data.frame(df_tot, df[[n]])
  }
  
  
  fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})
  
  df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
  
  pdf(paste("Results/Figures/GSEA_all_datasets_",DE,".pdf", sep=""), 8, 8)
  print(ggplot(df_tot, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
  save(fgsea_res, file=paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
}
save(dat_names, file="Results/RData/dat_names.RData")

##################
## collapsing GSEA p (no random drugs)
##################
#load gsea tests for all interventions and merge p as age-mimicking
p_dn<-c()
p_up<-c()
for(DE in c(DEGs_array, DEGs_seq)){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn<-c(p_dn, meta_dn(fgsea_res))
  p_up<-c(p_up, meta_up(fgsea_res))
  
}
#load gsea for chronic lsd and merge p as age-mimicking
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn<-c(p_dn, meta_dn(fgsea_res))
p_up<-c(p_up, meta_up(fgsea_res))

#p_dn[p_dn==0]<-10^(-200)

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_array, DEGs_seq, "LSD"))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")

#load gsea tests for all interventions and merge p for age-reversal
p_dn_rev<-c()
p_up_rev<-c()
for(DE in c(DEGs_array, DEGs_seq)){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
  p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))
  
}
#load gsea for chronic lsd and merge p for age-reversal
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))


df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c(DEGs_array, DEGs_seq, "LSD"))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_array, DEGs_seq, "LSD")))

dataset<-c("MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise", "Alcohol", "Alcohol", "Alcohol",
           "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "Exercise","EE (Young)", "EE (Old)", 
           "LSD (chronic)")
type<-rep("Psychoplastogen", length(dataset))
type[c(5, 18:20)]<-"Positive Control"
type[c(6,7,8)]<-"Negative Control"

df<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=dataset,
               type=type)

pdf("Results/Figures/All_psychedelicsAndCTRL.pdf",6.5,5.5)
ggplot(df, aes(x=pup_diff, y=pdn_diff, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()


###################
### adding rats aging as positive control
###################
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

#####genes correlation with age in humans across all datasets (already calculated above)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance

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
dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance

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

####define DEGs in rats aging
GSE75772<-read.csv("Data/GSE75772_DEseq-normalized.counts_mPFC.csv.gz", row.names = 1)

library(biomaRt)
ensembl.rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'rgd_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(GSE75772),
           mart = ensembl.rat)


GSE75772<-changenames(data=GSE75772, anno=ids)
GSE75772<-log2(GSE75772+1)

GSE75772_meta<-read.csv("Data/GSE75772_series_matrix.txt", sep="\t")

homologs<-rat_homologs

data<-get("GSE75772")
age<-unlist(GSE75772_meta[10,c(2:32)])

cond2<-data[,age=="age: Young"]#ctrl
cond1<-data[,age=="age: Aged"]#passive administration

DE<-DExpr(cond1, cond2)

genes_up<-rownames(DE)[which(DE$logFC>0 & DE$P.Value<0.05)]
genes_dn<-rownames(DE)[which(DE$logFC<0 & DE$P.Value<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_res<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
  
  p[[n]]<-ggplot(df[[n]], aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
  
}

pdn[[i]]<-meta_dn(fgsea_res)
pup[[i]]<-meta_up(fgsea_res)

df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}

fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/GSEA_all_datasets_GSE75772.pdf", 8, 8)
print(ggplot(df_tot, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
        scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
dev.off()
save(fgsea_res, file="Results/RData/fgsea_res_GSE75772.RData")

##################
## collapsing GSEA p (no random drugs)
##################
p_dn<-c()
p_up<-c()
for(DE in c(DEGs_array, DEGs_seq, "GSE75772")){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn<-c(p_dn, meta_dn(fgsea_res))
  p_up<-c(p_up, meta_up(fgsea_res))
  
}
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn<-c(p_dn, meta_dn(fgsea_res))
p_up<-c(p_up, meta_up(fgsea_res))

#p_dn[p_dn==0]<-10^(-200)

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_array, DEGs_seq, "GSE75772", "LSD"))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")


p_dn_rev<-c()
p_up_rev<-c()
for(DE in c(DEGs_array, DEGs_seq, "GSE75772")){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
  p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))
  
}
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))

df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c(DEGs_array, DEGs_seq, "GSE75772", "LSD"))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_array, DEGs_seq, "GSE75772", "LSD")))

dataset<-c("MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise", "Alcohol", "Alcohol", "Alcohol",
           "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "Exercise","EE (Young)", "EE (Old)", 
            "Rat aging", "LSD (chronic)")
type<-rep("Psychoplastogen", length(dataset))
type[c(5, 18:20)]<-"Positive Control"
type[c(6,7,8,21)]<-"Negative Control"

df_all<-data.frame(pup_rev=df_rev$pup, pdn_rev=df_rev$pdn, pup=df$pup, pdn=df$pdn, pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=dataset,
                   type=type)
df_all[df_all=="-Inf"]<-(-300)
df_all[df_all< (-300)]<-(-300)
df_all[df_all=="Inf"]<-(300)
df_all$pup_diff= -df$pup+df_rev$pup 
df_all$pdn_diff= -df$pdn+df_rev$pdn
library(openxlsx)
write.xlsx(df_all, file="Results/age_reversal_p_all.xlsx")

df<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=dataset,
               type=type)
df[df=="-Inf"]<-(-300)
df[df< (-300)]<-(-300)
df[df=="Inf"]<-(300)

pdf("Results/Figures/All_psychedelicsAndCTRL_rataging.pdf",6.5,5.5)
ggplot(df, aes(x=pup_diff, y=pdn_diff, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()


##save all tests' statistics in a table
library(data.table)

all_res <- list()
ind<-0
for (DE in c(DEGs_array, DEGs_seq , "GSE75772")) {
  ind<-ind+1
  
  load(paste0("Results/RData/fgsea_res_", DE, ".RData"))
  
  dt_DE <- rbindlist(
    lapply(seq_along(fgsea_res), function(i) {
      x <- as.data.table(fgsea_res[[i]])
      x[, `:=`(aging_dataset=dat_names[i], DE = DE, dataset = dataset[ind], type = type[ind])]
      x
    }),
    use.names = TRUE, fill = TRUE
  )
   setcolorder(dt_DE, c( "aging_dataset", "DE", "dataset", setdiff(names(dt_DE), c("DE", "dataset",  "aging_dataset"))))
  
  all_res[[DE]] <- dt_DE
}

ind<-ind+1
load(paste("Results/RData/fgsea_res.RData", sep=""))
dt_last <- rbindlist(
  lapply(seq_along(fgsea_res), function(i) {
    
    x <- as.data.table(fgsea_res[[i]])
    
    x[, `:=`(
      aging_dataset=dat_names[i],
      DE      = "DE_GSE179379",
      dataset = dataset[ind],
      type    = "Psychoplastogen"
    )]
    
    x
  }),
  use.names = TRUE,
  fill = TRUE
)
setcolorder(dt_last, c( "aging_dataset", "DE", "dataset", setdiff(names(dt_DE), c("DE", "dataset",  "aging_dataset"))))
all_res[["DE_GSE179379"]] <- dt_last

# tabella finale con TUTTI i DE
gsea_table <- rbindlist(all_res, use.names = TRUE, fill = TRUE)
gsea_table<-gsea_table[,c("aging_dataset", "DE",            "dataset",       "pathway",          "padj",         
                          "NES" ,          "size" ,         "leadingEdge" ,  "type")]


colnames(gsea_table)[2:4]<-c("dataset", "intervention", "direction")
library(openxlsx)
gsea_table$dataset<-gsub("DE_", "", gsea_table$dataset)

setDT(gsea_table)

# Definizione invert come prima
gsea_table[, invert_mimick := fifelse(
  direction == "DN", sign(NES) ==  1,
  fifelse(direction == "UP", sign(NES) == -1, NA)
)]

gsea_table[, invert_rev := fifelse(
  direction == "DN", sign(NES) == -1,
  fifelse(direction == "UP", sign(NES) ==  1, NA)
)]

# Creo colonne vuote
gsea_table[, p_two2one_meta := NA_real_]
gsea_table[, p_two2one_rev  := NA_real_]

# Indici validi (p dentro (0,1] e invert non NA)
valid_mimick <- which(
  !is.na(gsea_table$padj) &
    gsea_table$padj > 0 &
    gsea_table$padj <= 1 &
    !is.na(gsea_table$invert_mimick)
)

valid_rev <- which(
  !is.na(gsea_table$padj) &
    gsea_table$padj > 0 &
    gsea_table$padj <= 1 &
    !is.na(gsea_table$invert_rev)
)

# Applico two2one solo ai validi
gsea_table$p_two2one_meta[valid_mimick] <-
  metap::two2one(
    gsea_table$padj[valid_mimick],
    two = rep(TRUE, length(valid_mimick)),
    invert = gsea_table$invert_mimick[valid_mimick]
  )

gsea_table$p_two2one_rev[valid_rev] <-
  metap::two2one(
    gsea_table$padj[valid_rev],
    two = rep(TRUE, length(valid_rev)),
    invert = gsea_table$invert_rev[valid_rev]
  )

write.xlsx(gsea_table, file="Results/all_GSEA_aging.xlsx")



