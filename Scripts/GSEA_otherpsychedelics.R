################################
####### gsea lsd vs aging
###############################
library(fgsea)
library(ggplot2)
load("Results/RData/ExerciseAndAlcohol.RData")
load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Psychedelics metanalysis/Psychedelics_PFC.RData")
load("RData/Alldata_7Jul.RData")
rat_homologs<-read.csv("Data/Human rat homologs.txt")
mouse_homologs<-read.csv("Data/Human mouse homologs.txt", sep="\t")

mouse_data<-c("DE_GSE64607","DE_GSE38465_SAMP8","DE_GSE38465_SAMR1","DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57", "DE_MDMA", "DE_GSE161626","DE_GSE161626_48h","DE_GSE161626_7d","DE_GSE81672","DE_GSE26364","DE_GSE209859",
              "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks", "DE_GSE30880")
rat_data<-c("DE_GSE14720", "DE_GSE23728","DE_GSE179380","DE_DMT", "DE_pharm", "DE_harm")
##MDMA
DEGs_array<-c("DE_MDMA",
              
              "DE_GSE26364",
              #DOI
              
              "DE_GSE23728",
              #LSD
              "DE_GSE179380", 
              "DE_GSE64607","DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57",
              "DE_GSE30880")

DEGs_seq<-c("DE_DMT", "DE_GSE161626","DE_GSE161626_48h","DE_GSE161626_7d", "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
"DE_GSE81672", "DE_harm", "DE_pharm")

library(metap)

meta_dn_rev<-function(x){
  istwo <- rep(T, 14)
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
  istwo <- rep(T, 14)
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
  istwo <- rep(T, 14)
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
  istwo <- rep(T, 14)
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
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
dat<-setdiff(dat, c("GSE5388"))#outlier

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
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]
    
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
  
  p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
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
print(ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
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
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  dat<-setdiff(dat, c("GSE5388"))#outlier
  
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
    
    p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
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
  print(ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
  save(fgsea_res, file=paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
}

save(pup, file="Results/pup_otherpsychedelics.RData")
save(pdn, file="Results/pdn_otherpsychedelics.RData")

plot(-log10(unlist(pup)), -log10(unlist(pdn)))
df<-data.frame(pup=-log10(unlist(pup)), pdn=-log10(unlist(pdn)), dataset=c(DEGs_array, DEGs_seq))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label()


##################
## collapsing GSEA p (no random drugs)
##################
p_dn<-c()
p_up<-c()
for(DE in c(DEGs_array, DEGs_seq)){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn<-c(p_dn, meta_dn(fgsea_res))
  p_up<-c(p_up, meta_up(fgsea_res))
  
}
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn<-c(p_dn, meta_dn(fgsea_res))
p_up<-c(p_up, meta_up(fgsea_res))

#p_dn[p_dn==0]<-10^(-200)

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_array, DEGs_seq, "LSD"))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")


p_dn_rev<-c()
p_up_rev<-c()
for(DE in c(DEGs_array, DEGs_seq)){
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
  
  p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
  p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))
  
}
load(paste("Results/RData/fgsea_res.RData", sep=""))
p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))

df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c(DEGs_array, DEGs_seq, "LSD"))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_array, DEGs_seq, "LSD")))

dataset<-c("MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise", "Alcohol", "Alcohol", "Alcohol",
           "Enriched Environment", "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "LSD (chronic)")
type<-rep("Psychoplastogen", length(dataset))
type[c(5,9)]<-"Positive Control"
type[c(6,7,8)]<-"Negative Control"

df<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=dataset, direction="increase aging",
               type=type)

pdf("Results/Figures/All_psychedelicsAndCTRL.pdf",6.5,5.5)
ggplot(df, aes(x=pup_diff, y=pdn_diff, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()


#######################
### random drugs
######################
##fatto senza alcuni dataset, nel caso aggiungere DOI 48h e 7day
DEGs_array<-c("DE_MDMA",
              "DE_GSE47541",
              "DE_GSE26364",
              #DOI
              "DE_GSE14720",
              "DE_GSE23728",
              #LSD
              "DE_GSE179380")

DEGs_seq<-c("DE_DMT", "DE_GSE161626","DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
            "DE_GSE81672", "DE_harm", "DE_pharm")


load("Results/RData/CMAP_NPC_mean.RData")
load("Results/RData/compounds_NPC.RData")

for(DE in DEGs_array){
fgsea_res_drugs<-list()

if(DE %in% mouse_data){
  homologs<-mouse_homologs
} else if(DE %in% rat_data){
  homologs<-rat_homologs
} 

DEGs<-get(DE)
genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

for(comp in 1:ncol(data_mean_all)){
  print(paste(DE,comp))
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE5388"))
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  genes_up_rand<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=T)[1:length(genes_up)]]
  genes_dn_rand<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=F)[1:length(genes_dn)]]
  
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
    fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up_rand, DN=genes_dn_rand), forgesea)
    
  }
  
}
save(fgsea_res_drugs, file=paste("Results/RData/fgsea_res_drugs_",DE, ".RData", sep=""))
}



for(DE in DEGs_seq){
  fgsea_res_drugs<-list()
  
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } 
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  for(comp in 1:ncol(data_mean_all)){
    print(paste(DE,comp))
    region<-c("DLPFC")
    diagnosis<-"Healthy"
    dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
    dat<-setdiff(dat, c("GSE5388"))
    dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
    
    genes_up_rand<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=T)[1:length(genes_up)]]
    genes_dn_rand<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=F)[1:length(genes_dn)]]
    
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
      fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up_rand, DN=genes_dn_rand), forgesea)
      
    }
    
  }
  save(fgsea_res_drugs, file=paste("Results/RData/fgsea_res_drugs_",DE, ".RData", sep=""))
}


################################
###evaluation with nes
##############################
library(metap)

p_dn<-c()
p_up<-c()
for(DE in c(DEGs_array, DEGs_seq)){
  load(paste("Results/RData/fgsea_res_drugs_",DE, ".RData", sep=""))
  load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
nes_dn<-unlist(lapply(fgsea_res, function(x){x[1,6]}))
nes_dn_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,6]}))})

better_dn<-lapply(nes_dn_drug, function(x){x>nes_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset

nes_up<-unlist(lapply(fgsea_res, function(x){x[2,6]}))
nes_up_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,6]}))})

better_up<-lapply(nes_up_drug, function(x){x<nes_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
p_dn<-c(p_dn, sumlog(colSums(better_dn_mat, na.rm=T)/193)[[3]])
p_up<-c(p_up, sumlog(colSums(better_up_mat, na.rm=T)/193)[[3]])

}

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_array, DEGs_seq))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label()

#LSD RNA-seq
load("Results/RData/fgsea_res.RData")
load("Results/RData/fgsea_res_drugs_RAW.RData")
nes_dn<-unlist(lapply(fgsea_res, function(x){x[1,6]}))
nes_dn_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,6]}))})

better_dn<-lapply(nes_dn_drug, function(x){x>nes_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset

nes_up<-unlist(lapply(fgsea_res, function(x){x[2,6]}))
nes_up_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,6]}))})

better_up<-lapply(nes_up_drug, function(x){x<nes_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
p_dn<-c(p_dn, sumlog(colSums(better_dn_mat, na.rm=T)/193)[[3]])
p_up<-c(p_up, sumlog(colSums(better_up_mat, na.rm=T)/193)[[3]])

library(ggrepel)
df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_array, DEGs_seq, "LSD"))
pdf("Results/Figures/all_psychedelics_rand_drugs.pdf",5,5)
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")
dev.off()

######evaluation collapsing the pvalue

library(metap)

meta_dn<-function(x){
  istwo <- rep(T, 14)
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
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}


#############

####TODO
load(paste("Results/RData/fgsea_res_", DE, ".RData", sep=""))
load(paste("Results/RData/fgsea_res_drugs_",DE, ".RData", sep=""))

nes_dn<-unlist(lapply(fgsea_res, function(x){x[1,6]}))
nes_dn_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,6]}))})

df<-data.frame(NES=unlist(nes_dn_drug), dataset=rep(dat_names, each=length(nes_dn_drug)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(NES=nes_dn, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_nes_dn_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=NES))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=NES), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_bw()
dev.off()

better_dn<-lapply(nes_dn_drug, function(x){x>nes_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
colSums(better_dn_mat)/193

pdf("Results/Figures/Random_drugs_nes_dn_RAW_scatter.pdf", 4,4)
plot(-log10(pvals_dn), 1-colSums(better_dn_mat)/193, pch=19, ylab="% of tests with NES higher than LSD (dn)")
dev.off()

nes_up<-unlist(lapply(fgsea_res, function(x){x[2,6]}))
nes_up_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,6]}))})

df<-data.frame(NES=unlist(nes_up_drug), dataset=rep(dat_names, each=length(nes_up_drug)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(NES=nes_up, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_nes_up_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=NES))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=NES), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_bw()
dev.off()

####
better_up<-lapply(nes_up_drug, function(x){x<nes_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat)/193
pdf("Results/Figures/Random_drugs_nes_up_RAW_scatter.pdf", 4,4)
plot(-log10(pvals_up), 1-colSums(better_up_mat)/193, pch=19, ylab="% of tests with NES lower than LSD (up)")
dev.off()

pdf("Results/Figures/Random_drugs_nes_perc_RAW.pdf", 4,4)
plot(1-colSums(better_dn_mat, na.rm=T)/193, 1-colSums(better_up_mat, na.rm=T)/193, xlab="% of tests with higher NES than LSD (dn)",
     ylab="% of tests with NES lower than LSD (up)", pch=19)
dev.off()

plot(colSums(better_dn_mat, na.rm=T)/193, colSums(better_up_mat, na.rm=T)/193, xlab="% of tests with higher NES than LSD (dn)",
     ylab="% of tests with NES lower than LSD (up)", pch=19)

p_dn<-sumlog(colSums(better_dn_mat, na.rm=T)/193)
p_up<-sumlog(colSums(better_up_mat, na.rm=T)/193)

}



