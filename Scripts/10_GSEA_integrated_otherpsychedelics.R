rm(list=ls())
setwd("workdir")#workdir = working directory
#################################################
#########GSEA with the batch corrected dataset
##################################################
load("Data/ExerciseAndAlcohol.RData")
load("Results/RData/alldata_resid.RData")
load("Results/RData/age_all.RData")

load("Data/Psychedelics_PFC.RData")
load("Data/Alldata_20Sep.RData")
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

library(fgsea)
library(ggplot2)
fgsea_tot<-list()
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
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
genes_cor<-cor(t(alldata_resid), as.numeric(age_all))

universe<-unique(homologs[homologs[,2] %in% rownames(DEGs),3])

forgesea<-unlist(genes_cor)
names(forgesea)<-rownames(genes_cor)
forgesea<-forgesea[!is.na(forgesea)]
forgesea<-forgesea[names(forgesea) %in% universe]
fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
fgsea_tot[[i]]<-fgsea_res

print(DE)
print(fgsea_res)

p1<-plotEnrichment(genes_up, forgesea)
p2<-plotEnrichment(genes_dn, forgesea)

df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), 
                     data.frame(p2$data, dir=rep("dn", nrow(p2$data))))

p<-ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf(paste("Results/Figures/GSEA_all_batchcorrect_", DE, ".pdf", sep=""), 5,4)
print(ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red")))
dev.off()
}


for(DE in DEGs_seq){
  i<-i+1
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
  
  genes_cor<-cor(t(alldata_resid), as.numeric(age_all))
  
  universe<-unique(homologs[homologs[,2] %in% rownames(DEGs),3])
  
  forgesea<-unlist(genes_cor)
  names(forgesea)<-rownames(genes_cor)
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% universe]
  fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
  fgsea_tot[[i]]<-fgsea_res
  
  print(DE)
  print(fgsea_res)
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), 
                       data.frame(p2$data, dir=rep("dn", nrow(p2$data))))
  
  p<-ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
  
  pdf(paste("Results/Figures/GSEA_all_batchcorrect_", DE, ".pdf", sep=""), 5,4)
  print(ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red")))
  dev.off()
}

i<-i+1
load("Results/RData/DE_GSE179379.RData")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(rat_homologs[which(rat_homologs[,2] %in% genes_up),3])
genes_dn<-unique(rat_homologs[which(rat_homologs[,2] %in% genes_dn),3])

genes_cor<-cor(t(alldata_resid), as.numeric(age_all))

universe<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

forgesea<-unlist(genes_cor)
names(forgesea)<-rownames(genes_cor)
forgesea<-forgesea[!is.na(forgesea)]
forgesea<-forgesea[names(forgesea) %in% universe]
fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
fgsea_tot[[i]]<-fgsea_res

p_up<-c()
NES_up<-c()
p_dn<-c()
NES_dn<-c()
for(i in 1:length(fgsea_tot)){
p_up<-c(p_up, fgsea_tot[[i]]$pval[2])
NES_up<-c(NES_up, fgsea_tot[[i]]$NES[2])
p_dn<-c(p_dn, fgsea_tot[[i]]$pval[1])
NES_dn<-c(NES_dn, fgsea_tot[[i]]$NES[1])
}

library(ggrepel)
dataset<-c("MDMA", "Ketamine (long term)", "DOI (2h)", 
           "LSD (single)", "Exercise", "Alcohol", "Alcohol", "Alcohol",
           "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
           "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
           "Harmaline", "Pharmahuasca", "Exercise","EE (Young)", "EE (Old)", 
           "LSD (chronic)")
type<-rep("Psychoplastogen", length(dataset))
type[c(5, 18:20)]<-"Positive Control"
type[c(6,7,8)]<-"Negative Control"

df<-data.frame(NES_up=NES_up, NES_dn=NES_dn, dataset=dataset,
               type=type)

pdf("Results/Figures/All_psychedelics_integratedAndCTRL.pdf",6.5,5.5)
ggplot(df, aes(x=NES_up, y=NES_dn, label=dataset, colour=type))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()

