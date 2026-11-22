suppressMessages(library(DESeq2))  # DESeqResults objects need DESeq2/S4Vectors loaded (rownames() returns NULL otherwise)
rm(list=ls())
setwd("workdir")#workdir = working directory
set.seed(1234)  # fixed seed; where a random state is set after load() below, that state is used instead (README.md, "Random numbers")

load("Results/RData/DE_GSE179379.RData")
load("Data/original_RData_notclosed/Alldata_20Sep.RData")
assign(".Random.seed", readRDS("Data/original_RData_notclosed/rng_states.rds")$Alldata_20Sep, envir = .GlobalEnv)  # random state of the published runs (the .Random.seed saved in Data/original_RData_notclosed/Alldata_20Sep.RData), set explicitly
homologs<-read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022


library(GEOquery)
library(ggplot2)
#preprocessing
#change names in gene symbols

changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]
  
  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }
  
  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}


################################
##########MSBB
###############################
# MSBB (AD Knowledge Portal) is controlled-access data. Without Data/raw_data_closed/AD portal/, the
# gene-wise correlations with age and with the Braak score are read from
# Data/closed_data_intermediate/MSBB_gene_correlations.RData (written by this script
# when the MSBB files are present). The GSEA below is the same in both cases.
has_msbb <- file.exists("Data/raw_data_closed/AD portal/MSBB_Normalized_counts_(CQN).tsv")
if (!has_msbb) {
  load("Data/closed_data_intermediate/MSBB_gene_correlations.RData")  # cor_age, cor_braak
  message("Data/raw_data_closed/AD portal not found: MSBB gene correlations read from Data/closed_data_intermediate")
}
if (has_msbb) {
MSBB<-read.csv("Data/raw_data_closed/AD portal/MSBB_Normalized_counts_(CQN).tsv", sep="\t", row.names = 1)
MSBB_meta<-read.csv("Data/raw_data_closed/AD portal/RNAseq_Harmonization_MSBB_combined_metadata.csv")

#reorder metadata to mach order of counts
MSBB_meta<-MSBB_meta[match(colnames(MSBB), MSBB_meta[,1]),]

MSBB_sel<-MSBB[,which(MSBB_meta$tissue=="frontal pole")]
age<-MSBB_meta$ageDeath[which(MSBB_meta$tissue=="frontal pole")]
age[age=="90+"]<-91
age<-as.numeric(age)

# Ensembl 115 (the live release in Sep 2025) instead of live biomaRt
source("Scripts/data_preparation/ensembl_release_map.R")
ens115 <- ensembl_release_map(115, "hsapiens", "Data/ensembl")
ids <- ens115[ens115$ensembl_gene_id %in% rownames(MSBB_sel), ]


MSBB_sel<-changenames(data=MSBB_sel, anno=ids)

pca<-prcomp(t(MSBB_sel))
df<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2], MSBB_meta[which(MSBB_meta$tissue=="frontal pole"),])
ggplot(df, aes(PC1, PC2, colour=age))+geom_point()
cor_age<-cor(t(MSBB_sel), as.numeric(age))
}


library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])


genes_cor<-cor_age
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

forgesea<-unlist(genes_cor)
names(forgesea)<-rownames(genes_cor)
forgesea<-forgesea[!is.na(forgesea)]
forgesea<-forgesea[names(forgesea) %in% rat_homologs]
fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)

p1<-plotEnrichment(genes_up, forgesea)
p2<-plotEnrichment(genes_dn, forgesea)

df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data))))

p<-ggplot(df, aes(rank, ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf("Results/Figures/GSEA_MSBB.pdf", 4, 4)
print(p)
dev.off()


########
#cognitive performance

if (has_msbb) {
MSBB_sel<-MSBB[,which(MSBB_meta$tissue=="frontal pole" & MSBB_meta$ageDeath=="90+")]
Braak<-MSBB_meta$Braak[which(MSBB_meta$tissue=="frontal pole" & MSBB_meta$ageDeath=="90+")]

ids <- ens115[ens115$ensembl_gene_id %in% rownames(MSBB_sel), ]  # Ensembl 115


MSBB_sel<-changenames(data=MSBB_sel, anno=ids)
cor_braak<-cor(t(MSBB_sel), as.numeric(Braak))
save(cor_age, cor_braak, file="Results/RData/MSBB_gene_correlations.RData")  # gene-level summary, no sample data
}


library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])


genes_cor<-cor_braak
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

forgesea<-unlist(genes_cor)
names(forgesea)<-rownames(genes_cor)
forgesea<-forgesea[!is.na(forgesea)]
forgesea<-forgesea[names(forgesea) %in% rat_homologs]
fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)

p1<-plotEnrichment(genes_up, forgesea)
p2<-plotEnrichment(genes_dn, forgesea)

df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data))))

p<-ggplot(df, aes(rank, ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf("Results/Figures/GSEA_MSBB_Braak.pdf", 4, 4)
print(p)
dev.off()




