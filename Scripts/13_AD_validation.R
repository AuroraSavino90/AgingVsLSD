rm(list=ls())
setwd("workdir")#workdir = working directory

load("Results/RData/DE_GSE179379.RData")
load("Data/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")


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
MSBB<-read.csv("Data/AD portal/MSBB_Normalized_counts_(CQN).tsv", sep="\t", row.names = 1)
MSBB_meta<-read.csv("Data/AD portal/RNAseq_Harmonization_MSBB_combined_metadata.csv")

#reorder metadata to mach order of counts
MSBB_meta<-MSBB_meta[match(colnames(MSBB), MSBB_meta[,1]),]

MSBB_sel<-MSBB[,which(MSBB_meta$tissue=="frontal pole")]
age<-MSBB_meta$ageDeath[which(MSBB_meta$tissue=="frontal pole")]
age[age=="90+"]<-91
age<-as.numeric(age)

library(biomaRt)
ensembl.hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(MSBB_sel),
           mart = ensembl.hs)


MSBB_sel<-changenames(data=MSBB_sel, anno=ids)

pca<-prcomp(t(MSBB_sel))
df<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2], MSBB_meta[which(MSBB_meta$tissue=="frontal pole"),])
ggplot(df, aes(PC1, PC2, colour=age))+geom_point()


library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])


genes_cor<-cor(t(MSBB_sel), as.numeric(age))
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

MSBB_sel<-MSBB[,which(MSBB_meta$tissue=="frontal pole" & MSBB_meta$ageDeath=="90+")]
Braak<-MSBB_meta$Braak[which(MSBB_meta$tissue=="frontal pole" & MSBB_meta$ageDeath=="90+")]

library(biomaRt)
ensembl.hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(MSBB_sel),
           mart = ensembl.hs)


MSBB_sel<-changenames(data=MSBB_sel, anno=ids)


library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])


genes_cor<-cor(t(MSBB_sel), as.numeric(Braak))
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




