rm(list=ls())
setwd("workdir")#workdir = working directory

library(GEOquery)
library(ggplot2)
source("Scripts/utilities.R")
options(timeout = 5000)

################################
##########GSE36192
###############################

gds<-getGEO("GSE36192",destdir="Data/", AnnotGPL = TRUE)

GSE36192<-exprs(gds$GSE36192_series_matrix.txt.gz)
GSE36192_meta<-pData(gds$GSE36192_series_matrix.txt.gz)
GSE36192_anno<-fData(gds$GSE36192_series_matrix.txt.gz)

GSE36192<-changenames(data=GSE36192, anno=cbind(GSE36192_anno$ID, GSE36192_anno$`Gene symbol`))
GSE36192<-log2(GSE36192)

################################
##########GSE25219
###############################

gds<-getGEO("GSE25219",destdir="Data/", AnnotGPL = TRUE)

GSE25219<-exprs(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)
GSE25219_meta<-pData(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)
GSE25219_anno<-fData(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)

symbols<-gsub('([^/]*)//([^/]*)//.*', '\\2', GSE25219_anno$gene_assignment)
symbols<-gsub(' ', '', symbols)

GSE25219<-changenames(data=GSE25219, anno=cbind(GSE25219_anno$ID, symbols))

##########


library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])



#############

  sample<-GSE36192_meta$geo_accession[which(GSE36192_meta$`tissue:ch1`=="frontal cortex")]
  data<-GSE36192[, sample]
  age<-GSE36192_meta$`age (y):ch1`[which(GSE36192_meta$`tissue:ch1`=="frontal cortex")]
  age<-as.numeric(age)*365
  sample<-sample[which(age>=7300)]
  data<-data[, sample]
  age<-age[which(age>=7300)]

  genes_cor<-cor(t(data), as.numeric(age))
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

  pdf("Results/Figures/GSEA_GSE36192.pdf", 4, 4)
 print(p)
 dev.off()

  #############

  sample<-GSE25219_meta$geo_accession[which(GSE25219_meta$`region:ch1`=="DFC")]
  data<-GSE25219[, sample]
  age<-GSE25219_meta$`age:ch1`[which(GSE25219_meta$`region:ch1`=="DFC")]

  sample<-sample[grep(" Y", age)]
  data<-data[, sample]
  age<-age[grep(" Y", age)]

  age<-gsub(" Y", "", age)
  age<-as.numeric(age)*365
  sample<-sample[which(age>=7300)]
  data<-data[, sample]
  age<-age[which(age>=7300)]

  genes_cor<-cor(t(data), as.numeric(age))
  rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

  forgesea<-unlist(genes_cor)
  names(forgesea)<-rownames(genes_cor)
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea, nPermSimple=100000)

  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)

  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data))))

  p<-ggplot(df, aes(rank, ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

  pdf("Results/Figures/GSEA_GSE25219.pdf", 4, 4)
  print(p)
  dev.off()

  
