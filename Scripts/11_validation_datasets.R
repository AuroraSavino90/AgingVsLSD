rm(list=ls())
setwd("workdir")#workdir = working directory

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
##########GSE46706, gene-level
###############################GSE60862

gds<-getGEO("GSE60862",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging", AnnotGPL = TRUE)

GSE60862<-exprs(gds$GSE60862_series_matrix.txt.gz)
GSE60862_meta<-pData(gds$GSE60862_series_matrix.txt.gz)
GSE60862_anno<-fData(gds$GSE60862_series_matrix.txt.gz)

symbols<-gsub('([^/]*)//([^/]*)//.*', '\\2', GSE60862_anno$gene_assignment)
symbols<-gsub(' ', '', symbols)

GSE60862<-changenames(data=GSE60862, anno=cbind(GSE60862_anno$ID, symbols))

################################
##########GSE36192
###############################

gds<-getGEO("GSE36192",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging", AnnotGPL = TRUE)

GSE36192<-exprs(gds$GSE36192_series_matrix.txt.gz)
GSE36192_meta<-pData(gds$GSE36192_series_matrix.txt.gz)
GSE36192_anno<-fData(gds$GSE36192_series_matrix.txt.gz)

GSE36192<-changenames(data=GSE36192, anno=cbind(GSE36192_anno$ID, GSE36192_anno$`Gene symbol`))

################################
##########GSE25219
###############################

gds<-getGEO("GSE25219",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging", AnnotGPL = TRUE)

GSE25219<-exprs(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)
GSE25219_meta<-pData(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)
GSE25219_anno<-fData(gds$`GSE25219-GPL5175_series_matrix.txt.gz`)

symbols<-gsub('([^/]*)//([^/]*)//.*', '\\2', GSE25219_anno$gene_assignment)
symbols<-gsub(' ', '', symbols)

GSE25219<-changenames(data=GSE25219, anno=cbind(GSE25219_anno$ID, symbols))

##########
#probabilmente i campioni GSE60862 sono tutti in GSE36192, usare solo  GSE36192 che è più grande
match(paste(GSE60862_meta$`age at death (in years):ch1`,GSE60862_meta$`post-mortem interval (in hours):ch1`, tolower(GSE60862_meta$`Sex:ch1`))[GSE60862_meta$`brain region:ch1`=="frontal cortex"], paste(GSE36192_meta$`age (y):ch1`,GSE36192_meta$`pmi (hr):ch1`, GSE36192_meta$`gender:ch1`)[GSE36192_meta$`tissue:ch1`=="frontal cortex"])

GSE36192<-log2(GSE36192)


pca<-prcomp(t(GSE60862[,which(GSE60862_meta$`brain region:ch1`=="frontal cortex")]))
df<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2], GSE60862_meta[which(GSE60862_meta$`brain region:ch1`=="frontal cortex"),])
ggplot(df, aes(PC1, PC2, colour=brain.bank.ch1))+geom_point()

library(fgsea)
library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])



  sample<-GSE60862_meta$geo_accession[which(GSE60862_meta$`brain region:ch1`=="frontal cortex")]
  data<-GSE60862[, sample]
  age<-GSE60862_meta$`age at death (in years):ch1`[which(GSE60862_meta$`brain region:ch1`=="frontal cortex")]
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

  p<-ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
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

  p<-ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

  pdf("Results/Figures/GSEA_GSE25219.pdf", 4, 4)
  print(p)
  dev.off()

  
