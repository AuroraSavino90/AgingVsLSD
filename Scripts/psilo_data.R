psilo<-read.csv("../Data/GSE209859_counts.csv.gz")

library(DESeq2)
names_psilo<-psilo[,1]
psilo<-psilo[,-c(1:6)]
rownames(psilo)<-names_psilo
rownames(psilo)<-gsub("\\..*","",rownames(psilo))

library(biomaRt)
ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(psilo),
           mart = ensembl.mouse)

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

psilo<-changenames(data=psilo, anno=ids)

###normalization
#GSE179379<-GSE179379[which(rowSums(GSE179379>=5)>=5),]
region<-rep(NA, ncol(psilo))
region[grep("Hypothalamus", colnames(psilo))]<-"Hypothalamus"
region[grep("Frontal", colnames(psilo))]<-"Cortex"

treatment<-rep(NA, ncol(psilo))
treatment[grep("Psilocybin", colnames(psilo))]<-"Psilocybin"
treatment[grep("Saline", colnames(psilo))]<-"Saline"

time<-rep(NA, ncol(psilo))
time[grep("weeks", colnames(psilo))]<-"4weeks"
time[grep("hours", colnames(psilo))]<-"3hours"

anno<-data.frame(region=region,treatment=treatment, time=time)

dds_FC3hours <- DESeqDataSetFromMatrix(countData = psilo[,anno$time=="3hours" & anno$region=="Cortex"],
                              colData = anno[anno$time=="3hours" & anno$region=="Cortex",],
                              design= ~ treatment)
dds_FC3hours <- DESeq(dds_FC3hours)
DE_FC3hours <- results(dds_FC3hours, contrast = c("treatment","Psilocybin", "Saline"))

dds_FC4weeks <- DESeqDataSetFromMatrix(countData = psilo[,anno$time=="4weeks" & anno$region=="Cortex"],
                                       colData = anno[anno$time=="4weeks" & anno$region=="Cortex",],
                                       design= ~ treatment)
dds_FC4weeks <- DESeq(dds_FC4weeks)
DE_FC4weekss <- results(dds_FC4weeks, contrast = c("treatment","Psilocybin", "Saline"))

load("../RData/Alldata_7Jul.RData")
homologs<-read.csv("../Data/Human rat homologs.txt")
load("../RData/DE_GSE179379.RData")


genes_up_LSD<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn_LSD<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]

genes_up_psilo<-rownames(DE_FC3hours)[which(DE_FC3hours$log2FoldChange>0 & DE_FC3hours$pvalue<0.05)]
genes_dn_psilo<-rownames(DE_FC3hours)[which(DE_FC3hours$log2FoldChange<0 & DE_FC3hours$pvalue<0.05)]

genes_up_psilo<-rownames(DE_FC3hours)[which(DE_FC3hours$log2FoldChange>0 & DE_FC3hours$pvalue<0.05)]
genes_dn_psilo<-rownames(DE_FC3hours)[which(DE_FC3hours$log2FoldChange<0 & DE_FC3hours$pvalue<0.05)]


genes<-intersect(rownames(DE_GSE179379)[which(DE_GSE179379$pvalue<0.05)],
          rownames(DE_FC4weekss)[which(DE_FC4weekss$pvalue<0.05)]
)
plot(DE_GSE179379[genes, "log2FoldChange"], DE_FC4weekss[genes, "log2FoldChange"])

intersect(genes_up_LSD, genes_up_psilo)
intersect(genes_dn_LSD, genes_dn_psilo)

genes_up<-unique(homologs[which(homologs[,2] %in% genes_up_psilo),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn_psilo),3])


