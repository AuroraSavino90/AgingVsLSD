#####################################################################
##function to change gene names (e.g. from ENSEMBL to gene symbol)
###################################################################

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


##############Data

invitro_TPM<-read.csv("Data/In vitro RNA-seq_all/salmon.merged.gene_tpm.tsv", sep="\t")
rownames(invitro_TPM)<-invitro_TPM[,1]

anno<-invitro_TPM[,c(1:2)]
invitro_TPM<-changenames(invitro_TPM[,-c(1:2)], anno = anno)

invitro_TPM<-log2(invitro_TPM+1)


###########################
### quality checks: PCA
##########################
library(openxlsx)
library(FactoMineR)
library(ggplot2)
metadata<-read.xlsx("Data/In vitro RNA-seq_all/Metadata.xlsx", rowNames = T)

pca<-PCA(t(invitro_TPM))

metadata$Treatment<-as.factor(metadata$Treatment)
metadata$Replicate<-as.factor(metadata$Replicate)
metadata$Trt<-rep(c("Ct", "LSD", "Abeta", "Abeta", "Abeta",
                    "Abeta + LSD", "Abeta + LSD", "Abeta + LSD"),3)

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)

df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)

png("results/Figures/invitro_PCA_rep.png", res=300, 1800, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=Replicate))+geom_point(size=2)+theme_classic()
dev.off()

png("results/Figures/invitro_PCA_trt.png", res=300, 1800, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

colnames(invitro_TPM)<-paste(metadata$Treatment, metadata$Replicate)

#######################################################
########### DEGs
#########################################################

#from https://www.biostars.org/p/9570365/,
#from https://www.biostars.org/p/9581212/,
# tximport is scaling the counts from transcripts of Salmon of isoforms of different lengths
# and put together into gene-level. Since isoforms are different in length, tximport corrects for the
# biases in di different counts because of different lengths.
# nfcore/rnaseq does Salmon and then tximport. 
# So: use salmon.merged.gene_counts_length_scaled.tsv file (raw counts estimated after bias length correction)
# and feed this directly into DESeqDataSetFromMatrix() function and round() the numbers to get integers.
#read and prepare raw count table


invitro_counts=read.table("Data/In vitro RNA-seq_all/salmon.merged.gene_counts.tsv",header=TRUE)
rownames(invitro_counts)=invitro_counts$gene_id
anno<-invitro_counts[,c(1:2)]
invitro_counts<-changenames(invitro_counts[,-c(1:2)], anno = anno)
invitro_counts=round(invitro_counts)
save(invitro_counts, file="results/invitro_counts.RData")

rownames(metadata)<-colnames(invitro_counts)

library(sva)
adjusted <- ComBat_seq(as.matrix(invitro_counts), batch=metadata$Replicate, group=NULL)

invitro_RPM<-t(t(adjusted)/colSums(adjusted))*10000000

pca<-PCA(t(invitro_RPM))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)
df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)

png("results/Figures/invitro_PCA_rep_postcorr.png", res=300, 1800, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=Replicate))+geom_point(size=2)+theme_classic()
dev.off()

png("results/Figures/invitro_PCA_trt_postcorr12.png", res=300, 1800, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

png("results/Figures/invitro_PCA_trt_postcorr34.png", res=300, 1800, 1500)
ggplot(df, aes(x=PC3, y=PC4, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

invitro_counts<-adjusted[rowSums(adjusted>=10)>=3,]
invitro_TPM<-invitro_TPM[rownames(invitro_counts),]


############## DEGs
library(DESeq2)

invitro_counts_abeta<-invitro_counts[,metadata$Treatment %in% c("Abeta 1uM ","Abeta 500nM ","Abeta 100nM", "Ct"  )]

dds <- DESeqDataSetFromMatrix(countData = invitro_counts_abeta,
                                 colData = metadata[metadata$Treatment %in% c("Abeta 1uM ","Abeta 500nM ","Abeta 100nM", "Ct"  ),],
                                 design= ~ Treatment)
dds <- DESeq(dds)
ddsAbeta1000 <- results(dds, contrast = c("Treatment","Abeta 1uM ", "Ct"))
ddsAbeta500 <- results(dds, contrast = c("Treatment","Abeta 500nM ", "Ct"))
ddsAbeta100 <- results(dds, contrast = c("Treatment","Abeta 100nM", "Ct"))
write.xlsx(data.frame(ddsAbeta1000), file="results/Abeta1000DEGs.xlsx", rowNames=T)
write.xlsx(data.frame(ddsAbeta500), file="results/Abeta500DEGs.xlsx", rowNames=T)
write.xlsx(data.frame(ddsAbeta100), file="results/Abeta100DEGs.xlsx", rowNames=T)

thr<-0.01
Abeta_genes_dn<-unique(c(rownames(ddsAbeta100)[which(ddsAbeta100$padj<thr & ddsAbeta100$log2FoldChange<0)],
                  rownames(ddsAbeta500)[which(ddsAbeta500$padj<thr & ddsAbeta500$log2FoldChange<0)],
                  rownames(ddsAbeta1000)[which(ddsAbeta1000$padj<thr & ddsAbeta1000$log2FoldChange<0)]))
Abeta_genes_up<-unique(c(rownames(ddsAbeta100)[which(ddsAbeta100$padj<thr & ddsAbeta100$log2FoldChange>0)],
               rownames(ddsAbeta500)[which(ddsAbeta500$padj<thr & ddsAbeta500$log2FoldChange>0)],
               rownames(ddsAbeta1000)[which(ddsAbeta1000$padj<thr & ddsAbeta1000$log2FoldChange>0)]))

write.csv(Abeta_genes_up, file="results/Abeta_genes_up.csv")
write.csv(Abeta_genes_dn, file="results/Abeta_genes_dn.csv")

library(GSVA)
ssgsea_par <- ssgseaParam(as.matrix(invitro_RPM), list(up=Abeta_genes_up,dn=Abeta_genes_dn))  # all other values are default values
ssgsea <- gsva(ssgsea_par)

library(ggside)
df<-data.frame(up=ssgsea[1,], dn=ssgsea[2,], metadata)
df$Trt<-factor(df$Trt, levels=c("Ct", "LSD", "Abeta + LSD", "Abeta"))

png("results/Figures/invitro_Abeta_revert_sep_rep.png", res=300, 2500, 1600)
ggplot(df, aes(x=up, y=dn, colour=as.factor(Trt)))+geom_point(size=3)+facet_grid(~Replicate)+theme_bw()
dev.off()


png("results/Figures/invitro_Abeta_revert_sep.png", res=300, 2200, 1800)
ggplot(df, aes(x=up, y=dn, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()

invitro_RPM_sel<-invitro_RPM[c(Abeta_genes_up, Abeta_genes_dn),]
pca<-PCA(t(invitro_RPM_sel))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)
df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)

png("results/Figures/invitro_Abeta_revert_sep_PCA12.png", res=300, 2200, 1800)
ggplot(df, aes(x=PC1, y=PC2, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()

png("results/Figures/invitro_Abeta_revert_sep_PCA13.png", res=300, 2200, 1800)
ggplot(df, aes(x=PC1, y=PC3, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()