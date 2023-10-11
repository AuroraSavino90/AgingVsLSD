rm(list=ls())
library(biomaRt)
library(openxlsx)
library(ggplot2)
library(DESeq2)
library(FactoMineR)

meta<-read.xlsx("Data/Samples.xlsx",1)
meta$dose[meta$dose=="ctrl"]<-0

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

#load counts
filesL1<-list.files("Data/Counts", pattern="L001")
filesL2<-list.files("Data/Counts",pattern="L002")
L1mat<-read.csv(paste("Data/Counts", filesL1[1], sep="/"), sep="\t")
for(i in filesL1){
  L1mat<-cbind(L1mat, read.csv(paste("Data/Counts", i, sep="/"), sep="\t")[,2])
}
rownames(L1mat)<-L1mat[,1]
L1mat<-L1mat[,-c(1,2)]
L2mat<-read.csv(paste("Data/Counts", filesL2[1], sep="/"), sep="\t")
for(i in filesL2){
  L2mat<-cbind(L2mat, read.csv(paste("Data/Counts", i, sep="/"), sep="\t")[,2])
}
rownames(L2mat)<-L2mat[,1]
L2mat<-L2mat[,-c(1,2)]

counts<-L1mat+L2mat

#convert to gene symbol
ensembl.mouse<- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
ids<-getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(counts),
           mart = ensembl.mouse)
counts_mgi<-changenames(data=counts, anno=ids)

#filter genes with a minimum number of counts per sample
counts_mgi_filt<-counts_mgi[rowSums(counts_mgi>=10)>=3,]

#normalize (RPM)
RPM<-t(t(counts_mgi_filt)/colSums(counts_mgi_filt))*1000000
RPM_log<-log2(RPM+1)


###only first 3 cohorts
counts_mgi_filt1<-counts_mgi[,which(meta$cohort %in% c("C1", "C2", "C3"))]
counts_mgi_filt1<-counts_mgi_filt1[rowSums(counts_mgi_filt1>=10)>=3,]

RPM1<-t(t(counts_mgi_filt1)/colSums(counts_mgi_filt1))*1000000
vargenes<-apply(RPM1,1, var)
RPM1<-RPM1[which(vargenes>0.5),]

RPM_log1<-log2(RPM1+1)

meta_filt1<-meta[which(meta$cohort %in% c("C1", "C2", "C3")),]

pca<-PCA(t(RPM_log1))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3], PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],meta_filt1)

ggplot(df, aes(x=PC1, y=PC2, colour=as.numeric(dose)))+geom_point(size=3)+scale_color_gradient(low = "blue", high = "red")
ggplot(df, aes(x=PC3, y=PC4, colour=daytime))+geom_point(size=3)
ggplot(df, aes(x=PC1, y=PC2, colour=growth))+geom_point()

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_mgi_filt1[rownames(RPM_log1),],
                                      colData = data.frame(meta_filt1),
                                      design= ~ dose+cohort)
dds <- DESeq2::DESeq(dds)
DE_500 <- DESeq2::results(dds, contrast = c("dose", "500", "0"))
DE_200 <- DESeq2::results(dds, contrast = c("dose","200","0"))
DE_100 <- DESeq2::results(dds, contrast = c("dose", "100", "0"))
DE_50 <- DESeq2::results(dds, contrast = c("dose", "50", "0"))


library(fgsea)
library(ggrepel)

load("D:/OneDrive - Htechnopole/Documents/Work/Projects/Psychedelics/AgingVsLSD/Results/RData/Alldata_20Sep.RData")
rat_homologs<-read.csv("D:/OneDrive - Htechnopole/Documents/Work/Projects/Psychedelics/AgingVsLSD/Data/Human rat homologs.txt", sep="\t")
mouse_homologs<-read.csv("D:/OneDrive - Htechnopole/Documents/Work/Projects/Psychedelics/AgingVsLSD/Data/Human mouse homologs.txt", sep="\t")

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


######################### RNA-seq
DEGs_seq<-c("DE_500", "DE_200", "DE_100", "DE_50")
cor_up_tot<-list()
cor_dn_tot<-list()
i<-0
n<-0
pdn<-list()
pup<-list()

for(DE in DEGs_seq){
  i<-i+1
  n<-n+1
homologs<-mouse_homologs
 
  
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
  
  pdf(paste("GSEA_aging_",DE,".pdf", sep=""), 8, 8)
  print(ggplot(df_tot, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
  save(fgsea_res, file=paste("fgsea_res_", DE, ".RData", sep=""))
}



##################
## collapsing GSEA p (no random drugs)
##################
p_dn<-c()
p_up<-c()
for(DE in c(DEGs_seq)){
  load(paste("fgsea_res_", DE, ".RData", sep=""))
  
  p_dn<-c(p_dn, meta_dn(fgsea_res))
  p_up<-c(p_up, meta_up(fgsea_res))
  
}
#p_dn[p_dn==0]<-10^(-200)

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_seq))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")


p_dn_rev<-c()
p_up_rev<-c()
for(DE in c(DEGs_seq)){
  load(paste("fgsea_res_", DE, ".RData", sep=""))
  
  p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
  p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))
  
}

df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c( DEGs_seq))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_seq)))


df<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=DEGs_seq)

pdf("Results/All_DE_aging.pdf",5,5)
ggplot(df, aes(x=pup_diff, y=pdn_diff, label=dataset))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))
dev.off()

##########
## GO for LSD
#########

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down<-list()
ego_up<-list()
i<-0
for(DE in DEGs_seq){
  i<-i+1
  homologs<-mouse_homologs
  
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  
  universe<-rownames(DEGs)
  
  ego_down[[i]]<- enrichGO(gene = genes_dn,
                    keyType="SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    universe=universe, OrgDb="org.Mm.eg.db")

ego_up[[i]]<- enrichGO(gene = genes_up,
                  keyType="SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  universe=universe, OrgDb="org.Mm.eg.db")

}

save(ego_up, file="Results/ego_up_valid.RData")
save(ego_down, file="Results/ego_down_valid.RData")

######################
### MSIGDB
##########################

library(msigdbr)
library(clusterProfiler)
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol)

i<-0
REACTOME_down<-list()
REACTOME_up<-list()
for(DE in DEGs_seq){
  i<-i+1
  homologs<-mouse_homologs
  
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
universe<-rownames(DEGs)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

REACTOME_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  =1,
                   qvalueCutoff  = 1,
                   universe=universe, 
                   TERM2GENE = m_df)

REACTOME_up[[i]]<- enricher(gene = genes_up,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1,
                 universe=universe, TERM2GENE = m_df)

}

save(REACTOME_up, file="Results/REACTOME_up_valid.RData")
save(REACTOME_down, file="Results/REACTOME_down_valid.RData")

View(TF_up[[2]][rownames(toplot_sel),])

allpath<-unique(c(TF_up[[1]][,1], TF_up[[2]][,1], TF_up[[3]][,1], TF_up[[4]][,1]))
mat<-matrix(nrow=length(allpath), ncol=4)
rownames(mat)<-allpath
mat[TF_up[[1]][,1],1]<-TF_up[[1]][,"pvalue"]
mat[TF_up[[2]][,1],2]<-TF_up[[2]][,"pvalue"]
mat[TF_up[[3]][,1],3]<-TF_up[[3]][,"pvalue"]
mat[TF_up[[4]][,1],4]<-TF_up[[4]][,"pvalue"]
mat<-mat[rowSums(is.na(mat))<2,]
pheatmap(-log10(mat))
sort(rowSums(mat<0.05, na.rm=T))
