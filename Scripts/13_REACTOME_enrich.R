rm(list=ls())
load("Results/RData/Alldata_20Sep.RData")
load("Results/RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
load("Results/RData/fgsea_res.RData")

library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741", "GSE5388"))

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

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
dat<-setdiff(dat, c("GSE102741", "GSE5388"))

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
#LSD
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]

##########################
###### MSIGDB
###########################
library(msigdbr)
library(clusterProfiler)
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, gene_symbol)

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_MsigdbC2CP<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_MsigdbC2CP[[n]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
  
}

save(fgsea_MsigdbC2CP, file="Results/RData/fgsea_MsigdbREACTOME.RData")

load("Results/RData/fgsea_MsigdbREACTOME.RData")
m_df<-data.frame(m_df)
allp<-matrix(NA, nrow=length(unique(m_df[,1])), ncol=length(fgsea_MsigdbC2CP))
rownames(allp)<-unique(m_df[,1])
allNES<-matrix(NA, nrow=length(unique(m_df[,1])), ncol=length(fgsea_MsigdbC2CP))
rownames(allNES)<-unique(m_df[,1])
for(i in 1:length(fgsea_MsigdbC2CP)){
  allp[fgsea_MsigdbC2CP[[i]]$ID,i]<-fgsea_MsigdbC2CP[[i]]$p.adjust
  allNES[fgsea_MsigdbC2CP[[i]]$ID,i]<-fgsea_MsigdbC2CP[[i]]$NES
}
allNES[allp>0.05]<-NA
hist(rowSums(allNES>0, na.rm=T))
hist(rowSums(allNES<0, na.rm=T))

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

DEG_down_conv<-unique(homologs[which(homologs[,2] %in% DEG_down),3])
DEG_up_conv<-unique(homologs[which(homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down<- enricher(gene = DEG_down_conv,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe=universe, 
                   TERM2GENE = m_df)

TF_up<- enricher(gene = DEG_up_conv,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 universe=universe, TERM2GENE = m_df)


######evaluation collapsing the pvalue

library(metap)
allpaths<-c()
for(i in 1:length(fgsea_MsigdbC2CP)){
  allpaths<-union(allpaths, fgsea_MsigdbC2CP[[i]]$ID)
}

allp_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP))
rownames(allp_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP)){
  allp_mat[fgsea_MsigdbC2CP[[i]]$ID,i]<-fgsea_MsigdbC2CP[[i]]$pvalue
}

allNES_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP))
rownames(allNES_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP)){
  allNES_mat[fgsea_MsigdbC2CP[[i]]$ID,i]<-fgsea_MsigdbC2CP[[i]]$NES
}

DE_GO_LSD<-union(TF_up$ID, TF_down$ID)

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in TF_up$ID[which(TF_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, length(genes_cor))
  toinvert <- ifelse(sign(allNES_mat_sel[up, ])==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pup)<-TF_up$ID[which(TF_up$ID %in% rownames(allNES_mat_sel))]

pdn<-list()
ind<-0
for(dn in TF_down$ID[which(TF_down$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[dn, ]
  
  istwo <- rep(T, length(genes_cor))
  toinvert <- ifelse(sign(allNES_mat_sel[dn, ])==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pdn)<-TF_down$ID[which(TF_down$ID %in% rownames(allNES_mat_sel))]

path_sel<-c(names(pup[which(pup<0.05)]),names(pdn[which(pdn<0.05)]))
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

pdf("Results/Figures/REACTOMEpath.pdf",15,20)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)
dev.off()

toplot_sel<-toplot[names(pup[order(unlist(pup))[1:10]]),]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot_sel), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot_sel), na.rm=T)/paletteLength, max(unlist(toplot_sel), na.rm=T), length.out=floor(paletteLength/2)))

pdf("Results/Figures/REACTOMEup_top.pdf",15,8)
pheatmap(toplot_sel, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)
dev.off()
pheatmap(-log10(allp_mat[grep("E4F1", rownames(allp_mat)),])*sign(allNES_mat[grep("E4F1", rownames(allp_mat)),]), breaks=myBreaks, color = myColor)


