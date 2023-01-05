.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")

load("RData/DE_GSE179379.RData")
load("RData/Alldata_7Jul.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

##########################
#### GSEA GO
##########################

library(clusterProfiler)
library(fgsea)
library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE30272"))#no raw
dat<-setdiff(dat, c("GSE5388"))#outlier
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
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
gsea_GO<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  
  gsea_GO[[n]]<-gseGO(forgesea,
                      ont = "BP",
                      "org.Hs.eg.db",
                      keyType = "SYMBOL",
                      exponent = 1,
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 1e-10,
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      verbose = TRUE,
                      seed = FALSE,
                      by = "fgsea"
  )
  
}
save(gsea_GO, file="Results/RData/gsea_GO.RData")

##########
## GO for LSD
#########

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down<- enrichGO(gene = DEG_down,
                    keyType="SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")

ego_up<- enrichGO(gene = DEG_up,
                  keyType="SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")


##################################
##### showing paths in common
##################################
load("Results/RData/gsea_GO.RData")

allpaths<-c()
for(i in 1:length(gsea_GO)){
  allpaths<-union(allpaths, gsea_GO[[i]]$ID[gsea_GO[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allpaths_mat[gsea_GO[[i]]$ID[which(gsea_GO[[i]]$p.adjust<0.05)],i]<-1
}


DE_GO_LSD<-union(ego_up$ID, ego_down$ID)

allpaths_mat_sel<-allpaths_mat[which(rownames(allpaths_mat) %in% DE_GO_LSD),]

library(pheatmap)
shared_paths<-allpaths_mat_sel[rowSums(allpaths_mat_sel)>7, ]
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allpaths_mat[gsea_GO[[i]]$ID[which(gsea_GO[[i]]$p.adjust<0.05)],i]<-gsea_GO[[i]]$NES[which(gsea_GO[[i]]$p.adjust<0.05)]
}
allpaths_mat_sel<-allpaths_mat[which(rownames(allpaths_mat) %in% DE_GO_LSD),]

hist(rowSums(sign(allpaths_mat_sel)))

toplot<-allpaths_mat_sel[rowSums(sign(allpaths_mat_sel))<(-6),]
rownames(toplot)<-ego_up$Description[match(rownames(toplot), ego_up$ID)]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot,   cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T)


library(metap)
allpaths<-c()
for(i in 1:length(gsea_GO)){
  allpaths<-union(allpaths, gsea_GO[[i]]$ID)
}

allp_mat<-matrix(0,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allp_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allp_mat[gsea_GO[[i]]$ID,i]<-gsea_GO[[i]]$pvalue
}

allNES_mat<-matrix(0,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allNES_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allNES_mat[gsea_GO[[i]]$ID,i]<-gsea_GO[[i]]$NES
}

DE_GO_LSD<-union(ego_up$ID, ego_down$ID)

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in ego_up$ID[which(ego_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(allNES_mat_sel[up, ])==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pup)<-ego_up$ID[which(ego_up$ID %in% rownames(allNES_mat_sel))]

pdn<-list()
ind<-0
for(dn in ego_down$ID[which(ego_down$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[dn, ]
  
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(allNES_mat_sel[dn, ])==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pdn)<-ego_down$ID[which(ego_down$ID %in% rownames(allNES_mat_sel))]
pup<-p.adjust(pup, method="fdr")
pdn<-p.adjust(pdn, method="fdr")
pup[which(pup<0.05)]
pdn[which(pdn<0.05)]

path_sel<-c(names(pdn[which(pdn<10^(-2))]))
toplot<-allNES_mat[path_sel,]
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names
rownames(toplot)<-c(ego_up$Description, ego_down$Description)[match(rownames(toplot), c(ego_up$ID, ego_down$ID))]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

anno_p<-data.frame(p= -log10(pdn[path_sel]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order(anno_p$p, decreasing = T),]
anno_p<-data.frame(p= -log10(pdn[path_sel][order(anno_p$p, decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GO_shared_dn.pdf",8,5)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

path_sel<-c(names(pup[which(pup<10^(-26))]))
toplot<-allNES_mat[path_sel,]
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names
rownames(toplot)<-c(ego_up$Description, ego_down$Description)[match(rownames(toplot), c(ego_up$ID, ego_down$ID))]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

anno_p<-data.frame(p= -log10(pup[path_sel]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order(anno_p$p, decreasing = T),]
anno_p<-data.frame(p= -log10(pup[path_sel][order(anno_p$p, decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GO_shared_up.pdf",10,15)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

library(msigdbr)
m_df = msigdbr(species = "Homo sapiens", category = "C5") #Gets the list of functional categories. I never remember the nomenclature, but you find what “C6” is on the MSigDB site.
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") 
go_list <- getBM(attributes=c("go_id", "hgnc_symbol"),
                 filters = "go",
                 values = c("GO:0050803"),
                 mart=mart)

#"GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"
genes_go<-unique(m_list["GOBP_DENDRITIC_SPINE_DEVELOPMENT"])
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p<-list()
df<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  p[[n]]<-plotEnrichment(unlist(unique(genes_go)), forgesea)
  df[[n]]<-data.frame(scaled=range01(p[[n]]$data[,1]), p[[n]]$data, dir=rep("GO:0050807", nrow(p[[n]]$data)), dataset=rep(dat_names[n],nrow(p[[n]]$data)))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}

ggplot(df_tot, aes(scaled, y, colour=dataset))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())


###which genes of FOS and EGR1?
GO_mat<-matrix(0,nrow=length(unique(unlist(genes_go))), ncol=length(genes_cor))
rownames(GO_mat)<-unique(unlist(genes_go))
for(i in 1:length(genes_cor)){
  GO_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go))))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go)))),]
}

pheatmap(GO_mat)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot, cellwidth=5, cellheight=5, breaks=myBreaks, color = myColor, keep.dendro=T)
pheatmap(toplot)

DE_GSE179379_up<-DE_GSE179379
rownames(DE_GSE179379_up)<-toupper(rownames(DE_GSE179379))

cor.test(rowMeans(GO_mat), -log10(DE_GSE179379_up[rownames(GO_mat),"pvalue"])*sign(DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"]))
df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=-log10(DE_GSE179379_up[rownames(GO_mat),"pvalue"])*sign(DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"]))
df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

library(ggrepel)
pdf("Results/Figures/DENDRITE_scatter.pdf",6,6)
ggplot(df, aes(x=aging, y=LSD))+geom_point()+geom_label_repel(label=df$gene)+geom_smooth(method="lm")+theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
dev.off()


#"GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"
genes_go<-unique(m_list["GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"])

p<-list()
df<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  p[[n]]<-plotEnrichment(unlist(unique(genes_go)), forgesea)
  df[[n]]<-data.frame(scaled=range01(p[[n]]$data[,1]), p[[n]]$data, dir=rep("GO:0050807", nrow(p[[n]]$data)), dataset=rep(dat_names[n],nrow(p[[n]]$data)))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}

ggplot(df_tot, aes(scaled, y, colour=dataset))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())


###which genes of FOS and EGR1?
GO_mat<-matrix(0,nrow=length(unique(unlist(genes_go))), ncol=length(genes_cor))
rownames(GO_mat)<-unique(unlist(genes_go))
for(i in 1:length(genes_cor)){
  GO_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go))))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go)))),]
}

pheatmap(GO_mat)

cor.test(rowMeans(GO_mat), -log10(DE_GSE179379_up[rownames(GO_mat),"pvalue"])*sign(DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"]))
df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

library(ggrepel)
pdf("Results/Figures/SYNAPSE_scatter.pdf",,6)
ggplot(df, aes(x=aging, y=LSD))+geom_point()+geom_label_repel(label=df$gene)+geom_smooth(method="lm")+theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
dev.off()

genes_go<-unique(m_list["GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS"])

p<-list()
df<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  p[[n]]<-plotEnrichment(unlist(unique(genes_go)), forgesea)
  df[[n]]<-data.frame(scaled=range01(p[[n]]$data[,1]), p[[n]]$data, dir=rep("GO:0050807", nrow(p[[n]]$data)), dataset=rep(dat_names[n],nrow(p[[n]]$data)))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}

ggplot(df_tot, aes(scaled, y, colour=dataset))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())


###which genes of FOS and EGR1?
GO_mat<-matrix(0,nrow=length(unique(unlist(genes_go))), ncol=length(genes_cor))
rownames(GO_mat)<-unique(unlist(genes_go))
for(i in 1:length(genes_cor)){
  GO_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go))))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go)))),]
}

pheatmap(GO_mat)

cor.test(rowMeans(GO_mat), -log10(DE_GSE179379_up[rownames(GO_mat),"pvalue"])*sign(DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"]))
df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

library(ggrepel)
pdf("Results/Figures/PROTEASOME_scatter.pdf",5,6)
ggplot(df, aes(x=aging, y=LSD))+geom_point()+geom_label_repel(label=df$gene)+geom_smooth(method="lm")+theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
dev.off()

##############
#### enrichr
##############
#NOT RERUN AFTER CHANGING GTEX
encode <- read.gmt("Data/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")
#encode <- read.gmt("Data/ENCODE_TF_ChIP-seq_2015.txt")
#encode<-encode[grep("ENCODE", encode[,1]),]

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_TF<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_TF[[n]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=encode, pvalueCutoff = 1, maxGSSize = 10000)

}

#save(fgsea_TF, file="Results/RData/fgsea_TF_ENCODE.RData")
save(fgsea_TF, file="Results/RData/fgsea_TF.RData")

load("Results/RData/fgsea_TF.RData")
fgsea_TF<-fgsea_TF[-which(dat_names=="GSE5388")]

allp<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allp)<-unique(encode[,1])
allNES<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allNES)<-unique(encode[,1])
for(i in 1:length(fgsea_TF)){
  allp[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$p.adjust
  allNES[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$NES
}
allNES[allp>0.05]<-NA
hist(rowSums(allNES>0, na.rm=T))
hist(rowSums(allNES<0, na.rm=T))


DEG_down_conv<-unique(homologs[which(homologs[,2] %in% DEG_down),3])

DEG_up_conv<-unique(homologs[which(homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down<- enricher(gene = DEG_down_conv,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    universe=universe, 
                    TERM2GENE = encode)

TF_up<- enricher(gene = DEG_up_conv,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  universe=universe, TERM2GENE = encode)


######evaluation collapsing the pvalue

library(metap)
allpaths<-c()
for(i in 1:length(fgsea_TF)){
  allpaths<-union(allpaths, fgsea_TF[[i]]$ID)
}

allp_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_TF))
rownames(allp_mat)<-allpaths
for(i in 1:length(fgsea_TF)){
  allp_mat[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$pvalue
}

allNES_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_TF))
rownames(allNES_mat)<-allpaths
for(i in 1:length(fgsea_TF)){
  allNES_mat[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$NES
}

DE_GO_LSD<-union(TF_up$ID, TF_down$ID)

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in TF_up$ID[which(TF_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, 14)
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
  
  istwo <- rep(T, 14)
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
toplot<-allNES_mat[path_sel,]

#path_sel<-c("FOS ENCODE", "EGR1 ENCODE")
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names[-which(dat_names=="GSE5388")]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pdf("Results/Figures/TF_shared.pdf",8,4)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)
dev.off()
              

FOS<-encode$gene[encode$term=="FOS ENCODE"]
EGR1<-encode$gene[encode$term=="EGR1 ENCODE"]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p1<-list()
p2<-list()
df<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  p1[[n]]<-plotEnrichment(unique(FOS), forgesea)
  p2[[n]]<-plotEnrichment(unique(EGR1), forgesea)
  df[[n]]<-rbind.data.frame(data.frame(scaled=range01(p1[[n]]$data[,1]), p1[[n]]$data, path=rep("FOS", nrow(p1[[n]]$data)), dataset=rep(dat_names[n],nrow(p1[[n]]$data))),
                            data.frame(scaled=range01(p2[[n]]$data[,1]), p2[[n]]$data, path=rep("EGR1", nrow(p2[[n]]$data)), dataset=rep(dat_names[n],nrow(p2[[n]]$data))))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}
df_tot$dataset<-factor(df_tot$dataset, levels=dat_names)

rank_datasets<-dat_names[-5][order(allNES_mat["FOS ENCODE",]-allNES_mat["EGR1 ENCODE",], decreasing=T)]
df_tot<-df_tot[df_tot$dataset!="GSE5388",]
df_tot$dataset<-factor(df_tot$dataset, levels=rank_datasets)

pdf("Results/Figures/TF_shared_top_GSEA.pdf",8,8)
ggplot(df_tot, aes(scaled, y, colour=path))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  facet_wrap(.~dataset)
dev.off()

pdf("Results/Figures/TF_shared_top.pdf",8,4)
pheatmap(toplot[c("FOS ENCODE", "EGR1 ENCODE"), rank_datasets], cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, cluster_cols = F)
dev.off()

###which genes of FOS and EGR1?
FOS_mat<-matrix(0,nrow=length(unique(FOS)), ncol=length(genes_cor))
rownames(FOS_mat)<-unique(FOS)
for(i in 1:length(genes_cor)){
  FOS_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(FOS))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(FOS)),]
}

pheatmap(FOS_mat)

toplot<-FOS_mat[,-5]
toplot2<-apply(-toplot,2, rank)
rownames(toplot2)<-rownames(toplot)
colnames(toplot2)<-colnames(toplot)
toplot<-toplot2

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot, cellwidth=5, cellheight=5, breaks=myBreaks, color = myColor, keep.dendro=T)
pheatmap(toplot)


rankprod<-apply(toplot2,1,prod)
names(rankprod)<-rownames(toplot2)
sort(rankprod, decreasing=F)[1:10]
toplot2["AGFG2",]
FOS_mat["AGFG2",]


EGR1_mat<-matrix(0,nrow=length(unique(EGR1)), ncol=length(genes_cor))
rownames(EGR1_mat)<-unique(EGR1)
for(i in 1:length(genes_cor)){
  EGR1_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(EGR1))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(EGR1)),]
}

toplot<-EGR1_mat[,-5]
toplot2<-apply(toplot,2, rank)
rownames(toplot2)<-rownames(toplot)
colnames(toplot2)<-colnames(toplot)
toplot<-toplot2

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot, cellwidth=5, cellheight=5, breaks=myBreaks, color = myColor, keep.dendro=T)
pheatmap(toplot)

rankprod<-apply(toplot2,1,prod)
names(rankprod)<-rownames(toplot2)
sort(rankprod, decreasing=F)[1:10]


################
#### ENCODE
################
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

encode <- read.gmt("Data/ENCODE_TF_ChIP-seq_2015.txt")
#encode<-encode[grep("ENCODE", encode[,1]),]

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_TF<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_TF[[n]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=encode, pvalueCutoff = 1, maxGSSize = 10000)
  
}

save(fgsea_TF, file="Results/RData/fgsea_TF_ENCODE.RData")

load("Results/RData/fgsea_TF_ENCODE.RData")

allp<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allp)<-unique(encode[,1])
allNES<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allNES)<-unique(encode[,1])
for(i in 1:length(fgsea_TF)){
  allp[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$p.adjust
  allNES[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$NES
}
allNES[allp>0.05]<-NA
hist(rowSums(allNES>0, na.rm=T))
hist(rowSums(allNES<0, na.rm=T))


DEG_down_conv<-unique(homologs[which(homologs[,2] %in% DEG_down),3])

DEG_up_conv<-unique(homologs[which(homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down<- enricher(gene = DEG_down_conv,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe=universe, 
                   TERM2GENE = encode)

TF_up<- enricher(gene = DEG_up_conv,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 universe=universe, TERM2GENE = encode)


######evaluation collapsing the pvalue

library(metap)
allpaths<-c()
for(i in 1:length(fgsea_TF)){
  allpaths<-union(allpaths, fgsea_TF[[i]]$ID)
}

allp_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_TF))
rownames(allp_mat)<-allpaths
for(i in 1:length(fgsea_TF)){
  allp_mat[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$pvalue
}

allNES_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_TF))
rownames(allNES_mat)<-allpaths
for(i in 1:length(fgsea_TF)){
  allNES_mat[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$NES
}

DE_GO_LSD<-union(TF_up$ID, TF_down$ID)

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in TF_up$ID[which(TF_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, 15)
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
  
  istwo <- rep(T, 15)
  toinvert <- ifelse(sign(allNES_mat_sel[dn, ])==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pdn)<-TF_down$ID[which(TF_down$ID %in% rownames(allNES_mat_sel))]

genes<-c("CHD1", "FOXP2", "CEBPD", "TAF1", "SREBF2", "NFE2","NR3C1", "EP300", "GATA2",
         "NR2C2", "RXRA", "BCL3", "SUPT20H", "RFX5")




###which genes of FOS and EGR1?
genes_mat<-matrix(0,nrow=length(unique(genes)), ncol=length(genes_cor))
rownames(genes_mat)<-unique(genes)
for(i in 1:length(genes_cor)){
  genes_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(genes))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(genes)),]
}

pheatmap(genes_mat)


toplot<-genes_mat
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1


anno2<-data.frame(LSDdir=as.factor(sign(DE_GSE179379_up[genes, "log2FoldChange"])), LSDp=-log10(DE_GSE179379_up[genes, "pvalue"]))
rownames(anno2)<-genes

pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, cluster_rows = F, annotation_row = anno2)

path_sel<-c(names(pup[which(pup<0.05)]),names(pdn[which(pdn<0.05)]))
toplot<-allNES_mat[path_sel,]

#path_sel<-c("FOS ENCODE", "EGR1 ENCODE")
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pdf("Results/Figures/TF_shared_ENCODE.pdf",8,6)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)
dev.off()

anno<-data.frame(LSDsetp=c(-log10(TF_up$pvalue)[which(TF_up$ID %in% names(pup)[pup<0.05])], -log10(TF_down$pvalue)[which(TF_down$ID %in% names(pdn)[pdn<0.05])]),
                 AgingSetp= -log10(unlist(c(pup[pup<0.05], pdn[pdn<0.05]))),
                 LSDTFdir=as.factor(sign(DE_GSE179379_up[genes, "log2FoldChange"])), LSDpTF=-log10(DE_GSE179379_up[genes, "pvalue"]),
                 TFAge=rowSums(genes_mat>0))
rownames(anno)<-c(TF_up$ID[which(TF_up$ID %in% names(pup)[pup<0.05])], TF_down$ID[which(TF_down$ID %in% names(pdn)[pdn<0.05])])

pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, annotation_row = anno)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, annotation_row = anno)

toplot_sel<-toplot[names(pup[which(pup<0.05)]),]
pheatmap(toplot_sel, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, annotation_row = anno[names(pup[which(pup<0.05)]),])

toplot_sel<-toplot[names(pdn[which(pdn<0.05)]),]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot_sel), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot_sel), na.rm=T)/paletteLength, max(unlist(toplot_sel), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot_sel, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, annotation_row = anno[names(pdn[which(pdn<0.05)]),])



FOS<-encode$gene[encode$term=="TAF1 A549 hg19"]
EGR1<-encode$gene[encode$term=="BCL3 GM12878 hg19"]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p1<-list()
p2<-list()
df<-list()
genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-sort(forgesea, decreasing=T)
  p1[[n]]<-plotEnrichment(unique(FOS), forgesea)
  p2[[n]]<-plotEnrichment(unique(EGR1), forgesea)
  df[[n]]<-rbind.data.frame(data.frame(scaled=range01(p1[[n]]$data[,1]), p1[[n]]$data, path=rep("TAF1", nrow(p1[[n]]$data)), dataset=rep(dat_names[n],nrow(p1[[n]]$data))),
                            data.frame(scaled=range01(p2[[n]]$data[,1]), p2[[n]]$data, path=rep("BCL3", nrow(p2[[n]]$data)), dataset=rep(dat_names[n],nrow(p2[[n]]$data))))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}
df_tot$dataset<-factor(df_tot$dataset, levels=dat_names)

rank_datasets<-dat_names[order(-allNES_mat["TAF1 A549 hg19",]+allNES_mat["BCL3 GM12878 hg19",], decreasing=T)]
df_tot$dataset<-factor(df_tot$dataset, levels=rank_datasets)

ggplot(df_tot, aes(scaled, y, colour=path))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  facet_wrap(.~dataset)

pdf("Results/Figures/TF_shared_top_GSEA_ENCODE.pdf",8,8)
ggplot(df_tot, aes(scaled, y, colour=path))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("TAF1"="purple", "BCL3"="orange"))+facet_wrap(~dataset, scale="free_x")
dev.off()

pdf("Results/Figures/TF_shared_top_ENCODE.pdf",8,4)
pheatmap(toplot[c("TAF1 A549 hg19", "BCL3 GM12878 hg19"), rank_datasets], cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F, cluster_cols = F)
dev.off()

###which genes of FOS and EGR1?
FOS_mat<-matrix(0,nrow=length(unique(FOS)), ncol=length(genes_cor))
rownames(FOS_mat)<-unique(FOS)
for(i in 1:length(genes_cor)){
  FOS_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(FOS))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(FOS)),]
}

pheatmap(FOS_mat)

toplot<-FOS_mat
toplot2<-apply(-toplot,2, rank)
rownames(toplot2)<-rownames(toplot)
colnames(toplot2)<-colnames(toplot)
toplot<-toplot2

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot, cellwidth=5, cellheight=5, breaks=myBreaks, color = myColor, keep.dendro=T)
pheatmap(toplot)

##########################
###### MSIGDB
###########################
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% 
  dplyr::select(gs_name, gene_symbol)

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_MsigdbC3<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_MsigdbC3[[n]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
  
}

save(fgsea_MsigdbC3, file="Results/RData/fgsea_MsigdbC3TFT.RData")

load("Results/RData/fgsea_MsigdbC3TFT.RData")
m_df<-data.frame(m_df)
allp<-matrix(NA, nrow=length(unique(m_df[,1])), ncol=length(fgsea_MsigdbC3))
rownames(allp)<-unique(m_df[,1])
allNES<-matrix(NA, nrow=length(unique(m_df[,1])), ncol=length(fgsea_MsigdbC3))
rownames(allNES)<-unique(m_df[,1])
for(i in 1:length(fgsea_MsigdbC3)){
  allp[fgsea_MsigdbC3[[i]]$ID,i]<-fgsea_MsigdbC3[[i]]$p.adjust
  allNES[fgsea_MsigdbC3[[i]]$ID,i]<-fgsea_MsigdbC3[[i]]$NES
}
allNES[allp>0.05]<-NA
hist(rowSums(allNES>0, na.rm=T))
hist(rowSums(allNES<0, na.rm=T))


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
for(i in 1:length(fgsea_MsigdbC3)){
  allpaths<-union(allpaths, fgsea_MsigdbC3[[i]]$ID)
}

allp_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC3))
rownames(allp_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC3)){
  allp_mat[fgsea_MsigdbC3[[i]]$ID,i]<-fgsea_MsigdbC3[[i]]$pvalue
}

allNES_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC3))
rownames(allNES_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC3)){
  allNES_mat[fgsea_MsigdbC3[[i]]$ID,i]<-fgsea_MsigdbC3[[i]]$NES
}

DE_GO_LSD<-union(TF_up$ID, TF_down$ID)

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in TF_up$ID[which(TF_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, 14)
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
  
  istwo <- rep(T, 14)
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
toplot<-allNES_mat[path_sel,]

#path_sel<-c("FOS ENCODE", "EGR1 ENCODE")
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-dat_names

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)


toplot_sel<-toplot[names(pup[order(unlist(pup))[1:20]]),]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot_sel), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot_sel), na.rm=T)/paletteLength, max(unlist(toplot_sel), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot_sel, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)

pheatmap(-log10(allp_mat[grep("E4F1", rownames(allp_mat)),])*sign(allNES_mat[grep("E4F1", rownames(allp_mat)),]), breaks=myBreaks, color = myColor)

toplot_sel<-toplot[names(pdn[which(unlist(pdn)<0.05)]),]

DE_GSE179379_up<-DE_GSE179379
rownames(DE_GSE179379_up)<-toupper(rownames(DE_GSE179379))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot_sel), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot_sel), na.rm=T)/paletteLength, max(unlist(toplot_sel), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(toplot_sel, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F)


#########################
##########################
######### Is LDS reverting aging?
##########################
##########################

################################
####### gsea lsd vs aging
###############################


library(fgsea)
library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
dat<-setdiff(dat, c("GSE5388"))#outlier

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
rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_res<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){
  
  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea, nPermSimple=100000)
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
  
  p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
  
}


df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}


fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/GSEA_all_datasets_RAW.pdf", 8, 8)
ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x")
dev.off()

save(fgsea_res, file="Results/RData/fgsea_res.RData")

#############################
dat<-dat_names
binned<-matrix(nrow=length(dat), ncol=100)
rownames(binned)<-dat
age_median<-c()
age_min<-c()
age_max<-c()
age_var<-c()
for(dd in dat){
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens")]))
  age_median<-c(age_median, median(age, na.rm = T))
  age_min<-c(age_min, min(age, na.rm = T))
  age_max<-c(age_max, max(age, na.rm = T))
  age_var<-c(age_var, var(age, na.rm = T))
  for(bin in 1:100){
    age_bin<-c(((bin-1)*365+1):(bin*365))
    if(length(which(age %in% intersect(age_bin, age)))>0){
      binned[dd,bin]<-length(which(age %in% intersect(age_bin, age)))
    } else {
      binned[dd,bin]<-0
    }
    
  }
}
colnames(binned)<-c(1:100)

score_dn<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])})
score_up<-lapply(fgsea_res, function(x){log10(x[2,3])*sign(x[2,6])})
score_dn<-lapply(fgsea_res, function(x){x[1,6]})
score_up<-lapply(fgsea_res, function(x){x[2,6]})

plot(unlist(score_dn), age_max)
plot(unlist(score_up), age_max)
plot(unlist(score_up)+unlist(score_dn), age_max)



#########################################
#########################################
#########################################

age_range<-c()
for(dd in dat){
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd)]))
  age_range<- c(age_range, max(age, na.rm=T)-min(age, na.rm=T))
}
names(age_range)<-dat

##############################


#######################
### random drugs
######################

load("Results/RData/CMAP_NPC_mean.RData")
load("Results/RData/compounds_NPC.RData")

fgsea_res_drugs<-list()
for(comp in 1:ncol(data_mean_all)){
  print(comp)
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE5388"))
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
  genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  genes_up<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=T)[1:length(genes_up)]]
  genes_dn<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=F)[1:length(genes_dn)]]
  
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
    fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea, nPermSimple=100000)
    
  }
  
}
save(fgsea_res_drugs, file="fgsea_res_drugs_RAW.RData")

#########################
####### EVALUATION
########################
load("Results/RData/compounds_NPC.RData")
load("Results/RData/fgsea_res_drugs_RAW.RData")


fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})


pvals_dn<-unlist(lapply(fgsea_res, function(x){x[1,3]}))
pvals_dn_rand<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,3]*ifelse(sign(x[1,6])==(-1), NA, 1)}))})
pvals_dn_rand_mat<-matrix(nrow=length(pvals_dn_rand), unlist(pvals_dn_rand), byrow=T)

df<-data.frame(p=unlist(pvals_dn_rand), dataset=rep(dat_names, each=length(pvals_dn_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_dn, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_p_dn_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

better_dn<-lapply(pvals_dn_rand, function(x){x<pvals_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
(colSums(better_dn_mat, na.rm=T)/193)[order(unlist(fgsea_rank), decreasing=T)]
compounds[order(rowSums(better_dn_mat, na.rm=T), decreasing = T)][1:10]

pvals_up<-unlist(lapply(fgsea_res, function(x){x[2,3]}))
pvals_up_rand<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,3]*ifelse(sign(x[2,6])==(1), NA, 1)}))})

df<-data.frame(p=unlist(pvals_up_rand), dataset=rep(dat_names, each=length(pvals_up_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_up, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_p_up_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

better_up<-lapply(pvals_up_rand, function(x){x<pvals_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat, na.rm=T)/193
compounds[order(rowSums(better_up_mat, na.rm=T), decreasing = T)][1:10]

pdf("Results/Figures/Random_drugs_p_perc_RAW.pdf", 4,4)
plot(1-colSums(better_dn_mat, na.rm=T)/193, 1-colSums(better_up_mat, na.rm=T)/193, xlab="% of tests less significant than LSD (dn)",
     ylab="% of tests less significant than LSD (up)", pch=19)
dev.off()

################################
###evaluation with nes
##############################
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

pdn<-meta_dn(fgsea_res)
sump_dn<-lapply(fgsea_res_drugs, meta_dn)
sum(unlist(sump_dn)<2.2*10^(-16))
sum(unlist(sump_dn)<pdn)
compounds[which(unlist(sump_dn)<pdn)]

df<-data.frame(p=unlist(sump_dn), type=rep("random", length(sump_dn)))
df<-rbind.data.frame(df, data.frame(p=pdn, type="LSD"))

pdf("Results/Figures/Random_drugs_fisher_dn_RAW.pdf", 4,4)
ggplot(subset(df, type="random"), aes(x=-log10(p)))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=-log10(p)), colour="red", linetype="dashed")+theme_classic()
dev.off()

pup<-meta_up(fgsea_res)
sump_up<-lapply(fgsea_res_drugs, meta_up)
sum(unlist(sump_up)<2.2*10^(-16))
sum(unlist(sump_up)<pup)

df<-data.frame(p=unlist(sump_up), type=rep("random", length(sump_up)))
df<-rbind.data.frame(df, data.frame(p=pup, type="LSD"))

pdf("Results/Figures/Random_drugs_fisher_up_RAW.pdf", 4,4)
ggplot(subset(df, type="random"), aes(x=-log10(p)))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=-log10(p)), colour="red", linetype="dashed")+theme_classic()
dev.off()


#############################
#### random sets of genes ###
#############################
set.seed(4805736257)
fgsea_res_rand<-list()
for(iter in 1:1000){
  print(iter)
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE30272"))#tmp
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
  genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  genes_up<-unique(homologs[sample(1:nrow(homologs),length(genes_up)),3])
  genes_dn<-unique(homologs[sample(1:nrow(homologs),length(genes_dn)),3])
  
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
  rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])
  
  genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
  fgsea_res_rand[[iter]]<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){
    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    forgesea<-forgesea[names(forgesea) %in% rat_homologs]
    fgsea_res_rand[[iter]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
    
  }
  
}
save(fgsea_res_rand, file="Results/RData/fgsea_res_rand_RAW.RData")

load(file="Results/RData/fgsea_res_rand_RAW.RData")

pvals_dn<-unlist(lapply(fgsea_res, function(x){x[1,3]}))
pvals_dn_rand<-lapply(fgsea_res_rand, function(y){unlist(lapply(y, function(x){x[1,3]}))})

better_dn<-lapply(pvals_dn_rand, function(x){x<pvals_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
colSums(better_dn_mat)

pvals_up<-unlist(lapply(fgsea_res, function(x){x[2,3]}))
pvals_up_rand<-lapply(fgsea_res_rand, function(y){unlist(lapply(y, function(x){x[2,3]}))})

df<-data.frame(p=unlist(pvals_dn_rand), dataset=rep(dat_names, each=length(pvals_dn_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_dn, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
pdf("Results/Figures/Random_genes_p_dn_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

df<-data.frame(p=unlist(pvals_up_rand), dataset=rep(dat_names, each=length(pvals_up_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_up, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
pdf("Results/Figures/Random_genes_p_up_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

better_up<-lapply(pvals_up_rand, function(x){x<pvals_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat, na.rm=T)
compounds[order(rowSums(better_up_mat, na.rm=T), decreasing = T)][1:10]


###evaluation with nes
nes_dn<-unlist(lapply(fgsea_res, function(x){x[1,6]}))
nes_dn_rand<-lapply(fgsea_res_rand, function(y){unlist(lapply(y, function(x){x[1,6]}))})

better_dn<-lapply(nes_dn_rand, function(x){x>nes_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
colSums(better_dn_mat)

nes_up<-unlist(lapply(fgsea_res, function(x){x[2,6]}))
nes_up_rand<-lapply(fgsea_res_rand, function(y){unlist(lapply(y, function(x){x[2,6]}))})

####
better_up<-lapply(nes_up_rand, function(x){x<nes_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat)




###########Other psychedelics
DEGs_array<-c("DE_MDMA",
              "DE_GSE47541",
              "DE_GSE19914",#pups, remove
              
              #ketamine
              "DE_GSE73800_k1",
              "DE_GSE73800_k2",
              "DE_GSE73800_k4",
              "DE_GSE73800_k8",
              "DE_GSE26364",
              #PCP
              "DE_GSE73800_p1",
              "DE_GSE73800_p2",
              "DE_GSE73800_p4",
              "DE_GSE73800_p8",
              
              #DOI
              "DE_GSE14720",
              "DE_GSE23728",
              #LSD
              "DE_GSE179380")

library(fgsea)
library(ggplot2)
for(DE in DEGs_array){
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
  genes_up<-toupper(genes_up)
  genes_dn<-toupper(genes_dn)
  
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE30272"))#tmp
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  
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
  
  df_tot<-df[[1]]
  for(n in 2:length(genes_cor)){
    df_tot<-rbind.data.frame(df_tot, df[[n]])
  }
  
  
  fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})
  
  df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
  
  png(paste("GSEA_all_datasets", DE, "_RAW.png", sep=""), res=300, 4000, 4000)
  print(ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
}



DEGs_seq<-c("DE_GSE179379",
            "DE_GSE161626",
            "DE_DMT",
            "DE_pharm")



for(DE in DEGs_seq){
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  genes_up<-toupper(genes_up)
  genes_dn<-toupper(genes_dn)
  
  
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE30272"))#tmp
  dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
  
  
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
  
  df_tot<-df[[1]]
  for(n in 2:length(genes_cor)){
    df_tot<-rbind.data.frame(df_tot, df[[n]])
  }
  
  
  fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})
  
  df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
  
  png(paste("GSEA_all_datasets", DE, "_RAW.png", sep=""), res=300, 4000, 4000)
  print(ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
}

