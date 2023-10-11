rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")

load("Results/RData/DE_GSE179379.RData")
load("Results/RData/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
load("Results/ego_up_valid.RData")
load("Results/ego_down_valid.RData")

##########
## GO for LSD
#########

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down[[5]]<- enrichGO(gene = DEG_down,
                    keyType="SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")

ego_up[[5]]<- enrichGO(gene = DEG_up,
                  keyType="SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")


##################################
##### statistics across aging datasets
##################################
load("Results/RData/gsea_GO.RData")

library(metap)
allpaths<-c()
for(i in 1:length(gsea_GO)){
  allpaths<-union(allpaths, gsea_GO[[i]]$Description)
}

allp_mat<-matrix(NA,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allp_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allp_mat[gsea_GO[[i]]$Description,i]<-gsea_GO[[i]]$pvalue
}

allNES_mat<-matrix(NA,nrow=length(allpaths), ncol=length(gsea_GO))
rownames(allNES_mat)<-allpaths
for(i in 1:length(gsea_GO)){
  allNES_mat[gsea_GO[[i]]$Description,i]<-gsea_GO[[i]]$NES
}


##############################
#### paths in common across LSD datasets (Up)
#############################
library(pheatmap)
library(ggplot2)
library(metap)
library(ggrepel)

allpaths_LSD<-c()
for(i in 1:length(ego_up)){
  allpaths_LSD<-union(allpaths_LSD, ego_up[[i]]$Description)
}

allp_mat_LSD<-matrix(NA,nrow=length(allpaths_LSD), ncol=length(ego_up))
rownames(allp_mat_LSD)<-allpaths_LSD
for(i in 1:length(ego_up)){
  allp_mat_LSD[ego_up[[i]]$Description,i]<-ego_up[[i]]$pvalue
}

pup_LSD<-list()
ind<-0
for(up in allpaths_LSD){
  ind<-ind+1
  istwo <- rep(F, 5)
  toinvert <- rep(F, 5)
  missing<-which(is.na(allp_mat_LSD[up, ]))
  if(length(missing)==0){
    pup_LSD[[ind]]<-sumlog(two2one(allp_mat_LSD[up, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pup_LSD[[ind]]<-sumlog(two2one(allp_mat_LSD[up, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}

pup_LSD<-p.adjust(pup_LSD, method="fdr")
names(pup_LSD)<-allpaths_LSD

DE_GO_LSD<-names(pup_LSD)[which(pup_LSD<0.05)]

##################################
##### statistics across aging datasets for paths up-regulated by LSD
##################################

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in DE_GO_LSD[which(DE_GO_LSD %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  
  istwo <- rep(T, length(gsea_GO))
  toinvert <- ifelse(sign(allNES_mat_sel[up, ])==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pup[[ind]]<-sumlog(two2one(allp_mat_sel[up, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pup)<-DE_GO_LSD[which(DE_GO_LSD %in% rownames(allNES_mat_sel))]

pup<-p.adjust(pup, method="fdr")

path_sel<-names(pup[which(pup<0.05)])
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-names(gsea_GO)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

anno_p<-data.frame(p= -log10(pup[path_sel]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order(anno_p$p, decreasing = T),]
anno_p<-data.frame(p= -log10(pup[path_sel][order(anno_p$p, decreasing = T)]))
rownames(anno_p)<-rownames(toplot)


graphics.off()
pdf("Results/Figures/GO_shared_up.pdf",15,50)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

path_sel2<-names(sort(-log10(unlist(pup_LSD[path_sel]))-log10(unlist(pup[path_sel])), decreasing=T))[1:20]
toplot<-allNES_mat[path_sel2,]
toplot<- -log10(allp_mat[path_sel2,])*sign(allNES_mat[path_sel2,])
colnames(toplot)<-names(gsea_GO)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

anno_p<-data.frame(p= -log10(pup[path_sel2]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order((-log10(unlist(pup_LSD[path_sel2]))-log10(unlist(pup[path_sel2]))),decreasing = T),]
anno_p<-data.frame(p= -log10(pup[path_sel2][order((-log10(unlist(pup_LSD[path_sel2]))-log10(unlist(pup[path_sel2]))),decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GO_shared_up_top.pdf",10,15)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

############################
### LSD plots up
##############################
path_sel<-names(pup[which(pup<0.05)])

toplot<- -log10(allp_mat_LSD[path_sel,])
colnames(toplot)<-c("50ug/kg", "100ug/kg", "200ug/kg", "500ug/kg", "Chronic")


paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

pdf("Results/Figures/GOup_allLSD_drugonly.pdf",15,100)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F,  cluster_cols = F)
dev.off()

toplot<-toplot[path_sel2,]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

toplot<-toplot[order((-log10(unlist(pup_LSD[path_sel2]))-log10(unlist(pup[path_sel2]))),decreasing = T),]
anno_p<-data.frame(p= -log10(pup_LSD[path_sel2][order((-log10(unlist(pup_LSD[path_sel2]))-log10(unlist(pup[path_sel2]))),decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GOup_top_allLSD_drugonly.pdf",15,8)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F,  cluster_cols = F)
dev.off()


df<-data.frame(logpLSD=-log10(unlist(pup_LSD[path_sel])),logpAging=-log10(unlist(pup[path_sel])),
               path=path_sel)

pdf("Results/Figures/GOup_LSDvsAging.pdf",7,7)
ggplot(df, aes(x=logpLSD, y=logpAging, label=gsub("GO_", "", path)))+geom_point()+geom_label_repel(label.size = 0.15,max.overlaps =20)+theme_classic()
dev.off()


##############################
#### paths in common across LSD datasets (Down)
#############################

library(metap)
allpaths_LSD<-c()
for(i in 1:length(ego_down)){
  allpaths_LSD<-union(allpaths_LSD, ego_down[[i]]$Description)
}

allp_mat_LSD<-matrix(NA,nrow=length(allpaths_LSD), ncol=length(ego_down))
rownames(allp_mat_LSD)<-allpaths_LSD
for(i in 1:length(ego_down)){
  allp_mat_LSD[ego_down[[i]]$Description,i]<-ego_down[[i]]$pvalue
}

pdn_LSD<-list()
ind<-0
for(dn in allpaths_LSD){
  ind<-ind+1
  istwo <- rep(F, 5)
  toinvert <- rep(F, 5)
  missing<-which(is.na(allp_mat_LSD[dn, ]))
  if(length(missing)==0){
    pdn_LSD[[ind]]<-sumlog(two2one(allp_mat_LSD[dn, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pdn_LSD[[ind]]<-sumlog(two2one(allp_mat_LSD[dn, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}

pdn_LSD<-p.adjust(pdn_LSD, method="fdr")
names(pdn_LSD)<-allpaths_LSD

DE_GO_LSD<-names(pdn_LSD)[which(pdn_LSD<0.05)]

##################################
##### statistics across aging datasets for paths dn-regulated by LSD
##################################

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pdn<-list()
ind<-0
for(dn in DE_GO_LSD[which(DE_GO_LSD %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  
  istwo <- rep(T, length(gsea_GO))
  toinvert <- ifelse(sign(allNES_mat_sel[dn, ])==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ], two = istwo, invert = toinvert))[[3]]
  } else {
    pdn[[ind]]<-sumlog(two2one(allp_mat_sel[dn, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
  }
}
names(pdn)<-DE_GO_LSD[which(DE_GO_LSD %in% rownames(allNES_mat_sel))]

pdn<-p.adjust(pdn, method="fdr")

path_sel<-names(pdn[which(pdn<0.05)])
toplot<- -log10(allp_mat[path_sel,])*sign(allNES_mat[path_sel,])
colnames(toplot)<-names(gsea_GO)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

anno_p<-data.frame(p= -log10(pdn[path_sel]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order(anno_p$p, decreasing = T),]
anno_p<-data.frame(p= -log10(pdn[path_sel][order(anno_p$p, decreasing = T)]))
rownames(anno_p)<-rownames(toplot)


graphics.off()
pdf("Results/Figures/GO_shared_dn.pdf",15,50)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

path_sel2<-names(sort(-log10(unlist(pdn_LSD[path_sel]))-log10(unlist(pdn[path_sel])), decreasing=T))[1:20]
toplot<-allNES_mat[path_sel2,]
toplot<- -log10(allp_mat[path_sel2,])*sign(allNES_mat[path_sel2,])
colnames(toplot)<-names(gsea_GO)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

anno_p<-data.frame(p= -log10(pdn[path_sel2]))
rownames(anno_p)<-rownames(toplot)
toplot<-toplot[order((-log10(unlist(pdn_LSD[path_sel2]))-log10(unlist(pdn[path_sel2]))),decreasing = T),]
anno_p<-data.frame(p= -log10(pdn[path_sel2][order((-log10(unlist(pdn_LSD[path_sel2]))-log10(unlist(pdn[path_sel2]))),decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GO_shared_dn_top.pdf",10,15)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

############################
### LSD plots dn
##############################
path_sel<-names(pdn[which(pdn<0.05)])

toplot<- -log10(allp_mat_LSD[path_sel,])
colnames(toplot)<-c("50ug/kg", "100ug/kg", "200ug/kg", "500ug/kg", "Chronic")


paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

pdf("Results/Figures/GOdn_allLSD_drugonly.pdf",15,100)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T,  cluster_rows = F,  cluster_cols = F)
dev.off()

toplot<-toplot[path_sel2,]

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))

toplot<-toplot[order((-log10(unlist(pdn_LSD[path_sel2]))-log10(unlist(pdn[path_sel2]))),decreasing = T),]
anno_p<-data.frame(p= -log10(pdn_LSD[path_sel2][order((-log10(unlist(pdn_LSD[path_sel2]))-log10(unlist(pdn[path_sel2]))),decreasing = T)]))
rownames(anno_p)<-rownames(toplot)

pdf("Results/Figures/GOdn_top_allLSD_drugonly.pdf",15,8)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F,  cluster_cols = F)
dev.off()


df<-data.frame(logpLSD=-log10(unlist(pdn_LSD[path_sel])),logpAging=-log10(unlist(pdn[path_sel])),
               path=path_sel)

pdf("Results/Figures/GOdn_LSDvsAging.pdf",7,7)
ggplot(df, aes(x=logpLSD, y=logpAging, label=gsub("GO_", "", path)))+geom_point()+geom_label_repel(label.size = 0.15,max.overlaps =20)+theme_classic()
dev.off()

#########################################

load(file="Results/genes_cor.RData")

library(msigdbr)
library(biomaRt)
library(fgsea)

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
  df[[n]]<-data.frame(scaled=range01(p[[n]]$data[,1]), p[[n]]$data, dir=rep("GO:0050807", nrow(p[[n]]$data)), dataset=rep(names(gsea_GO)[n],nrow(p[[n]]$data)))
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

DE_GSE179379_up<-DE_GSE179379
DE_GSE179379_up<-DE_GSE179379_up[!is.na(homologs[match(rownames(DE_GSE179379_up), homologs[,2]),3]),]
rownames(DE_GSE179379_up)<-homologs[match(rownames(DE_GSE179379_up), homologs[,2]),3]

df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

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
  df[[n]]<-data.frame(scaled=range01(p[[n]]$data[,1]), p[[n]]$data, dir=rep("GO:0050807", nrow(p[[n]]$data)), dataset=rep(names(gsea_GO)[n],nrow(p[[n]]$data)))
}
df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}


GO_mat<-matrix(0,nrow=length(unique(unlist(genes_go))), ncol=length(genes_cor))
rownames(GO_mat)<-unique(unlist(genes_go))
for(i in 1:length(genes_cor)){
  GO_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go))))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go)))),]
}

df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])


pdf("Results/Figures/SYNAPSE_scatter.pdf",5,6)
ggplot(df, aes(x=aging, y=LSD))+geom_point()+geom_label_repel(label=df$gene)+geom_smooth(method="lm")+theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
dev.off()

