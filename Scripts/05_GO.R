rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")

load("Results/RData/DE_GSE179379.RData")
load("Results/RData/Alldata_20Sep.RData")
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
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier
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
dat<-setdiff(dat, c("GSE5388", "GSE102741"))#outlier

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

allNES_mat_sel<-allNES_mat[which(rownames(allNES_mat) %in% DE_GO_LSD),]
allp_mat_sel<-allp_mat[which(rownames(allp_mat) %in% DE_GO_LSD),]
pup<-list()
ind<-0
for(up in ego_up$ID[which(ego_up$ID %in% rownames(allNES_mat_sel))]){
  ind<-ind+1
  allNES_mat_sel[up, ]
  
  istwo <- rep(T, length(gsea_GO))
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
  
  istwo <- rep(T, length(gsea_GO))
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

dev.off()
pdf("Results/Figures/GO_shared_dn.pdf",8,5)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T, annotation_row = anno_p, cluster_rows = F)
dev.off()

path_sel<-c(names(pup[which(pup<10^(-20))]))
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

DE_GSE179379_up<-DE_GSE179379
rownames(DE_GSE179379_up)<-toupper(rownames(DE_GSE179379))

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

df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

library(ggrepel)
pdf("Results/Figures/SYNAPSE_scatter.pdf",5,6)
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

GO_mat<-matrix(0,nrow=length(unique(unlist(genes_go))), ncol=length(genes_cor))
rownames(GO_mat)<-unique(unlist(genes_go))
for(i in 1:length(genes_cor)){
  GO_mat[rownames(genes_cor[[i]])[which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go))))],i]<-genes_cor[[i]][which(rownames(genes_cor[[i]]) %in% unique(unique(unlist(genes_go)))),]
}

df<-data.frame(gene=rownames(GO_mat), aging=rowMeans(GO_mat), LSD=DE_GSE179379_up[rownames(GO_mat),"log2FoldChange"])

library(ggrepel)
pdf("Results/Figures/PROTEASOME_scatter.pdf",5,6)
ggplot(df, aes(x=aging, y=LSD))+geom_point()+geom_label_repel(label=df$gene)+geom_smooth(method="lm")+theme_classic()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
dev.off()
