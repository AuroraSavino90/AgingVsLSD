rm(list=ls())
setwd("workdir")#workdir = working directory

library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
load("Data/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
####################
### GO categories
###################
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

pdf("Results/Figures/DEG_up.pdf", 6,5)
dotplot(ego_up, showCategory=10)
dev.off()
pdf("Results/Figures/DEG_dn.pdf", 6,5)
dotplot(ego_down, showCategory=10)
dev.off()
