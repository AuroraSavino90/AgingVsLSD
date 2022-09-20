rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD")

library(ggplot2)
load("Results/RData/DE_GSE179379.RData")
load("Results/RData/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
####################
### GO categories
###################
forgesea<-sort(forgesea, decreasing=T)
gsea_GO<-gseGO(forgesea,
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
