rm(list=ls())
setwd("workdir")#workdir = working directory
load("Data/Psychedelics_PFC.RData")

library(DESeq2)

###normalization
#GSE179379<-GSE179379[which(rowSums(GSE179379>=5)>=5),]
dds <- DESeqDataSetFromMatrix(countData = GSE179379,
                              colData = data.frame(treatment=rep(c("LSD", "Saline"), each=10)),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE179379 <- results(dds, contrast = c("treatment","LSD", "Saline"))

save(DE_GSE179379, file="Results/RData/DE_GSE179379.RData")
