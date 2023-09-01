load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Psychedelics metanalysis/Psychedelics_alldata.RData")

library(DESeq2)
library(openxlsx)

###normalization
#GSE179379<-GSE179379[which(rowSums(GSE179379>=5)>=5),]
dds <- DESeqDataSetFromMatrix(countData = GSE179379,
                              colData = data.frame(treatment=rep(c("LSD", "Saline"), each=10)),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE179379 <- results(dds, contrast = c("treatment","LSD", "Saline"))

save(DE_GSE179379, file="~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD/RData/DE_GSE179379.RData")
