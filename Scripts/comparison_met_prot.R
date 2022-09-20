library(openxlsx)
met<-read.xlsx("../Supplementary Table 1- DM CpG Sites.xlsx",1)
length(unique(met$gene.symbol[met$CpG.DM.FDR<0.05]))

met_up<-unique(met_promoter$gene.symbol[as.numeric(met$promoter.overlaps.mean_pvalue)<0.05 & met$CpG.methylation.difference>0  & met$feature.overlap.type=="promoter"])
met_dn<-unique(met_promoter$gene.symbol[as.numeric(met$promoter.overlaps.mean_pvalue)<0.05 & met$CpG.methylation.difference<0 & met$feature.overlap.type=="promoter"])
met_up<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference>0])
met_dn<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference<0])
met_up<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference>0 & met$feature.overlap.type=="promoter"])
met_dn<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference<0 & met$feature.overlap.type=="promoter"])
met_up<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference>0 & met$feature.overlap.type=="transcript"])
met_dn<-unique(met_promoter$gene.symbol[met$CpG.DM.FDR<0.05 & met$CpG.methylation.difference<0 & met$feature.overlap.type=="transcript"])

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]


prot_up<-read.xlsx("../Supplementary Table 2- Proteins With p Minor 0.05.xlsx",1)
prot_up<-prot_up[,3]

prot_dn<-read.xlsx("../Supplementary Table 2- Proteins With p Minor 0.05.xlsx",2)
prot_dn<-prot_dn[,3]

DE_GSE179379$methylation<-NA
DE_GSE179379$methylation[rownames(DE_GSE179379) %in% met_up]<-"up"
DE_GSE179379$methylation[rownames(DE_GSE179379) %in% met_dn]<-"dn"


DE_GSE179379$prot<-NA
DE_GSE179379$prot[rownames(DE_GSE179379) %in% prot_up]<-"up"
DE_GSE179379$prot[rownames(DE_GSE179379) %in% prot_dn]<-"dn"

DE_GSE179379<-data.frame(DE_GSE179379)
DE_GSE179379<-DE_GSE179379[which(DE_GSE179379$pvalue<0.05),]
ggplot(DE_GSE179379, aes(y=-log10(padj), x=log2FoldChange, colour=methylation))+geom_point()
ggplot(DE_GSE179379, aes(y=log2FoldChange, x=methylation))+geom_boxplot()

t.test(DE_GSE179379$log2FoldChange[DE_GSE179379$methylation=="up" & DE_GSE179379$pvalue<0.05], DE_GSE179379$log2FoldChange[is.na(DE_GSE179379$methylation) & DE_GSE179379$pvalue<0.05])
t.test(DE_GSE179379$log2FoldChange[DE_GSE179379$methylation=="up" & DE_GSE179379$pvalue<0.05], DE_GSE179379$log2FoldChange[DE_GSE179379$methylation=="dn" & DE_GSE179379$pvalue<0.05])
t.test(DE_GSE179379$log2FoldChange[DE_GSE179379$prot=="up"], DE_GSE179379$log2FoldChange[DE_GSE179379$prot=="dn"])
