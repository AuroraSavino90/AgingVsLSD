rm(list=ls())
setwd("workdir")#workdir = working directory
source("Scripts/utilities.R")


library(data.table)
GTEx_PFC <- fread('Data/gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz',data.table=FALSE)
rownames(GTEx_PFC)<-GTEx_PFC[,2]

# Ensembl 115 instead of live biomaRt: the published Fig S16 was made on 1 Oct 2025 with the live
# Ensembl site, then release 115 (with 115 the histograms of S16A/E are identical to the published ones).
# This script only makes Fig S16: the GTEx_PFC used in the analyses (Alldata_20Sep) is built in
# Scripts/data_preparation with Ensembl 107, as in 2022.
source("Scripts/data_preparation/ensembl_release_map.R")
genes_conv <- ensembl_release_map(115, "hsapiens", "Data/ensembl")

metadata<-read.csv('Data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', sep="\t")
symbol<-genes_conv[match(gsub("\\..*","",rownames(GTEx_PFC)), genes_conv[,1]),2]


GTEx_PFC<-changenames(GTEx_PFC[,-c(1:3)], anno=cbind(rownames(GTEx_PFC), symbol))

GTEx_PFC<-t(t(GTEx_PFC)/colSums(GTEx_PFC))*1000000
GTEx_PFC<-log2(GTEx_PFC+1)

#####metadata

gtex_url <- 'https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/'  # local copies in Data/ are read when present
subj <- fread(if (file.exists('Data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')) 'Data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt' else paste0(gtex_url, 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'))
samp <- fread(if (file.exists('Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')) 'Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt' else paste0(gtex_url, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))

samp[, SUBJID := gsub('([^-]*)-([^-]*)-.*', '\\1-\\2', SAMPID)]

PFC_subj<-samp$SUBJID[which(samp$SAMPID %in% colnames(GTEx_PFC))]
PFC_samp<-samp$SAMPID[which(samp$SAMPID %in% colnames(GTEx_PFC))]

age<-subj$AGE[match(PFC_subj, subj$SUBJID)]
age_bin<-age
age_bin[age_bin=="20-29"]<-1
age_bin[age_bin=="30-39"]<-2
age_bin[age_bin=="40-49"]<-3
age_bin[age_bin=="50-59"]<-4
age_bin[age_bin=="60-69"]<-5
age_bin[age_bin=="70-79"]<-6

MT_genes<-colSums(2^GTEx_PFC[grep("MT-", rownames(GTEx_PFC)),]-1)

#FIG2
pdf("Results/Figures/GTExB.pdf")
hist(MT_genes/1000000)
dev.off()

#FIG4
boxplot(MT_genes~age_bin)
plot(MT_genes~age_bin)
cor.test(MT_genes,as.numeric(age_bin))

cor_age<-cor(t(GTEx_PFC), as.numeric(age_bin), method="p", use="pairwise.complete.obs")
cor_MT_genes<-cor(t(GTEx_PFC), log2(as.numeric(MT_genes)+1), method="p", use="pairwise.complete.obs")

#FIG 1
pdf("Results/Figures/GTExA.pdf")
hist(cor_age)
abline(v=0, col="red",  lty=c(2), lwd=c(3))
dev.off()

rownames(GTEx_PFC)[order(cor_age, decreasing = T)[1:10]]

#FIG3
pdf("Results/Figures/GTExC.pdf",5,5)
plot(GTEx_PFC["NMU",], log2(MT_genes+1), pch=19)
dev.off()
cor.test(GTEx_PFC["NMU",], log2(MT_genes+1))

plot(age_bin, GTEx_PFC["NMU",])
cor.test(as.numeric(age_bin), log2(MT_genes+1), method="p", use="pairwise.complete.obs")
#FIG5
pdf("Results/Figures/GTExE.pdf")
plot(cor_age, cor_MT_genes, pch=19)
dev.off()
cor.test(cor_age, cor_MT_genes, method="p", use="pairwise.complete.obs")

cor_NMU<-cor(t(GTEx_PFC), GTEx_PFC["NMU",], method="p", use="pairwise.complete.obs")
rownames(GTEx_PFC)[order(cor_NMU, decreasing = F)[1:10]]

expr_genes<-names(sort(rowSums(2^GTEx_PFC-1), decreasing=T)[1:10])
apply(2^(GTEx_PFC[expr_genes,])-1,1, var)
vars<-apply(2^(GTEx_PFC)-1,1, var)
hist(vars)

GTEx_PFC<-GTEx_PFC[,which(MT_genes<quantile(MT_genes,0.9))]

cor_age_filt<-cor(t(GTEx_PFC), as.numeric(age_bin)[which(MT_genes<quantile(MT_genes,0.9))], method="p", use="pairwise.complete.obs")

#FIG6
pdf("Results/Figures/GTExF.pdf")
hist(cor_age_filt)
abline(v=0, col="red",  lty=c(2), lwd=c(3))
dev.off()

GTEx_PFC<-GTEx_PFC[rowSums(GTEx_PFC>0)>=10,]

save(GTEx_PFC, file="Results/RData/GTEx_RPM_PFC_filt.RData")
