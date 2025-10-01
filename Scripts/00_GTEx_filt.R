rm(list=ls())
setwd("workdir")#workdir = working directory
source("Scripts/utilities.R")


library(data.table)
GTEx_PFC <- fread('Data/gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz',data.table=FALSE)
rownames(GTEx_PFC)<-GTEx_PFC[,2]

library("biomaRt")
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="https://useast.ensembl.org")
genes_conv<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                  values = gsub("\\..*","",rownames(GTEx_PFC)), 
                  mart = ensembl)

metadata<-read.csv('Data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', sep="\t")
symbol<-genes_conv[match(gsub("\\..*","",rownames(GTEx_PFC)), genes_conv[,1]),2]


GTEx_PFC<-changenames(GTEx_PFC[,-c(1:3)], anno=cbind(rownames(GTEx_PFC), symbol))

GTEx_PFC<-t(t(GTEx_PFC)/colSums(GTEx_PFC))*1000000
GTEx_PFC<-log2(GTEx_PFC+1)

#####metadata
#TODO mostrare che la distribuzione delle corr è skewed verso valori negativi,
#che i geni MT sono altamente espressi (fino al 65% dei reads)
#che i geni MT correlano con l'età
#che i geni MT sono alti dove i geni negativamente correlati sono bassi (esempio NMU)
#che i geni negativamente correlati con l'età sono espressi a bassi livelli dove sono alti gli MT
#ragionevole ipotizzare un bias dovuto alla normalizzazione, su cui influiscono moltissimo i valori molto alti degli MT
#eliminando i campioni con altissimi MT si riduce il problema iniziale

subj <- fread('https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
samp <- fread('https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')

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
