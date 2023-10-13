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

symbol<-genes_conv[match(gsub("\\..*","",rownames(GTEx_PFC)), genes_conv[,1]),2]


GTEx_PFC<-changenames(GTEx_PFC[,-c(1:3)], anno=cbind(rownames(GTEx_PFC), symbol))

GTEx_PFC<-t(t(GTEx_PFC)/colSums(GTEx_PFC))*1000000
GTEx_PFC<-log2(GTEx_PFC+1)

MT_genes<-colSums(2^GTEx_PFC[grep("MT-", rownames(GTEx_PFC)),]-1)
GTEx_PFC<-GTEx_PFC[,which(MT_genes<quantile(MT_genes,0.9))]

GTEx_PFC<-GTEx_PFC[rowSums(GTEx_PFC>0)>=10,]

save(GTEx_PFC, file="Results/RData/GTEx_RPM_PFC_filt.RData")
