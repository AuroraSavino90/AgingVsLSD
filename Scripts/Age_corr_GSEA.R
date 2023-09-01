library(fgsea)
library(ggplot2)
library(metap)

load("Data/Alldata_20Sep.RData")
load("Data/DE_GSE179379.RData")
rat_homologs<-read.csv("Data/Human rat homologs.txt")
source("functions.R")

######Define example lists
genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(rat_homologs[which(rat_homologs[,2] %in% genes_up),3])
genes_dn<-unique(rat_homologs[which(rat_homologs[,2] %in% genes_dn),3])

######data preprocessing
metadata$Age<-as.numeric(metadata$Age)
#TODO
#fix sample names for GSE101521
#add two validation datasets
#fix smoke, medication, exted region_simpl and remove region

######input data (FROM THE USER)
user_data<-list(Sample=NA, Dataset=NA, Age=c(20,106), Gender=NA, Race=NA, 
                Suicide=NA, Diagnosis="Healthy", Alcohol=NA, Drugs=NA,
                Smoke=NA, Medication=NA, Organism="Homo sapiens",
                Region=NA, Platform=NA, Region_simpl="DLPFC")
input_genes<-c("BDNF", "PTEN", "SPARC", "ARC", "FOXC1", "FOXM1")
#convert age to days
user_data[["Age"]]<-user_data[["Age"]]*365

#select datasets and samples
dat_and_samples<-dataset_select(user_data, metadata, 10) 
datasets<-dat_and_samples[[1]]
samples<-dat_and_samples[[2]]

##compute correlation with age (NOT NECESSARY FOR GSEA)
corrs_all<-lapply(datasets, function(x){corr_age(dataset=x, samples=samples, input_genes=genes_up)})
names(corrs_all)<-datasets

##GSEA
#TODO make it possibe to do for multiple lists in parallel
#TODO make it possible to use GO categories
GSEA_all<-lapply(datasets, function(x){GSEA_one(dataset=x, samples=samples, input_genes=genes_up)})
names(GSEA_all)<-datasets

##compute overall pvalue
meta_age(datasets, GSEA_all, direction="reversal")
meta_age(datasets, GSEA_all, direction="mimick")

###combined plot
plots<-lapply(GSEA_all, function(x){ x[[1]] })
library(gridExtra)
n <- length(plots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plots, ncol=nCol))



