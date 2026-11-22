rm(list=ls())
setwd("workdir")#workdir = working directory
# CMap (LINCS phase II, GSE70138) input of 08_2_random_drugs.R: mean level-4 z-score of every
# compound in NPC cells. Reads the local gctx file with rhdf5 (no clue.io key, no internet;
# the original version used slinky, now removed from Bioconductor). Writes the same objects
# as the original version: compounds_NPC.RData (compounds) and CMAP_NPC_mean.RData (data_mean_all).
library(rhdf5)

gctx <- "Data/raw_data_notclosed/GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n345976x12328.gctx"
info <- "Data/raw_data_notclosed/GSE70138_Broad_LINCS_inst_info.txt"

##load info about perturbagens
sig<-read.csv(info, sep="\t")
gene_info<-read.csv("Data/raw_data_notclosed/GSE92742_Broad_LINCS_gene_info.txt", sep="\t")

compounds<-unique(sig$pert_iname[which(sig$cell_id=="NPC")])
save(compounds, file="Results/RData/compounds_NPC.RData")

# column (instance) and row (gene) ids of the gctx matrix (genes x instances)
col_ids <- h5read(gctx, "/0/META/COL/id")
row_ids <- trimws(h5read(gctx, "/0/META/ROW/id"))

data_mean_all<-matrix(nrow=length(row_ids), ncol=length(compounds))
for(comp in 1:length(compounds)){
col.ix <- sig[which(sig$cell_id=="NPC" & sig$pert_iname==compounds[comp]),1]
data <- h5read(gctx, "/0/DATA/0/matrix", index = list(NULL, which(col_ids %in% col.ix)))
data_mean<-rowMeans(data)
names(data_mean)<-gene_info$pr_gene_symbol[match(row_ids, gene_info$pr_gene_id)]
data_mean_all[,comp]<-data_mean
}
rownames(data_mean_all)<-names(data_mean)
save(data_mean_all, file="Results/RData/CMAP_NPC_mean.RData")
