library(slinky)

# update following lines with your details:
user_key <- "7212be05eb7400944e28addf97deeeb7"
gctx <- "~/Work/Projects/GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n345976x12328.gctx"
info <- "~/Work/Projects/GSE70138_Broad_LINCS_inst_info.txt"
sl <- Slinky(user_key, gctx, info)


##load info about perturbagens
sig<-read.csv("~/Work/Projects/GSE70138_Broad_LINCS_inst_info.txt", sep="\t")

compounds<-unique(sig$pert_iname[which(sig$cell_id=="NPC")])
save(compounds, file="~/Work/Projects/Metanalyses/Aging/compounds_NPC.RData")
col.ix <- sig[which(sig$cell_id=="NPC" & sig$pert_iname==compounds[1]),1]
data <- readGCTX(sl[, which(colnames(sl) %in% col.ix)])

data_mean_all<-matrix(nrow=nrow(data), ncol=length(compounds))
for(comp in 1:length(compounds)){
col.ix <- sig[which(sig$cell_id=="NPC" & sig$pert_iname==compounds[comp]),1]


data <- readGCTX(sl[, which(colnames(sl) %in% col.ix)])
gene_info<-read.csv("~/Work/Projects/GSE92742_Broad_LINCS_gene_info.txt", sep="\t")
data_mean<-rowMeans(data)
names(data_mean)<-gene_info$pr_gene_symbol[match(names(data_mean), gene_info$pr_gene_id)]
data_mean_all[,comp]<-data_mean
}
rownames(data_mean_all)<-names(data_mean)
save(data_mean_all, file="CMAP_NPC_mean.RData")
sort(data_mean)[1:10]
sort(data_mean, decreasing=T)[1:10]

library(cmapR)
my_ds <- parse_gctx("~/Work/Projects/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
