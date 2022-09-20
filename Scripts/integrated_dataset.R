rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")

library(ggplot2)
load("RData/DE_GSE179379.RData")
load("RData/Alldata_27May.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, "GSE11512")
dat<-setdiff(dat, c("GSE30272"))#no raw
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
dat<-setdiff(dat, c("GSE5388"))#batch PC1>40% of variance

dat_names<-c()
age_all<-c()
dataset_all<-c()
for(dd in dat){

  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]
   dat_names<-c(dat_names, dd)
   age_all<-c(age_all, age)
   dataset_all<-c(dataset_all, rep(dd, ncol(data)))
  }
}

all<-c()
for(dd in dat_names){
  data<-get(dd)
  all<-c(all, rownames(data))
  #print(length(unione))
}
sum(table(all)>=9)
selgenes<-names(which(table(all)>=9))

############################ALLDATA
data<-get(dat_names[1])
sample<-metadata$Sample[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
data<-data[, sample]
age<-metadata$Age[which(metadata$Dataset==dat_names[1] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
age<-as.numeric(age)
if(length(which(age>6570))>=10){
  sample<-sample[which(age>6570)]
  data<-data[, sample]
alldata<-data[selgenes,]
}

for(i in 2:length(dat_names)){
  data<-get(dat_names[i])
  sample<-metadata$Sample[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dat_names[i] & metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]


  alldata<-cbind(alldata, data[selgenes,])
  }
}

library(FactoMineR)
pca<-PCA(t(alldata))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], age=age_all/365, dataset=dataset_all)

pdf("Results/Figures/PCA_DLPF_dataset.pdf",7,5)
ggplot(df, aes(PC1, PC2, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_DLPF_age.pdf",6,5)
ggplot(df, aes(PC1, PC2, colour=age))+geom_point(size=1)+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()
pdf("Results/Figures/PC1_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC1, PC4, colour=dataset))+geom_point()+theme_classic()
ggplot(df, aes(PC1, PC4, colour=age))+geom_point()+theme_classic()


batch = dataset_all

# parametric adjustment
library(sva)
mod = model.matrix(~age_all, data=data.frame(alldata))
combat_edata1 = ComBat(dat=alldata, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
combat_edata1 = ComBat(dat=alldata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

pca<-PCA(t(combat_edata1))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC1, PC2, colour=dataset))+geom_point()+theme_classic()
ggplot(df, aes(PC1, PC2, colour=age))+geom_point()+theme_classic()

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC1, PC4, colour=dataset))+geom_point()
ggplot(df, aes(PC2, PC4, colour=age))+geom_point()

mod0 = model.matrix(~1,data=data.frame(t(alldata)))
n.sv = num.sv(data.frame(alldata),mod,method="leek")
svobj = sva(as.matrix(alldata),mod,mod0,n.sv=8)
plot(age_all, svobj[[1]][,1], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,2], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,3], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,4], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,5], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,6], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,7], col=as.factor(dataset_all))
plot(age_all, svobj[[1]][,8], col=as.factor(dataset_all))

alldata_resid<-apply(alldata,1, function(x){residuals(lm(unlist(x)~svobj[[1]][,1]+svobj[[1]][,2]+svobj[[1]][,3]+svobj[[1]][,4]+svobj[[1]][,5]+svobj[[1]][,6]+svobj[[1]][,7]+svobj[[1]][,8]))})
alldata_resid<-t(alldata_resid)

pca<-PCA(t(alldata_resid))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
pdf("Results/Figures/PCA_SVA_dataset.pdf",7,5)
ggplot(df, aes(PC1, PC2, colour=dataset))+geom_point()+theme_classic()
dev.off()
pdf("Results/Figures/PCA_SVA_age.pdf",6,5)
ggplot(df, aes(PC1, PC2, colour=age))+geom_point()+ scale_color_gradient(low = "orange", high = "blue")+theme_classic()
dev.off()

ggplot(df, aes(PC3, PC4, colour=age))+geom_point()

pdf("Results/Figures/PC1_SVA_age.pdf",7,5)
ggplot(df, aes(age, PC1, colour=dataset))+geom_point()+theme_classic()
dev.off()

save(alldata_resid, file="Results/RData/alldata_resid.RData")

#################################################
#########GSEA with the batch corrected dataset
##################################################

library(fgsea)
library(ggplot2)

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

genes_cor<-cor(t(alldata_resid), as.numeric(age_all))

rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DE_GSE179379),3])

  forgesea<-unlist(genes_cor)
  names(forgesea)<-rownames(genes_cor)
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res<-fgsea(list(UP=genes_up, DN=genes_dn), sort(forgesea))
  
  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)
  
  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data))), 
                       data.frame(p2$data, dir=rep("dn", nrow(p2$data))))
  
  p<-ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

pdf("Results/Figures/GSEA_all_batchcorrect.pdf", 5,4)
ggplot(df, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))
dev.off()

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

##################################
##### showing paths in common
##################################

DE_GO_LSD<-union(ego_up$ID, ego_down$ID)

allpaths_mat<-data.frame(gsea_GO)
allpaths_mat_up<-allpaths_mat[which(rownames(allpaths_mat) %in% ego_up$ID),]
allpaths_mat_up$p.adjust<-p.adjust(allpaths_mat_up$pvalue, method="fdr")
allpaths_mat_dn<-allpaths_mat[which(rownames(allpaths_mat) %in% ego_down$ID),]
allpaths_mat_dn$p.adjust<-p.adjust(allpaths_mat_dn$pvalue, method="fdr")


############
#### TF
#############

encode <- read.gmt("Data/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

fgsea_TF<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=encode, pvalueCutoff = 1, maxGSSize = 10000)


DEG_down_conv<-unique(homologs[which(homologs[,2] %in% DEG_down),3])
DEG_up_conv<-unique(homologs[which(homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down<- enricher(gene = DEG_down_conv,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   universe=universe, 
                   TERM2GENE = encode)

TF_up<- enricher(gene = DEG_up_conv,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 universe=universe, TERM2GENE = encode)

allpaths_mat<-data.frame(fgsea_TF)
allpaths_mat_up<-allpaths_mat[which(rownames(allpaths_mat) %in% TF_up$ID),]
allpaths_mat_up$p.adjust<-p.adjust(allpaths_mat_up$pvalue, method="fdr")
allpaths_mat_dn<-allpaths_mat[which(rownames(allpaths_mat) %in% TF_down$ID),]
allpaths_mat_dn$p.adjust<-p.adjust(allpaths_mat_dn$pvalue, method="fdr")

#######Top correlated genes
sort(genes_cor[which(rownames(genes_cor) %in% genes_up),])[1:10]
sort(genes_cor[which(rownames(genes_cor) %in% genes_dn),], decreasing = T)[1:10]


cor_age<-function(organism, diagnosis, region){
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))

  #compare aging datasets
  data<-get(dat[1])
  unione<-rownames(data)
  for(dd in dat){
    data<-get(dd)
    unione<-union(unione, rownames(data))
    print(length(unione))
  }

  cor_age<-matrix(nrow=length(unione), ncol=length(dat))
  colnames(cor_age)<-dat
  rownames(cor_age)<-unione
  for(dd in dat){
    data<-get(dd)
    genes<-intersect(rownames(data), unione)
    sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl%in% region)]
    data<-data[, sample]
    age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd &metadata$Organism==organism & metadata$Diagnosis==diagnosis & metadata$Region_simpl%in% region)]))
    if(length(which(age>6570))>=10){
      sample<-sample[which(age>6570)]
      data<-data[, sample]
      age<-age[which(age>6570)]
      cor_age[genes,dd]<-cor(t(data[genes,]), age, method="p", use="pairwise.complete.obs")
    }
  }
  cor_age<-cor_age[,colSums(is.na(cor_age))!=nrow(cor_age)]
  return(cor_age)

}
cor_DLPFC<-cor_age(organism="Homo sapiens", region="DLPFC", diagnosis="Healthy")
goodgenes<-names(sort(abs(rowSums(cor_DLPFC, na.rm=T))[1:1000]))
#goodgenes<-names(sort(abs(rowSums((cor_DLPFC), na.rm=T))[1:100]))


pca<-PCA(t(alldata[rownames(alldata) %in% goodgenes,]))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2],PC3=pca$ind$coord[,3],PC4=pca$ind$coord[,4], age=age_all/365, dataset=dataset_all)
ggplot(df, aes(PC3, PC4, colour=dataset))+geom_point()
ggplot(df, aes(PC3, PC4, colour=age))+geom_point()

