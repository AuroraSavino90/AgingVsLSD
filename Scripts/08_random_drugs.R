rm(list=ls())
setwd("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD")

load("Results/RData/DE_GSE179379.RData")
load("Results/RData/Alldata_20Sep.RData")
homologs<-read.csv("Data/Human rat homologs.txt")
#######################
### random drugs
######################

load("Results/RData/CMAP_NPC_mean.RData")
load("Results/RData/compounds_NPC.RData")

fgsea_res_drugs<-list()
for(comp in 1:ncol(data_mean_all)){
  print(comp)
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE5388", "GSE102741"))
  
  genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
  genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  genes_up<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=T)[1:length(genes_up)]]
  genes_dn<-rownames(data_mean_all)[order(data_mean_all[,comp], decreasing=F)[1:length(genes_dn)]]
  
  dat_names<-c()
  genes_cor<-list()
  n<-0
  for(dd in dat){
    n<-n+1
    data<-get(dd)
    sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    age<-as.numeric(age)
    if(length(which(age>=7300))>=10){
      sample<-sample[which(age>=7300)]
      data<-data[, sample]
      age<-age[which(age>=7300)]
      
      genes_cor[[n]]<-cor(t(data), as.numeric(age))
      dat_names<-c(dat_names, dd)
    }
  }
  
  region<-"PFC"
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE5388", "GSE102741"))
  
  for(dd in dat){
    n<-n+1
    data<-get(dd)
    sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    data<-data[, sample]
    age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
    age<-as.numeric(age)
    if(length(which(age>=7300))>=10){
      sample<-sample[which(age>=7300)]
      data<-data[, sample]
      age<-age[which(age>=7300)]
      
      genes_cor[[n]]<-cor(t(data), as.numeric(age))
      dat_names<-c(dat_names, dd)
    }
  }
  #LSD
  #rat_homologs<-unique(homologs[homologs[,2] %in% rownames(GSE179379)[which(rowSums(GSE179379>10)>1)],3])
  
  genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
  fgsea_res_drugs[[comp]]<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){
    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    fgsea_res_drugs[[comp]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
    
  }
  
}
save(fgsea_res_drugs, file="fgsea_res_drugs_RAW.RData")

load("Results/RData/fgsea_res.RData")

fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})


pvals_dn<-unlist(lapply(fgsea_res, function(x){x[1,3]}))
pvals_dn_rand<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,3]*ifelse(sign(x[1,6])==(-1), NA, 1)}))})
pvals_dn_rand_mat<-matrix(nrow=length(pvals_dn_rand), unlist(pvals_dn_rand), byrow=T)

df<-data.frame(p=unlist(pvals_dn_rand), dataset=rep(dat_names, each=length(pvals_dn_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_dn, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_p_dn_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

better_dn<-lapply(pvals_dn_rand, function(x){x<pvals_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
(colSums(better_dn_mat, na.rm=T)/193)[order(unlist(fgsea_rank), decreasing=T)]
compounds[order(rowSums(better_dn_mat, na.rm=T), decreasing = T)][1:10]

pvals_up<-unlist(lapply(fgsea_res, function(x){x[2,3]}))
pvals_up_rand<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,3]*ifelse(sign(x[2,6])==(1), NA, 1)}))})

df<-data.frame(p=unlist(pvals_up_rand), dataset=rep(dat_names, each=length(pvals_up_rand)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(p=pvals_up, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_p_up_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=(-log10(p))))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=(-log10(p))), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_classic()
dev.off()

better_up<-lapply(pvals_up_rand, function(x){x<pvals_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat, na.rm=T)/193
compounds[order(rowSums(better_up_mat, na.rm=T), decreasing = T)][1:10]

pdf("Results/Figures/Random_drugs_p_perc_RAW.pdf", 4,4)
plot(1-colSums(better_dn_mat, na.rm=T)/193, 1-colSums(better_up_mat, na.rm=T)/193, xlab="% of tests less significant than LSD (dn)",
     ylab="% of tests less significant than LSD (up)", pch=19)
dev.off()

################################
###evaluation with nes
##############################
nes_dn<-unlist(lapply(fgsea_res, function(x){x[1,6]}))
nes_dn_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[1,6]}))})

df<-data.frame(NES=unlist(nes_dn_drug), dataset=rep(dat_names, each=length(nes_dn_drug)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(NES=nes_dn, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_nes_dn_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=NES))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=NES), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_bw()
dev.off()

better_dn<-lapply(nes_dn_drug, function(x){x>nes_dn})
better_dn_mat<-matrix(nrow=length(better_dn), unlist(better_dn), byrow=T)
##number of random lists performing better in each dataset
colSums(better_dn_mat)/193

pdf("Results/Figures/Random_drugs_nes_dn_RAW_scatter.pdf", 4,4)
plot(-log10(pvals_dn), 1-colSums(better_dn_mat)/193, pch=19, ylab="% of tests with NES higher than LSD (dn)")
dev.off()

nes_up<-unlist(lapply(fgsea_res, function(x){x[2,6]}))
nes_up_drug<-lapply(fgsea_res_drugs, function(y){unlist(lapply(y, function(x){x[2,6]}))})

df<-data.frame(NES=unlist(nes_up_drug), dataset=rep(dat_names, each=length(nes_up_drug)))
df$type<-rep("random", nrow(df))
df<-rbind.data.frame(df, data.frame(NES=nes_up, dataset=dat_names, type="LSD"))
df$dataset<-factor(df$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

pdf("Results/Figures/Random_drugs_nes_up_RAW.pdf", 8,8)
ggplot(subset(df, type="random"), aes(x=NES))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=NES), colour="red", linetype="dashed")+facet_wrap(.~dataset)+theme_bw()
dev.off()

####
better_up<-lapply(nes_up_drug, function(x){x<nes_up})
better_up_mat<-matrix(nrow=length(better_up), unlist(better_up), byrow=T)
##number of random lists performing better in each dataset
colSums(better_up_mat)/193
pdf("Results/Figures/Random_drugs_nes_up_RAW_scatter.pdf", 4,4)
plot(-log10(pvals_up), 1-colSums(better_up_mat)/193, pch=19, ylab="% of tests with NES lower than LSD (up)")
dev.off()

pdf("Results/Figures/Random_drugs_nes_perc_RAW.pdf", 4,4)
plot(1-colSums(better_dn_mat, na.rm=T)/193, 1-colSums(better_up_mat, na.rm=T)/193, xlab="% of tests with higher NES than LSD (dn)",
     ylab="% of tests with NES lower than LSD (up)", pch=19)
dev.off()

plot(colSums(better_dn_mat, na.rm=T)/193, colSums(better_up_mat, na.rm=T)/193, xlab="% of tests with higher NES than LSD (dn)",
     ylab="% of tests with NES lower than LSD (up)", pch=19)

p_dn<-sumlog(colSums(better_dn_mat, na.rm=T)/193)
p_up<-sumlog(colSums(better_up_mat, na.rm=T)/193)

######evaluation collapsing the pvalue

library(metap)

meta_dn<-function(x){
  istwo <- rep(T, length(x))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up<-function(x){
  istwo <- rep(T, length(x))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

pdn<-meta_dn(fgsea_res)
sump_dn<-lapply(fgsea_res_drugs, meta_dn)
sum(unlist(sump_dn)<2.2*10^(-16))
sum(unlist(sump_dn)<pdn)
compounds[which(unlist(sump_dn)<pdn)]

df<-data.frame(p=unlist(sump_dn), type=rep("random", length(sump_dn)))
df<-rbind.data.frame(df, data.frame(p=pdn, type="LSD"))

pdf("Results/Figures/Random_drugs_fisher_dn_RAW.pdf", 4,4)
ggplot(subset(df, type="random"), aes(x=-log10(p)))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=-log10(p)), colour="red", linetype="dashed")+theme_classic()
dev.off()

pup<-meta_up(fgsea_res)
sump_up<-lapply(fgsea_res_drugs, meta_up)
sum(unlist(sump_up)<2.2*10^(-16))
sum(unlist(sump_up)<pup)

df<-data.frame(p=unlist(sump_up), type=rep("random", length(sump_up)))
df<-rbind.data.frame(df, data.frame(p=pup, type="LSD"))

pdf("Results/Figures/Random_drugs_fisher_up_RAW.pdf", 4,4)
ggplot(subset(df, type="random"), aes(x=-log10(p)))+geom_density()+
  geom_vline(data=subset(df, type=="LSD"), aes(xintercept=-log10(p)), colour="red", linetype="dashed")+theme_classic()
dev.off()


