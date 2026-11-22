# ============================================================
# THRESHOLD-SENSITIVITY COPY of Scripts/05_DEvsAging_validation.R
# Changes vs the original: (1) in vivo LSD genes selected with the column
# set by SELECTION_MODE (P_SEQ: pvalue or padj); (2) outputs written to OUT
# (Threshold_sensitivity_invivo/Results_<mode>); (3) see invivo_setup.R.
# Run via Threshold_sensitivity_invivo/run_all.R
# ============================================================
rm(list=ls())
source("Threshold_sensitivity_invivo/Scripts/invivo_setup.R")  # SELECTION_MODE, OUT, P_SEQ, shims (replaces setwd)
# Published results: Bioconductor 3.21 annotations (org.Mm.eg.db, GO.db 3.21.0) and Ensembl 115.
set.seed(1234)  # fixed seed for fgsea/gseGO permutations (the published GSEA values come from unseeded runs; see REPRODUCIBILITY.md)
library(openxlsx)
library(ggplot2)
library(DESeq2)

source("Scripts/utilities.R")

# 15 samples deposited in GEO (GSE263662): cohorts C1-C3, one mouse per dose (0, 50, 100, 200, 500 ug/kg) per cohort
meta<-read.xlsx("Data/raw_data_notclosed/Invivo RNAseq/Samples_LSD.xlsx",1)
meta$dose[meta$dose=="ctrl"]<-0

#load counts (two lanes per sample, files ordered as the samples in meta)
filesL1<-list.files("Data/raw_data_notclosed/Invivo RNAseq", pattern="L001")
filesL2<-list.files("Data/raw_data_notclosed/Invivo RNAseq",pattern="L002")
filesL1<-filesL1[match(paste0("_S", meta$sample, "_"), regmatches(filesL1, regexpr("_S[0-9]+_", filesL1)))]
filesL2<-filesL2[match(paste0("_S", meta$sample, "_"), regmatches(filesL2, regexpr("_S[0-9]+_", filesL2)))]
L1mat<-read.csv(paste("Data/raw_data_notclosed/Invivo RNAseq", filesL1[1], sep="/"), sep="\t")
for(i in filesL1){
  L1mat<-cbind(L1mat, read.csv(paste("Data/raw_data_notclosed/Invivo RNAseq", i, sep="/"), sep="\t")[,2])
}
rownames(L1mat)<-L1mat[,1]
L1mat<-L1mat[,-c(1,2)]
L2mat<-read.csv(paste("Data/raw_data_notclosed/Invivo RNAseq", filesL2[1], sep="/"), sep="\t")
for(i in filesL2){
  L2mat<-cbind(L2mat, read.csv(paste("Data/raw_data_notclosed/Invivo RNAseq", i, sep="/"), sep="\t")[,2])
}
rownames(L2mat)<-L2mat[,1]
L2mat<-L2mat[,-c(1,2)]

counts<-L1mat+L2mat

#convert to gene symbol
# Ensembl 115 (the live release in Sep 2025) instead of live biomaRt
source("Scripts/data_preparation/ensembl_release_map.R")
ids <- ensembl_release_map(115, "mmusculus", "Data/ensembl")
ids <- ids[ids$ensembl_gene_id %in% rownames(counts), ]
counts_mgi<-changenames(data=counts, anno=ids)

#filter genes with at least 10 reads in at least 3 samples
counts_mgi_filt<-counts_mgi[rowSums(counts_mgi>=10)>=3,]

#remove genes with variance < 0.5 after normalisation as reads per million (RPM)
RPM<-t(t(counts_mgi_filt)/colSums(counts_mgi_filt))*1000000
vargenes<-apply(RPM,1, var)
counts_mgi_filt<-counts_mgi_filt[which(vargenes>0.5),]

#differential expression vs saline, adjusting for cohort
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_mgi_filt,
                                      colData = data.frame(meta),
                                      design= ~ dose+cohort)
dds <- DESeq2::DESeq(dds)
DE_500 <- DESeq2::results(dds, contrast = c("dose", "500", "0"))
DE_200 <- DESeq2::results(dds, contrast = c("dose","200","0"))
DE_100 <- DESeq2::results(dds, contrast = c("dose", "100", "0"))
DE_50 <- DESeq2::results(dds, contrast = c("dose", "50", "0"))


library(fgsea)
library(ggrepel)

load("Data/original_RData/Alldata_20Sep.RData")
assign(".Random.seed", readRDS("Data/original_RData_notclosed/rng_states.rds")$Alldata_20Sep, envir = .GlobalEnv)  # random state of the published runs (the .Random.seed saved in Data/original_RData/Alldata_20Sep.RData), set explicitly
rat_homologs<-read.csv("Data/homologs/Human rat homologs.txt", sep="\t")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022
mouse_homologs<-read.csv("Data/homologs/Human mouse homologs.txt", sep="\t")  # Ensembl 107: BioMart export, human genes (GRCh38.p13) with mouse homologs, August 2022

library(metap)

meta_dn_rev<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up_rev<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

meta_dn<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
    
  }
  
  return(p[[3]])
}

meta_up<-function(x){
  istwo <- rep(T, length(dat_names))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}


######################### RNA-seq
DEGs_seq<-c("DE_500", "DE_200", "DE_100", "DE_50")
cor_up_tot<-list()
cor_dn_tot<-list()
i<-0
n<-0
pdn<-list()
pup<-list()

for(DE in DEGs_seq){
  i<-i+1
  n<-n+1
homologs<-mouse_homologs
 
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs[[P_SEQ]]<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs[[P_SEQ]]<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  region<-c("DLPFC")
  diagnosis<-"Healthy"
  dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
  dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance
  
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
  dat<-setdiff(dat, c("GSE102741", "GSE5388"))#batch PC1>40% of variance
  
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
  
  genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
  fgsea_res<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){
    
    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)
    
    p1<-plotEnrichment(genes_up, forgesea)
    p2<-plotEnrichment(genes_dn, forgesea)
    
    df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))
    
    p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
      scale_colour_manual(values=c("dn"="blue", "up"="red"))
    
  }
  
  pdn[[i]]<-meta_dn_rev(fgsea_res)
  pup[[i]]<-meta_up_rev(fgsea_res)
  
  df_tot<-df[[1]]
  for(n in 2:length(genes_cor)){
    df_tot<-rbind.data.frame(df_tot, df[[n]])
  }
  
  
  fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})
  
  df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])
  
  pdf(paste(file.path(OUT, "Figures/GSEA_aging_"),DE,".pdf", sep=""), 8, 8)
  print(ggplot(df_tot, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
          scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))
  dev.off()
  
  save(fgsea_res, file=paste(file.path(OUT, "RData/fgsea_res_"), DE, ".RData", sep=""))
}



##################
## collapsing GSEA p (no random drugs)
##################
p_dn<-c()
p_up<-c()
for(DE in c(DEGs_seq)){
  load(paste(file.path(OUT, "RData/fgsea_res_"), DE, ".RData", sep=""))
  
  p_dn<-c(p_dn, meta_dn(fgsea_res))
  p_up<-c(p_up, meta_up(fgsea_res))
  
}
#p_dn[p_dn==0]<-10^(-200)

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_seq))
ggplot(df, aes(x=pup, y=pdn, label=dataset))+geom_point()+geom_label_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= -log10(0.05)), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= -log10(0.05)), colour="red", linetype="dashed")


p_dn_rev<-c()
p_up_rev<-c()
for(DE in c(DEGs_seq)){
  load(paste(file.path(OUT, "RData/fgsea_res_"), DE, ".RData", sep=""))
  
  p_dn_rev<-c(p_dn_rev, meta_dn_rev(fgsea_res))
  p_up_rev<-c(p_up_rev, meta_up_rev(fgsea_res))
  
}

df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c( DEGs_seq))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_seq)))


df<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=DEGs_seq)

pdf(file.path(OUT, "Figures/All_DE_aging.pdf"),5,5)
ggplot(df, aes(x=pup_diff, y=pdn_diff, label=dataset))+geom_point()+geom_text_repel()+
  theme_classic()+
  geom_vline(aes(xintercept= 0), colour="red", linetype="dashed")+
  geom_hline(aes(yintercept= 0), colour="red", linetype="dashed")+
  scale_colour_manual(values=c("Positive Control"="#8E24AA", "Negative Control"="#F57C00", "Psychoplastogen"="black"))+
  xlim(-50,50)+
  ylim(-250, 250)
dev.off()



library(data.table)
library(metap)

gsea_dt<-data.frame()
for(DE in DEGs_seq){
  load(file=paste(file.path(OUT, "RData/fgsea_res_"), DE, ".RData", sep=""))
  gsea_tmp <- rbindlist(
    lapply(seq_along(fgsea_res), function(k) {
      x <- as.data.table(fgsea_res[[k]])
      x[, `:=`(
        DE = DE,
        dataset_id = dat_names[k],
        direction = as.character(pathway)   # "UP" / "DN"
      )]
      x
    }),
    use.names = TRUE, fill = TRUE
  )
  gsea_dt<-rbind(gsea_dt, gsea_tmp)
}

p <- as.numeric(gsea_dt$padj)

# effetto per decidere inversione: in fgsea hai NES
eff <- as.numeric(gsea_dt$NES)

inv_std <- fifelse(
  gsea_dt$direction == "UP", eff < 0,
  fifelse(gsea_dt$direction == "DN", eff > 0, NA)
)

inv_inv <- fifelse(
  gsea_dt$direction == "UP", eff > 0,
  fifelse(gsea_dt$direction == "DN", eff < 0, NA)
)

gsea_dt[, `:=`(p_1t = NA_real_, p_1t_inv = NA_real_,
               inverted_1t = as.logical(inv_std),
               inverted_1t_inv = as.logical(inv_inv))]

ok_std <- !is.na(p) & !is.na(inv_std)
ok_inv <- !is.na(p) & !is.na(inv_inv)

gsea_dt$p_1t[ok_std] <- metap::two2one(
  p[ok_std],
  two    = rep(TRUE, sum(ok_std)),
  invert = inv_std[ok_std]
)

gsea_dt$p_1t_inv[ok_inv] <- metap::two2one(
  p[ok_inv],
  two    = rep(TRUE, sum(ok_inv)),
  invert = inv_inv[ok_inv]
)

write.xlsx(gsea_dt, file=file.path(OUT, "all_GSEA_invivo.xlsx"))



df_rev<-data.frame(pup=-log10(unlist(p_up_rev)), pdn=-log10(unlist(p_dn_rev)), dataset=c( DEGs_seq))

df_diff<-data.frame(pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_seq)))

df<-data.frame(pup=-log10(unlist(p_up)), pdn=-log10(unlist(p_dn)), dataset=c(DEGs_seq))


df_all<-data.frame(pup_rev=df_rev$pup, pdn_rev=df_rev$pdn, pup=df$pup, pdn=df$pdn, pup_diff= -df$pup+df_rev$pup, pdn_diff= -df$pdn+df_rev$pdn, dataset=gsub("DE_", "", c(DEGs_seq)))
df_all$pup_diff= -df$pup+df_rev$pup 
df_all$pdn_diff= -df$pdn+df_rev$pdn
df_all[df_all=="-Inf"]<-(-300)
df_all[df_all< (-300)]<-(-300)
df_all[df_all=="Inf"]<-(300)
df_all$average_APR <- (df_all$pup_diff + df_all$pdn_diff)/2

library(openxlsx)
write.xlsx(df_all, file=file.path(OUT, "invivo_metap_all.xlsx"))


##########
## GO for LSD
#########

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down<-list()
ego_up<-list()
i<-0
for(DE in DEGs_seq){
  i<-i+1
  homologs<-mouse_homologs
  
  
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs[[P_SEQ]]<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs[[P_SEQ]]<0.05)]
  
  universe<-rownames(DEGs)
  
  ego_down[[i]]<- enrichGO(gene = genes_dn,
                    keyType="SYMBOL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    universe=universe, OrgDb="org.Mm.eg.db")

ego_up[[i]]<- enrichGO(gene = genes_up,
                  keyType="SYMBOL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  universe=universe, OrgDb="org.Mm.eg.db")

}

save(ego_up, file=file.path(OUT, "RData/ego_up_valid.RData"))
save(ego_down, file=file.path(OUT, "RData/ego_down_valid.RData"))