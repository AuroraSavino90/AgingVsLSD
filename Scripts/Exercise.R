
library(GEOquery)
library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD")
#preprocessing
#change names in gene symbols

changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]

  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }

  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
      }
    }

    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }

  rownames(data)=annotation_sel
  return(data)
}

library(limma)
DExpr<-function(cond1, cond2){
  cond1<-data.frame(cond1)
  cond2<-data.frame(cond2)
  gs <- factor(c(rep("Cond1",ncol(cond1)),rep("Cond2",ncol(cond2))))
  #}
  design <- model.matrix(~gs + 0, cbind.data.frame(cond1, cond2))
  #colnames(design) <- levels(gs)

  fit <- lmFit(cbind.data.frame(cond1, cond2), design)  # fit linear model

  # set up contrasts of interest and recalculate model coefficients
  #cts <- paste(levels(gs)[1], levels(gs)[2], sep="-")
  cts<-"gsCond1-gsCond2"
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)

  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  return(tT)
}


################################
##########GSE132930 #NOT CLEAR WHAT THEY DID
###############################

gds<-getGEO("GSE132930",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE132930<-exprs(gds$GSE132930_series_matrix.txt.gz)
GSE132930_meta<-pData(gds$GSE132930_series_matrix.txt.gz)
GSE132930_anno<-fData(gds$GSE132930_series_matrix.txt.gz)

GSE132930_Y<-read.csv("Data/Exercise/GSE132930_GeneCount_C_young.txt.gz", sep="\t", row.names=1)
colnames(GSE132930_Y)<-unlist(strsplit(colnames(GSE132930_Y), "_"))[seq(2,ncol(GSE132930_Y)*7,7)]
GSE132930_O<-read.csv("Data/Exercise/GSE132930_GeneCount_C_midlife.txt.gz", sep="\t", row.names=1)
colnames(GSE132930_O)<-unlist(strsplit(colnames(GSE132930_O), "_"))[seq(2,ncol(GSE132930_O)*7,7)]

match(colnames(GSE132930_O), GSE132930_meta$title)


DEGs<-read.csv("Data/Exercise/GSE132930_Early_Diet_cortex_8wkrun_CD.8weeks.Run.C_vs_CD.8weeks.Sed.C.txt", sep="\t", row.names=1)

################################
##########GSE64607
###############################

gds<-getGEO("GSE64607",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE64607<-exprs(gds$GSE64607_series_matrix.txt.gz)
GSE64607_meta<-pData(gds$GSE64607_series_matrix.txt.gz)
GSE64607_anno<-fData(gds$GSE64607_series_matrix.txt.gz)

GSE64607<-changenames(data=GSE64607, anno=cbind(GSE64607_anno$ID, GSE64607_anno$`Gene symbol`))

ctrl<-which(GSE64607_meta$`tissue:ch1`=="lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1`=="saline vehicle" & GSE64607_meta$`timepoint:ch1`=="day 28")
trt<-which(GSE64607_meta$`tissue:ch1`=="lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1`=="voluntary running" & GSE64607_meta$`timepoint:ch1`=="day 28")

DE_GSE64607<-DExpr(GSE64607[,trt], GSE64607[,ctrl])#FC positivi sono più alti nel trattato

exercise_dn<-rownames(DE_GSE64607)[DE_GSE64607$logFC<0 & DE_GSE64607$P.Value<0.05]
exercise_up<-rownames(DE_GSE64607)[DE_GSE64607$logFC>0 & DE_GSE64607$P.Value<0.05]

inboth<-intersect(rownames(DE_GSE64607), rownames(DE_GSE179379))
plot(DE_GSE179379[inboth,"log2FoldChange"],DE_GSE64607[inboth,"logFC"])

################################
##########GSE38465
###############################

gds<-getGEO("GSE38465",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE38465<-exprs(gds$GSE38465_series_matrix.txt.gz)
GSE38465_meta<-pData(gds$GSE38465_series_matrix.txt.gz)
GSE38465_anno<-fData(gds$GSE38465_series_matrix.txt.gz)

GSE38465<-GSE38465[rowSums(is.na(GSE38465))==0,]
GSE38465<-changenames(data=GSE38465, anno=cbind(GSE38465_anno$ID, GSE38465_anno$`Gene symbol`))

ctrl<-which(GSE38465_meta$`strain:ch1`=="SAMP8" & GSE38465_meta$`treatment:ch1`=="sedentary")
trt<-which(GSE38465_meta$`strain:ch1`=="SAMP8" & GSE38465_meta$`treatment:ch1`=="exercised")

ctrl<-which(GSE38465_meta$`strain:ch1`=="SAMR1" & GSE38465_meta$`treatment:ch1`=="sedentary")
trt<-which(GSE38465_meta$`strain:ch1`=="SAMR1" & GSE38465_meta$`treatment:ch1`=="exercised")

DE_GSE38465<-DExpr(GSE38465[,trt], GSE38465[,ctrl])#FC positivi sono più alti nel trattato

exercise_dn<-rownames(DE_GSE38465)[DE_GSE38465$logFC<0 & DE_GSE38465$P.Value<0.05]
exercise_up<-rownames(DE_GSE38465)[DE_GSE38465$logFC>0 & DE_GSE38465$P.Value<0.05]

################################
##########GSE205907
###############################

gds<-getGEO("GSE205907",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE205907<-read.xlsx("Data/Exercise/GSE205907_fpkm.xlsx")
GSE205907_meta<-pData(gds$GSE205907_series_matrix.txt.gz)
GSE205907_anno<-fData(gds$GSE205907_series_matrix.txt.gz)

GSE205907<-changenames(data=GSE205907[,-1], anno=cbind(rownames(GSE205907), GSE205907[,1]))

ctrl<-which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP")
trt<-which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP+exercise")

DE_GSE205907<-DExpr(GSE205907[,trt], GSE205907[,ctrl])#FC positivi sono più alti nel trattato

exercise_dn<-rownames(DE_GSE205907)[DE_GSE205907$logFC<0 & DE_GSE205907$P.Value<0.05]
exercise_up<-rownames(DE_GSE205907)[DE_GSE205907$logFC>0 & DE_GSE205907$P.Value<0.05]

################################
##########GSE124954 VISUAL CORTEX
###############################

gds<-getGEO("GSE124954",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE124954<-exprs(gds$GSE124954_series_matrix.txt.gz)
GSE124954_meta<-pData(gds$GSE124954_series_matrix.txt.gz)
GSE124954_anno<-fData(gds$GSE124954_series_matrix.txt.gz)

###########HIPPOCAMPUS

gds<-getGEO("GSE30880",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE30880<-exprs(gds$GSE30880_series_matrix.txt.gz)
GSE30880_meta<-pData(gds$GSE30880_series_matrix.txt.gz)
GSE30880_anno<-fData(gds$GSE30880_series_matrix.txt.gz)

GSE30880<-changenames(data=GSE30880, anno=cbind(GSE30880_anno$ID, GSE30880_anno$`Gene symbol`))

####
gds<-getGEO("GSE105453",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)
GSE105453_meta<-pData(gds$GSE105453_series_matrix.txt.gz)

GSE105453<-read.csv("Data/Exercise/GSE105453_NormalizedCounts.txt", sep="\t")
rownames(GSE105453)<-GSE105453[,1]
GSE105453<-changenames(data=GSE105453[,-c(1:3)], anno=cbind(rownames(GSE105453), GSE105453[,3]))


gds<-getGEO("GSE111273",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)
GSE111273_meta<-pData(gds$GSE111273_series_matrix.txt.gz)

GSE111273<-read.csv("Data/Exercise/GSE111273_sample_counts.txt", sep="\t")
rownames(GSE111273)<-GSE111273[,1]

library(biomaRt)
ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(GSE111273),
           mart = ensembl.mouse)

GSE111273<-changenames(data=GSE111273[,-c(1)], anno=ids)

gds<-getGEO("GSE164798",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)
GSE164798_meta<-pData(gds$GSE164798_series_matrix.txt.gz)

GSE164798<-read.csv("Data/Exercise/GSE164798_Raw_gene_counts_matrix.txt", sep="\t")
rownames(GSE164798)<-GSE164798[,1]
GSE164798<-GSE164798[,-c(1:6)]
################################
##########GSE74010 #modafinil
###############################

gds<-getGEO("GSE74010",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE74010<-exprs(gds$GSE74010_series_matrix.txt.gz)
GSE74010_meta<-pData(gds$GSE74010_series_matrix.txt.gz)
GSE74010_anno<-fData(gds$GSE74010_series_matrix.txt.gz)

GSE74010<-changenames(data=GSE74010, anno=cbind(GSE74010_anno$ID, GSE74010_anno$`Gene symbol`))

ctrl<-which(GSE74010_meta$`tissue subtype:ch1`=="prefrontal cortex (PFC)" & GSE74010_meta$`injected with:ch1`=="saline during conditioning")
trt<-which(GSE74010_meta$`tissue subtype:ch1`=="prefrontal cortex (PFC)" & GSE74010_meta$`injected with:ch1`=="modafinil 65 mg/kg during conditioning")

DE_GSE74010<-DExpr(GSE74010[,trt], GSE74010[,ctrl])#FC positivi sono più alti nel trattato

exercise_dn<-rownames(DE_GSE74010)[DE_GSE74010$logFC<0 & DE_GSE74010$P.Value<0.05]
exercise_up<-rownames(DE_GSE74010)[DE_GSE74010$logFC>0 & DE_GSE74010$P.Value<0.05]


load("RData/DE_GSE179379.RData")
homologs<-read.csv("Data/Human rat homologs.txt")

genes_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
genes_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]

inboth<-intersect(DEGs$Associated.Gene.Name, rownames(DE_GSE179379))

plot(DE_GSE179379[inboth,"log2FoldChange"],DEGs$logFC_CD.8weeks.Run.C_vs_CD.8weeks.Sed.C[match(inboth, DEGs$Associated.Gene.Name)])
exercise_dn<-DEGs$Associated.Gene.Name[DEGs$logFC_CD.8weeks.Run.C_vs_CD.8weeks.Sed.C<0 & DEGs$FDR<0.05]
exercise_up<-DEGs$Associated.Gene.Name[DEGs$logFC_CD.8weeks.Run.C_vs_CD.8weeks.Sed.C>0 & DEGs$FDR<0.05]

LSD_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]
LSD_dn<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]


a<-length(intersect(exercise_dn, LSD_dn))
b<-length(setdiff(exercise_dn, LSD_dn))
c<-length(setdiff(LSD_dn, exercise_dn))
d<-length(setdiff(intersect(rownames(DE_GSE179379), DEGs$Associated.Gene.Name),
                  union(exercise_dn, LSD_dn)))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

a<-length(intersect(exercise_up, LSD_up))
b<-length(setdiff(exercise_up, LSD_up))
c<-length(setdiff(LSD_up, exercise_up))
d<-length(setdiff(intersect(rownames(DE_GSE179379), DEGs$Associated.Gene.Name),
                  union(exercise_up, LSD_up)))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

a<-length(intersect(exercise_dn, LSD_up))
b<-length(setdiff(exercise_dn, LSD_up))
c<-length(setdiff(LSD_up, exercise_dn))
d<-length(setdiff(intersect(rownames(DE_GSE179379), DEGs$Associated.Gene.Name),
                  union(exercise_dn, LSD_up)))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")


intersect(exercise_dn, LSD_up)
intersect(exercise_up, LSD_dn)

########################################
load("RData/Alldata_7Jul.RData")
homologs<-read.csv("Data/Human mouse homologs.txt", sep="\t")
genes_up<-exercise_up
genes_dn<-exercise_dn
genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

library(fgsea)
library(ggplot2)
region<-c("DLPFC")
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance
dat<-setdiff(dat, c("GSE5388"))#outlier

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
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]

    genes_cor[[n]]<-cor(t(data), as.numeric(age))
    dat_names<-c(dat_names, dd)
  }
}

region<-"PFC"
diagnosis<-"Healthy"
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region]))
dat<-setdiff(dat, c("GSE102741"))#batch PC1>40% of variance

for(dd in dat){
  n<-n+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  data<-data[, sample]
  age<-metadata$Age[which(metadata$Dataset==dd &metadata$Organism=="Homo sapiens" & metadata$Diagnosis==diagnosis & metadata$Region_simpl %in% region)]
  age<-as.numeric(age)
  if(length(which(age>6570))>=10){
    sample<-sample[which(age>6570)]
    data<-data[, sample]
    age<-age[which(age>6570)]

    genes_cor[[n]]<-cor(t(data), as.numeric(age))
    dat_names<-c(dat_names, dd)
  }
}
#LSD
rat_homologs<-unique(homologs[homologs[,2] %in% DEGs$Associated.Gene.Name,3])

genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
fgsea_res<-list()
p<-list()
df<-list()
for(n in 1:length(genes_cor)){

  forgesea<-unlist(genes_cor[[n]])
  names(forgesea)<-rownames(genes_cor[[n]])
  forgesea<-forgesea[!is.na(forgesea)]
  forgesea<-forgesea[names(forgesea) %in% rat_homologs]
  fgsea_res[[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)

  p1<-plotEnrichment(genes_up, forgesea)
  p2<-plotEnrichment(genes_dn, forgesea)

  df[[n]]<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dat_names[n],nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dat_names[n],nrow(p2$data))))

  p[[n]]<-ggplot(df[[n]], aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))

}


df_tot<-df[[1]]
for(n in 2:length(genes_cor)){
  df_tot<-rbind.data.frame(df_tot, df[[n]])
}


fgsea_rank<-lapply(fgsea_res, function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank), decreasing=T)])

ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x")

genes_cor_up<-matrix(nrow=length(genes_up), ncol=length(genes_cor))
rownames(genes_cor_up)<-genes_up
for(n in 1:length(genes_cor)){
  genes_cor_up[rownames(genes_cor[[n]])[which(rownames(genes_cor[[n]])%in% genes_up)],n]<-genes_cor[[n]][which(rownames(genes_cor[[n]])%in% genes_up),]
}

genes_cor_up<-genes_cor_up[rowSums(is.na(genes_cor_up))<5,]

pheatmap(genes_cor_up)
sort(rowMeans(genes_cor_up, na.rm=T))

rownames(DE_GSE64607)<-toupper(rownames(DE_GSE64607))
plot(DE_GSE64607[rownames(genes_cor_up),"logFC"],rowMeans(genes_cor_up, na.rm=T))

##################
meta_dn<-function(x){
  istwo <- rep(T, 14)
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
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

pup<-meta_up(fgsea_res)
pdn<-meta_dn(fgsea_res)
