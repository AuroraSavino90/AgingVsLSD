################################
##########GSE72507 #alcohol
###############################

gds<-getGEO("GSE72507",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE72507<-exprs(gds$GSE72507_series_matrix.txt.gz)
GSE72507_meta<-pData(gds$GSE72507_series_matrix.txt.gz)
GSE72507_anno<-fData(gds$GSE72507_series_matrix.txt.gz)

GSE72507<-changenames(data=GSE72507, anno=cbind(GSE72507_anno$ID, GSE72507_anno$`Gene symbol`))

ctrl<-intersect(grep("control", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="0 hour post final vapor chamber session"))
trt<-intersect(grep("CIE", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="0 hour post final vapor chamber session"))

DE_GSE72507_0<-DExpr(GSE72507[,trt], GSE72507[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-intersect(grep("control", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="8 hours post final vapor chamber session"))
trt<-intersect(grep("CIE", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="8 hours post final vapor chamber session"))

DE_GSE72507_8<-DExpr(GSE72507[,trt], GSE72507[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-intersect(grep("control", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="7 days post final vapor chamber session"))
trt<-intersect(grep("CIE", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="7 days post final vapor chamber session"))

DE_GSE72507_7<-DExpr(GSE72507[,trt], GSE72507[,ctrl])#FC positivi sono più alti nel trattato

################################
##########GSE143419 #alcohol
###############################

gds<-getGEO("GSE143419",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE143419<-exprs(gds$GSE143419_series_matrix.txt.gz)
GSE143419_meta<-pData(gds$GSE143419_series_matrix.txt.gz)
GSE143419_anno<-fData(gds$GSE143419_series_matrix.txt.gz)

GSE143419<-changenames(data=GSE143419, anno=cbind(GSE143419_anno$ID, GSE143419_anno$`Gene symbol`))

ctrl<-which(GSE143419_meta$`brain region:ch1`=="prefrontal cortex" & GSE143419_meta$`vapor treatment:ch1`=="Control"&GSE143419_meta$`drinking treatment:ch1`=="Nondrinker")
trt<-which(GSE143419_meta$`brain region:ch1`=="prefrontal cortex" & GSE143419_meta$`vapor treatment:ch1`=="CIE"&GSE143419_meta$`drinking treatment:ch1`=="Drinker")

DE_GSE143419<-DExpr(GSE143419[,trt], GSE143419[,ctrl])#FC positivi sono più alti nel trattato

################################
##########GSE60676 #alcohol
###############################

gds<-getGEO("GSE60676",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE60676<-exprs(gds$GSE60676_series_matrix.txt.gz)
GSE60676_meta<-pData(gds$GSE60676_series_matrix.txt.gz)
GSE60676_anno<-fData(gds$GSE60676_series_matrix.txt.gz)

GSE60676<-changenames(data=GSE60676, anno=cbind(GSE60676_anno$ID, GSE60676_anno$`Gene symbol`))

ctrl<-intersect(grep("PFC_0hr_control", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )
trt<-intersect(grep("PFC_0hr_treated", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )

DE_GSE60676_0<-DExpr(GSE60676[,trt], GSE60676[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-intersect(grep("PFC_8hr_control", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )
trt<-intersect(grep("PFC_8hr_treated", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )

DE_GSE60676_8<-DExpr(GSE60676[,trt], GSE60676[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-intersect(grep("PFC_120hr_control", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )
trt<-intersect(grep("PFC_120hr_treated", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )

DE_GSE60676_120<-DExpr(GSE60676[,trt], GSE60676[,ctrl])#FC positivi sono più alti nel trattato


################################
##########GSE176122 #TODO
###############################

gds<-getGEO("GSE176122",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE176122<-exprs(gds$GSE176122_series_matrix.txt.gz)
GSE176122_meta<-pData(gds$GSE176122_series_matrix.txt.gz)
GSE176122_anno<-fData(gds$GSE176122_series_matrix.txt.gz)

GSE176122<-changenames(data=GSE176122, anno=cbind(GSE176122_anno$ID, GSE176122_anno$`Gene symbol`))

ctrl<-which(GSE93311_meta$`brain region:ch1`=="prefrontal cortex" & GSE93311_meta$`vapor treatment:ch1`=="Control"&GSE93311_meta$`drinking treatment:ch1`=="Nondrinker")
trt<-which(GSE93311_meta$`brain region:ch1`=="prefrontal cortex" & GSE93311_meta$`vapor treatment:ch1`=="CIE"&GSE93311_meta$`drinking treatment:ch1`=="Drinker")

DE_GSE93311<-DExpr(GSE93311[,trt], GSE93311[,ctrl])#FC positivi sono più alti nel trattato

################################
##########GSE28515 #alcohol
###############################

gds<-getGEO("GSE28515",destdir="/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVsLSD", AnnotGPL = TRUE)

GSE28515<-exprs(gds$GSE28515_series_matrix.txt.gz)
GSE28515_meta<-pData(gds$GSE28515_series_matrix.txt.gz)
GSE28515_anno<-fData(gds$GSE28515_series_matrix.txt.gz)

GSE28515<-changenames(data=GSE28515, anno=cbind(GSE28515_anno$ID, GSE28515_anno$`Gene symbol`))

ctrl<-which(GSE28515_meta$source_name_ch1=="saline_DBA")
trt<-which(GSE28515_meta$source_name_ch1=="etoh_DBA")

DE_GSE28515_DBA<-DExpr(GSE28515[,trt], GSE28515[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-which(GSE28515_meta$source_name_ch1=="saline_BXD")
trt<-which(GSE28515_meta$source_name_ch1=="etoh_BXD")

DE_GSE28515_BXD<-DExpr(GSE28515[,trt], GSE28515[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-which(GSE28515_meta$source_name_ch1=="saline_C57")
trt<-which(GSE28515_meta$source_name_ch1=="etoh_C57")

DE_GSE28515_C57<-DExpr(GSE28515[,trt], GSE28515[,ctrl])#FC positivi sono più alti nel trattato

######################



GSE56249
GSE123114
#########################




GSE72507_dn<-rownames(DE_GSE72507)[DE_GSE72507$logFC<0 & DE_GSE72507$P.Value<0.05]
GSE72507_up<-rownames(DE_GSE72507)[DE_GSE72507$logFC>0 & DE_GSE72507$P.Value<0.05]

GSE143419_dn<-rownames(DE_GSE143419)[DE_GSE143419$logFC<0 & DE_GSE143419$P.Value<0.05]
GSE143419_up<-rownames(DE_GSE143419)[DE_GSE143419$logFC>0 & DE_GSE143419$P.Value<0.05]

GSE60676_dn<-rownames(DE_GSE60676)[DE_GSE60676$logFC<0 & DE_GSE60676$P.Value<0.05]
GSE60676_up<-rownames(DE_GSE60676)[DE_GSE60676$logFC>0 & DE_GSE60676$P.Value<0.05]

GSE28515_dn<-rownames(DE_GSE28515)[DE_GSE28515$logFC<0 & DE_GSE28515$P.Value<0.05]
GSE28515_up<-rownames(DE_GSE28515)[DE_GSE28515$logFC>0 & DE_GSE28515$P.Value<0.05]

Reduce(intersect, list(GSE72507_dn, GSE143419_dn, GSE60676_dn))

Reduce(intersect, list(GSE72507_up, GSE143419_up, GSE60676_up))


inall<-c(rownames(DE_GSE72507_0), rownames(DE_GSE72507_8), rownames(DE_GSE72507_7),
         rownames(DE_GSE143419), rownames(DE_GSE60676_0), rownames(DE_GSE60676_8),
         rownames(DE_GSE60676_120),
         rownames(DE_GSE28515_DBA), rownames(DE_GSE28515_BXD), rownames(DE_GSE28515_C57))

allDE<-cbind(DE_GSE72507_0[inall,"logFC"], DE_GSE72507_8[inall,"logFC"],
             DE_GSE72507_7[inall,"logFC"], DE_GSE143419[inall,"logFC"],
             DE_GSE60676_0[inall,"logFC"],DE_GSE60676_8[inall,"logFC"],
             DE_GSE60676_120[inall,"logFC"],DE_GSE28515_DBA[inall,"logFC"],
             DE_GSE28515_BXD[inall,"logFC"],DE_GSE28515_C57[inall,"logFC"])

corrs<-cor(allDE, use="pairwise.complete.obs")
pheatmap(corrs)

rownames(corrs)<-c("DE_GSE72507_0", "DE_GSE72507_8", "DE_GSE72507_7",
                   "DE_GSE143419", "DE_GSE60676_0", "DE_GSE60676_8",
                   "DE_GSE60676_120", "DE_GSE28515_DBA",
                   "DE_GSE28515_BXD", "DE_GSE28515_C57")
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(corrs), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(corrs), na.rm=T)/paletteLength, max(unlist(corrs), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1
pheatmap(corrs,   cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, keep.dendro=T)

##################
#################
#non tornano DBA,GSE60676_8,GSE143419
library(fgsea)
library(ggplot2)
fgsea_res<-list()
i<-0
for(DE in c("DE_GSE72507_0", "DE_GSE72507_8", "DE_GSE72507_7",
            "DE_GSE143419", "DE_GSE60676_0", "DE_GSE60676_8",
            "DE_GSE60676_120", "DE_GSE28515_DBA",
            "DE_GSE28515_BXD", "DE_GSE28515_C57")){
  i<-i+1
  DEGs<-get(DE)
  exercise_dn<-rownames(DEGs)[DEGs$logFC<0 & DEGs$P.Value<0.05]
  exercise_up<-rownames(DEGs)[DEGs$logFC>0 & DEGs$P.Value<0.05]

  genes_up<-exercise_up
  genes_dn<-exercise_dn
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])

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
  rat_homologs<-unique(homologs[homologs[,2] %in% rownames(DEGs),3])

  genes_cor<-genes_cor[!unlist(lapply(genes_cor, is.null))]
  fgsea_res[[i]]<-list()
  p<-list()
  df<-list()
  for(n in 1:length(genes_cor)){

    forgesea<-unlist(genes_cor[[n]])
    names(forgesea)<-rownames(genes_cor[[n]])
    forgesea<-forgesea[!is.na(forgesea)]
    forgesea<-forgesea[names(forgesea) %in% rat_homologs]
    fgsea_res[[i]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea)

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


  fgsea_rank<-lapply(fgsea_res[[i]], function(x){-log10(x[1,3])*sign(x[1,6])+log10(x[2,3])*sign(x[2,6])})

  df_tot$dataset<-factor(df_tot$dataset, levels=dat_names[order(unlist(fgsea_rank[[n]]), decreasing=T)])

  print(ggplot(df_tot, aes(x, y, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))+facet_wrap(~dataset, scale="free_x"))


}


meta_dn_rev<-function(x){
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))

  }

  return(p[[3]])
}

meta_up_rev<-function(x){
  istwo <- rep(T, 14)
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[2,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[2,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))
  }
  return(p[[3]])
}

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

pup_rev<-list()
pdn_rev<-list()
for(n in 1:length(fgsea_res)){

  pup_rev[[n]]<-meta_up_rev(fgsea_res[[n]])
  pdn_rev[[n]]<-meta_dn_rev(fgsea_res[[n]])

}

pup<-list()
pdn<-list()
for(n in 1:length(fgsea_res)){

  pup[[n]]<-meta_up(fgsea_res[[n]])
  pdn[[n]]<-meta_dn(fgsea_res[[n]])

}

plot(-log10(unlist(pup))+log10(unlist(pup_rev)),
-log10(unlist(pdn))+log10(unlist(pdn_rev)))

#GSE143419 troppe variabili nel trattamento. Considerati solo al tempo 0 e geneticamente "normali", non DBA o incroci con DBA.

(-log10(unlist(pup))+log10(unlist(pup_rev)))[c(1,5,10)]
(-log10(unlist(pdn))+log10(unlist(pdn_rev)))[c(1,5,10)]


c("DE_GSE72507_0", "DE_GSE72507_8", "DE_GSE72507_7",
  "DE_GSE143419", "DE_GSE60676_0", "DE_GSE60676_8",
  "DE_GSE60676_120", "DE_GSE28515_DBA",
  "DE_GSE28515_BXD", "DE_GSE28515_C57")
#

