##EXERCISE

ctrl<-which(GSE64607_meta$`tissue:ch1`=="lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1`=="saline vehicle" & GSE64607_meta$`timepoint:ch1`=="day 28")
trt<-which(GSE64607_meta$`tissue:ch1`=="lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1`=="voluntary running" & GSE64607_meta$`timepoint:ch1`=="day 28")

DE_GSE64607<-DExpr(GSE64607[,trt], GSE64607[,ctrl])#FC positivi sono più alti nel trattato

GSE64607_dn<-rownames(DE_GSE64607)[DE_GSE64607$logFC<0 & DE_GSE64607$P.Value<0.05]
GSE64607_up<-rownames(DE_GSE64607)[DE_GSE64607$logFC>0 & DE_GSE64607$P.Value<0.05]

####TOGLIERE
ctrl<-which(GSE38465_meta$`strain:ch1`=="SAMP8" & GSE38465_meta$`treatment:ch1`=="sedentary")
trt<-which(GSE38465_meta$`strain:ch1`=="SAMP8" & GSE38465_meta$`treatment:ch1`=="exercised")

DE_GSE38465_SAMP8<-DExpr(GSE38465[,trt], GSE38465[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-which(GSE38465_meta$`strain:ch1`=="SAMR1" & GSE38465_meta$`treatment:ch1`=="sedentary")
trt<-which(GSE38465_meta$`strain:ch1`=="SAMR1" & GSE38465_meta$`treatment:ch1`=="sedentary")

DE_GSE38465_SAMR1<-DExpr(GSE38465[,trt], GSE38465[,ctrl])#FC positivi sono più alti nel trattato


GSE38465_SAMP8_dn<-rownames(DE_GSE38465_SAMP8)[DE_GSE38465_SAMP8$logFC<0 & DE_GSE38465_SAMP8$P.Value<0.05]
GSE38465_SAMP8_up<-rownames(DE_GSE38465_SAMP8)[DE_GSE38465_SAMP8$logFC>0 & DE_GSE38465_SAMP8$P.Value<0.05]

GSE38465_SAMR1_dn<-rownames(DE_GSE38465_SAMR1)[DE_GSE38465_SAMR1$logFC<0 & DE_GSE38465_SAMR1$P.Value<0.05]
GSE38465_SAMR1_up<-rownames(DE_GSE38465_SAMR1)[DE_GSE38465_SAMR1$logFC>0 & DE_GSE38465_SAMR1$P.Value<0.05]
####FIN QUI

ctrl<-grep("Wild-type standard cage", GSE30880_meta$title)
trt<-grep("Wild-type Env Enrichment", GSE30880_meta$title)

DE_GSE30880<-DExpr(GSE30880[,trt], GSE30880[,ctrl])#FC positivi sono più alti nel trattato


ctrl<-intersect(grep("control", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="0 hour post final vapor chamber session"))
trt<-intersect(grep("CIE", GSE72507_meta[,1]),which(GSE72507_meta$`sacrifice time:ch1`=="0 hour post final vapor chamber session"))

DE_GSE72507_0<-DExpr(GSE72507[,trt], GSE72507[,ctrl])#FC positivi sono più alti nel trattato

DE_GSE60676_0<-DExpr(GSE60676[,trt], GSE60676[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-intersect(grep("PFC_8hr_control", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )
trt<-intersect(grep("PFC_8hr_treated", GSE60676_meta[,1]), which(GSE60676_meta$`tissue:ch1`=="prefrontal cortex") )

ctrl<-which(GSE28515_meta$source_name_ch1=="saline_C57")
trt<-which(GSE28515_meta$source_name_ch1=="etoh_C57")

DE_GSE28515_C57<-DExpr(GSE28515[,trt], GSE28515[,ctrl])#FC positivi sono più alti nel trattato

GSE105453<-log2(GSE105453+1)

ctrl<-which(GSE105453_meta$`age:ch1`=="young" & GSE105453_meta$`cognitive performance:ch1`=="nonperforming")
trt<-which(GSE105453_meta$`age:ch1`=="young" & GSE105453_meta$`cognitive performance:ch1`=="performing")

DE_GSE105453_Y<-DExpr(GSE105453[,trt], GSE105453[,ctrl])#FC positivi sono più alti nel trattato

ctrl<-which(GSE105453_meta$`age:ch1`=="aged" & GSE105453_meta$`cognitive performance:ch1`=="nonperforming")
trt<-which(GSE105453_meta$`age:ch1`=="aged" & GSE105453_meta$`cognitive performance:ch1`=="performing")

DE_GSE105453_O<-DExpr(GSE105453[,trt], GSE105453[,ctrl])#FC positivi sono più alti nel trattato

library(DESeq2)
GSE111273_age<-GSE111273_meta$`age:ch1`
GSE111273_trt<-GSE111273_meta$`perturbation:ch1`

dds <- DESeqDataSetFromMatrix(countData = GSE111273[,GSE111273_age %in% c("5 months")],
                              colData = data.frame(treatment=GSE111273_trt[GSE111273_age %in% c("5 months")]),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE111273_Y <- results(dds, contrast = c("treatment","enriched environment", "none"))

dds <- DESeqDataSetFromMatrix(countData = GSE111273[,GSE111273_age %in% c("5 months")],
                              colData = data.frame(treatment=GSE111273_trt[GSE111273_age %in% c("24 months")]),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE111273_O <- results(dds, contrast = c("treatment","enriched environment", "none"))


GSE164798_trt<-GSE164798_meta$`treatment:ch1`
GSE164798_strain<-GSE164798_meta$`strain:ch1`

dds <- DESeqDataSetFromMatrix(countData = GSE164798[,GSE164798_strain %in% c("C57BL/6")],
                              colData = data.frame(treatment=GSE164798_trt[GSE164798_strain %in% c("C57BL/6")]),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE164798_chronic <- results(dds, contrast = c("treatment","30 days voluntary saucer wheel running", "untreated"))
DE_GSE164798_acute <- results(dds, contrast = c("treatment","one day voluntary saucer wheel running", "untreated"))



library(fgsea)
library(ggplot2)
fgsea_res<-list()
i<-0
for(DE in c("DE_GSE64607","DE_GSE38465_SAMP8","DE_GSE38465_SAMR1","DE_GSE72507_0", "DE_GSE60676_0",  "DE_GSE28515_C57")){
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
    fgsea_res[[i]][[n]]<-fgsea(list(UP=genes_up, DN=genes_dn), forgesea, nPermSimple=100000)

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


