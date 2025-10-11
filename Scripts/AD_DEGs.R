library(metap)
library(limma)
library(clusterProfiler)
library(fgsea)
Abeta_genes_up<-read.csv(file="results/Abeta_genes_up.csv")[,2]
Abeta_genes_dn<-read.csv(file="results/Abeta_genes_dn.csv")[,2]

mouse_homologs<-read.csv("Data/Human mouse homologs.txt", sep="\t")
Abeta_genes_up<-unique(mouse_homologs[which(mouse_homologs[,2] %in% Abeta_genes_up),3])
Abeta_genes_dn<-unique(mouse_homologs[which(mouse_homologs[,2] %in% Abeta_genes_dn),3])



load("Data/Dementia_alldata_meta.RData")

region<-"DLPFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))

top_genes<-list()
i<-0
sets<-c()
p<-list()
fgsea_res<-list()
for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
  sets<-c(sets, dd)
  
  design <- model.matrix(~ 0+disease)
  colnames(design) <- c("Healthy", "Alzheimer")
  cont_matrix <- makeContrasts(ADvsControl = Alzheimer-Healthy, levels=design)
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(data, design)
  # Compute contrast
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  # Bayes statistics of differential expression
  # *There are several options to tweak!*
  fit_contrast <- eBayes(fit_contrast)
  # Generate a vocalno plot to visualize differential expression
  volcanoplot(fit_contrast)
  # Generate a list of top 100 differentially expressed genes
  top_genes[[i]] <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
  logFC<-topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")$logFC
names(logFC)<-rownames(topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH"))
  

fgsea_res[[i]]<-fgsea(list(UP=Abeta_genes_up, DN=Abeta_genes_dn), logFC)

p1<-plotEnrichment(Abeta_genes_up, logFC)
p2<-plotEnrichment(Abeta_genes_dn, logFC)

df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dd,nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dd,nrow(p2$data))))

p[[i]]<-ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
  scale_colour_manual(values=c("dn"="blue", "up"="red"))



}

region<-"PFC"
diagnosis<-c("Healthy", "Alzheimer")
dat<-na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis & metadata$Region %in% region]))

for(dd in dat){
  i<-i+1
  data<-get(dd)
  sample<-metadata$Sample[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  data<-data[, sample]
  disease<-metadata$Diagnosis[which(metadata$Dataset==dd  & metadata$Diagnosis %in% diagnosis & metadata$Region %in% region)]
  disease<-factor(disease, levels=c("Healthy", "Alzheimer"))
  sets<-c(sets, dd)
  
  design <- model.matrix(~ 0+disease)
  colnames(design) <- c("Healthy", "Alzheimer")
  cont_matrix <- makeContrasts(ADvsControl = Alzheimer-Healthy, levels=design)
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(data, design)
  # Compute contrast
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  # Bayes statistics of differential expression
  fit_contrast <- eBayes(fit_contrast)
  # Generate a vocalno plot to visualize differential expression
  volcanoplot(fit_contrast)
  top_genes[[i]] <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
  logFC<-topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")$logFC
  names(logFC)<-rownames(topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH"))
  
  
  fgsea_res[[i]]<-fgsea(list(UP=Abeta_genes_up, DN=Abeta_genes_dn), logFC)
  
  p1<-plotEnrichment(Abeta_genes_up, logFC)
  p2<-plotEnrichment(Abeta_genes_dn, logFC)
  
  df<-rbind.data.frame(data.frame(p1$data, dir=rep("up", nrow(p1$data)), dataset=rep(dd,nrow(p1$data))), data.frame(p2$data, dir=rep("dn", nrow(p2$data)), dataset=rep(dd,nrow(p2$data))))
  
  p[[i]]<-ggplot(df, aes(x=rank, y=ES, colour=dir))+geom_line()+ geom_hline(yintercept = 0, size=0.5)+ theme(panel.background = element_blank())+
    scale_colour_manual(values=c("dn"="blue", "up"="red"))
}


for(i in 1:length(p)){
  pdf(paste("results/Figures/AD",i,"_Abeta_GSEA.pdf", sep=""), 7,7)
  print(p[[i]])
  dev.off()
}