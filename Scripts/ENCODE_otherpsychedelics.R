################
#### ENCODE
################
encode <- read.gmt("Data/ENCODE_TF_ChIP-seq_2015.txt")
#encode<-encode[grep("ENCODE", encode[,1]),]


load("Results/RData/fgsea_TF_ENCODE.RData")

allp<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allp)<-unique(encode[,1])
allNES<-matrix(NA, nrow=length(unique(encode[,1])), ncol=length(fgsea_TF))
rownames(allNES)<-unique(encode[,1])
for(i in 1:length(fgsea_TF)){
  allp[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$p.adjust
  allNES[fgsea_TF[[i]]$ID,i]<-fgsea_TF[[i]]$NES
}

mouse_data<-c("DE_MDMA", "DE_GSE161626","DE_GSE161626_48h","DE_GSE161626_7d","DE_GSE81672","DE_GSE26364","DE_GSE209859",
              "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks")
rat_data<-c("DE_GSE14720", "DE_GSE23728","DE_GSE179380","DE_DMT", "DE_pharm", "DE_harm")
##MDMA
DEGs_array<-c("DE_MDMA",
              
              "DE_GSE26364",
              #DOI
              
              "DE_GSE23728",
              #LSD
              "DE_GSE179380")

DEGs_seq<-c("DE_DMT", "DE_GSE161626","DE_GSE161626_48h","DE_GSE161626_7d", "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
            "DE_GSE81672", "DE_harm", "DE_pharm")

TF_down<-list()
TF_up<-list()
i<-0
for(DE in DEGs_array){
  i<-i+1
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$logFC>0 & DEGs$P.Value<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$logFC<0 & DEGs$P.Value<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
universe<-rownames(DEGs)
universe<-unique(homologs[which(homologs[,2] %in% universe),3])

TF_down[[i]]<- enricher(gene = genes_dn,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=universe, 
                   TERM2GENE = encode)

TF_up[[i]]<- enricher(gene = genes_up,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1,
                 universe=universe, TERM2GENE = encode)
}

for(DE in DEGs_seq){
  i<-i+1
  if(DE %in% mouse_data){
    homologs<-mouse_homologs
  } else if(DE %in% rat_data){
    homologs<-rat_homologs
  } else {
    print(DE)
  }
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$pvalue<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$pvalue<0.05)]
  genes_up<-unique(homologs[which(homologs[,2] %in% genes_up),3])
  genes_dn<-unique(homologs[which(homologs[,2] %in% genes_dn),3])
  
  universe<-rownames(DEGs)
  universe<-unique(homologs[which(homologs[,2] %in% universe),3])
  
  TF_down[[i]]<- enricher(gene = genes_dn,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 1,
                          qvalueCutoff  = 1,
                          universe=universe, 
                          TERM2GENE = encode)
  
  TF_up[[i]]<- enricher(gene = genes_up,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        universe=universe, TERM2GENE = encode)
}

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

DEG_down_conv<-unique(rat_homologs[which(rat_homologs[,2] %in% DEG_down),3])
DEG_up_conv<-unique(rat_homologs[which(rat_homologs[,2] %in% DEG_up),3])

universe<-rownames(DE_GSE179379)
universe<-unique(rat_homologs[which(rat_homologs[,2] %in% universe),3])

i<-i+1
TF_down[[i]]<- enricher(gene = DEG_down_conv,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   universe=universe, 
                   TERM2GENE = encode)

TF_up[[i]]<- enricher(gene = DEG_up_conv,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1,
                 universe=universe, TERM2GENE = encode)



names(TF_down)<-c(DEGs_array, DEGs_seq, "LSD")
names(TF_up)<-c(DEGs_array, DEGs_seq, "LSD")

TF<-"TAF1"
p_TAF1<-c()
dat_TAF1<-c()
for(i in 1:length(TF_up)){
  p_TAF1<-c(p_TAF1, TF_up[[i]]$pvalue[grep(TF, TF_up[[i]]$ID)])
  dat_TAF1<-c(dat_TAF1, rep(names(TF_up)[i], length(TF_up[[i]]$pvalue[grep(TF, TF_up[[i]]$ID)])))
}
names(p_TAF1)<-dat_TAF1
which(p_TAF1<0.05)
