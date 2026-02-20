#####################################################################
##function to change gene names (e.g. from ENSEMBL to gene symbol)
###################################################################

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



###########################
### quality checks: PCA
##########################
library(openxlsx)
library(FactoMineR)
library(ggplot2)
library(cowplot)

metadata<-read.xlsx("Data/In vitro RNA-seq_all/Metadata.xlsx", rowNames = T)

invitro_counts=read.table("Data/In vitro RNA-seq_all/salmon.merged.gene_counts.tsv",header=TRUE)
rownames(invitro_counts)=invitro_counts$gene_id
anno<-invitro_counts[,c(1:2)]
invitro_counts<-changenames(invitro_counts[,-c(1:2)], anno = anno)
invitro_counts=round(invitro_counts)
save(invitro_counts, file="results/invitro_counts.RData")

rownames(metadata)<-colnames(invitro_counts)

invitro_RPM<-log2(t(t(invitro_counts)/colSums(invitro_counts))*1000000 + 1)

pca<-PCA(t(invitro_RPM))

metadata$Treatment<-as.factor(metadata$Treatment)
metadata$Replicate<-as.factor(metadata$Replicate)
metadata$Trt<-rep(c("Ct", "LSD", "Abeta", "Abeta", "Abeta",
                    "Abeta + LSD", "Abeta + LSD", "Abeta + LSD"),3)

df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)

df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)

pdf("results/Figures/invitro_PCA_rep.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC2, colour=Replicate))+geom_point(size=2)+theme_classic()
dev.off()

pdf("results/Figures/invitro_PCA_trt.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC2, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

#######################################################
########### DEGs
#########################################################

#from https://www.biostars.org/p/9570365/,
#from https://www.biostars.org/p/9581212/,
# tximport is scaling the counts from transcripts of Salmon of isoforms of different lengths
# and put together into gene-level. Since isoforms are different in length, tximport corrects for the
# biases in di different counts because of different lengths.
# nfcore/rnaseq does Salmon and then tximport. 
# So: use salmon.merged.gene_counts_length_scaled.tsv file (raw counts estimated after bias length correction)
# and feed this directly into DESeqDataSetFromMatrix() function and round() the numbers to get integers.
#read and prepare raw count table



library(sva)
adjusted <- ComBat_seq(as.matrix(invitro_counts), batch=metadata$Replicate, group=NULL)

invitro_RPM<-log2(t(t(adjusted)/colSums(adjusted))*1000000 + 1)

pca<-PCA(t(invitro_RPM))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)
df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)

pdf("results/Figures/invitro_PCA_rep_postcorr.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC2, colour=Replicate))+geom_point(size=2)+theme_classic()
dev.off()

pdf("results/Figures/invitro_PCA_trt_postcorr12.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC2, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

pdf("results/Figures/invitro_PCA_trt_postcorr34.pdf", 6, 6)
ggplot(df, aes(x=PC3, y=PC4, colour=Trt))+geom_point(size=2)+theme_classic()
dev.off()

invitro_counts<-adjusted[rowSums(adjusted>=10)>=3,]


############## DEGs
library(DESeq2)

invitro_counts_abeta<-invitro_counts[,metadata$Treatment %in% c("Abeta 1uM ","Abeta 500nM ","Abeta 100nM", "Ct"  )]

dds <- DESeqDataSetFromMatrix(countData = invitro_counts_abeta,
                                 colData = metadata[metadata$Treatment %in% c("Abeta 1uM ","Abeta 500nM ","Abeta 100nM", "Ct"  ),],
                                 design= ~ Treatment)
dds <- DESeq(dds)
ddsAbeta1000 <- results(dds, contrast = c("Treatment","Abeta 1uM ", "Ct"))
ddsAbeta500 <- results(dds, contrast = c("Treatment","Abeta 500nM ", "Ct"))
ddsAbeta100 <- results(dds, contrast = c("Treatment","Abeta 100nM", "Ct"))
write.xlsx(data.frame(ddsAbeta1000), file="results/Abeta1000DEGs.xlsx", rowNames=T)
write.xlsx(data.frame(ddsAbeta500), file="results/Abeta500DEGs.xlsx", rowNames=T)
write.xlsx(data.frame(ddsAbeta100), file="results/Abeta100DEGs.xlsx", rowNames=T)

save(dds, file="Results/dds_Abeta.RData")

DEGs_seq<-c("ddsAbeta100", "ddsAbeta500", "ddsAbeta1000")

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down_Abeta_sep<-list()
ego_up_Abeta_sep<-list()
i<-0
for(DE in DEGs_seq){
  i<-i+1
  DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$padj<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$padj<0.05)]
  
  universe<-rownames(DEGs)
  
  ego_down_Abeta_sep[[i]]<- enrichGO(gene = genes_dn,
                           keyType="SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           universe=universe, OrgDb="org.Mm.eg.db")
  
  ego_up_Abeta_sep[[i]]<- enrichGO(gene = genes_up,
                         keyType="SYMBOL",
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         universe=universe, OrgDb="org.Mm.eg.db")
  
}

save(ego_up_Abeta_sep, file="Results/ego_up_Abeta_sep.RData")
save(ego_down_Abeta_sep, file="Results/ego_down_Abeta_sep.RData")


thr<-0.05
Abeta_genes_dn<-unique(c(rownames(ddsAbeta100)[which(ddsAbeta100$padj<thr & ddsAbeta100$log2FoldChange<0)],
                  rownames(ddsAbeta500)[which(ddsAbeta500$padj<thr & ddsAbeta500$log2FoldChange<0)],
                  rownames(ddsAbeta1000)[which(ddsAbeta1000$padj<thr & ddsAbeta1000$log2FoldChange<0)]))
Abeta_genes_up<-unique(c(rownames(ddsAbeta100)[which(ddsAbeta100$padj<thr & ddsAbeta100$log2FoldChange>0)],
               rownames(ddsAbeta500)[which(ddsAbeta500$padj<thr & ddsAbeta500$log2FoldChange>0)],
               rownames(ddsAbeta1000)[which(ddsAbeta1000$padj<thr & ddsAbeta1000$log2FoldChange>0)]))

write.csv(Abeta_genes_up, file="results/Abeta_genes_up_005.csv")
write.csv(Abeta_genes_dn, file="results/Abeta_genes_dn_005.csv")

library(GSVA)
ssgsea_par <- ssgseaParam(as.matrix(invitro_RPM), list(up=Abeta_genes_up,dn=Abeta_genes_dn))  # all other values are default values
ssgsea <- gsva(ssgsea_par)

library(ggside)
df<-data.frame(up=ssgsea[1,], dn=ssgsea[2,], metadata)
df$Trt<-factor(df$Trt, levels=c("Ct", "LSD", "Abeta + LSD", "Abeta"))

pdf("results/Figures/invitro_Abeta_revert_sep_rep_005.pdf", 15, 7)
ggplot(df, aes(x=up, y=dn, colour=as.factor(Trt)))+geom_point(size=3)+facet_grid(~Replicate)+theme_bw()
dev.off()


pdf("results/Figures/invitro_Abeta_revert_sep_005.pdf", 9, 8)
ggplot(df, aes(x=up, y=dn, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()

invitro_RPM_sel<-invitro_RPM[c(Abeta_genes_up, Abeta_genes_dn),]
pca<-PCA(t(invitro_RPM_sel))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)
df$Treatment<-as.factor(df$Treatment)
df$Replicate<-as.factor(df$Replicate)
df$Trt<-as.factor(df$Trt)


pc1_var <- pca$eig[1, 2]
pc2_var <- pca$eig[2, 2]

xlab_txt <- sprintf("PC1 (%.1f%%)", pc1_var)
ylab_txt <- sprintf("PC2 (%.1f%%)", pc2_var)

x <- df$PC1; y <- df$PC2
cx <- mean(range(x)); cy <- mean(range(y))
half_range <- max(diff(range(x)), diff(range(y))) / 2
xlim <- c(cx - half_range, cx + half_range)
ylim <- c(cy - half_range, cy + half_range)

cols <- c(
  "Ct"         = "#00000080",
  "LSD"         = "#E69F00",
  "Abeta"       = "#2F5DA8FF",
  "Abeta + LSD" = "#8E44ADFF"
)

df$Trt <- factor(df$Trt, levels = c("Ct", "LSD", "Abeta", "Abeta + LSD"))

p_main <- ggplot(df, aes(PC1, PC2, colour = Trt, fill = Trt)) +
  geom_point(size = 3) +   # ⬅️ rimosso alpha
  geom_xsidedensity(alpha = 0.25, linewidth = 0.4) +
  geom_ysidedensity(alpha = 0.25, linewidth = 0.4) +
  scale_colour_manual(values = cols, drop = FALSE) +
  scale_fill_manual(values = cols, drop = FALSE) +
  labs(x = xlab_txt, y = ylab_txt) +
  theme_classic(base_size = 13) +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = 0.3
  )


# -------- plot solo per estrarre legenda
p_leg <- ggplot(df, aes(PC1, PC2, colour = Trt, fill = Trt)) +
  geom_point(size = 3) +
  scale_colour_manual(values = cols, drop = FALSE) +
  scale_fill_manual(values = cols, drop = FALSE) +
  theme_void(base_size = 13) +
  theme(legend.title = element_blank())

legend <- get_legend(p_leg)

# -------- composizione finale
pdf("results/Figures/invitro_Abeta_revert_sep_PCA12_005.pdf", width = 7.5, height = 6)

plot_grid(
  p_main,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.22)
)

dev.off()


pdf("results/Figures/invitro_Abeta_revert_sep_PCA34_005.pdf", 6, 6)
ggplot(df, aes(x=PC3, y=PC4, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()

pdf("results/Figures/invitro_Abeta_revert_sep_PCA13_005.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC3, colour=Trt, fill=Trt))+geom_point(size=3)+geom_xsidedensity(alpha = .3) +
  geom_ysidedensity(alpha = .3)+theme_classic()+
  theme(
    ggside.panel.border = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.axis.line = element_blank(),
    ggside.panel.scale = .3
  ) 
dev.off()


############################
#### Number of DEGs
############################

dds_all <- DESeqDataSetFromMatrix(countData = invitro_counts,
                              colData = metadata,
                              design= ~ Replicate+Treatment)
dds_all <- DESeq(dds_all)

save(dds_all, file="Results/dds_all.RData")


ddsLSD <- results(dds_all, contrast = c("Treatment","LSD 10uM ", "Ct"))
ddsAbeta1000 <- results(dds_all, contrast = c("Treatment","Abeta 1uM ", "Ct"))
ddsAbeta500 <- results(dds_all, contrast = c("Treatment","Abeta 500nM ", "Ct"))
ddsAbeta100 <- results(dds_all, contrast = c("Treatment","Abeta 100nM", "Ct"))
ddsAbeta1000LSD <- results(dds_all, contrast = c("Treatment","Abeta 1uM +LSD ", "Ct"))
ddsAbeta500LSD <- results(dds_all, contrast = c("Treatment","Abeta 500nM + LSD ", "Ct"))
ddsAbeta100LSD <- results(dds_all, contrast = c("Treatment","Abeta 100nM + LSD", "Ct"))

length(rownames(ddsLSD)[which(ddsLSD$padj<0.05)])
n100<-length(rownames(ddsAbeta100)[which(ddsAbeta100$padj<0.05)])
n500<-length(rownames(ddsAbeta500)[which(ddsAbeta500$padj<0.05)])
n1000<-length(rownames(ddsAbeta1000)[which(ddsAbeta1000$padj<0.05)])
n100LSD<-length(rownames(ddsAbeta100LSD)[which(ddsAbeta100LSD$padj<0.05)])
n500LSD<-length(rownames(ddsAbeta500LSD)[which(ddsAbeta500LSD$padj<0.05)])
n1000LSD<-length(rownames(ddsAbeta1000LSD)[which(ddsAbeta1000LSD$padj<0.05)])

library(ggplot2)


ctrl <- c(n100, n500, n1000)
lsd  <- c(n100LSD, n500LSD, n1000LSD)

df <- data.frame(
  Condition = factor(
    c(rep("A-beta", length(ctrl)),
      rep("A-beta + LSD", length(lsd))),
    levels = c("A-beta", "A-beta + LSD")
  ),
  N_DEGs = c(ctrl, lsd)
)


wt <- wilcox.test(ctrl, lsd, paired = TRUE, alternative = "greater")
pval <- wt$p.value
p_lab <- paste0("paired Wilcoxon\np = ", signif(pval, 2))

pdf("Results/Figures/N_DEGs.pdf", width = 4.5, height = 7)

ggplot(df, aes(Condition, N_DEGs, fill = Condition)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2,
    shape = 21,
    fill = "white",
    color = "grey20"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(df$N_DEGs) * 1.08,
    label = p_lab,
    size = 4
  ) +
  scale_fill_manual(values = c("grey80", "#E69F00")) +
  ylab("Number of DEGs") +
  xlab("") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 15, 10)
  ) +
  coord_cartesian(ylim = c(NA, max(df$N_DEGs) * 1.15))

dev.off()


#######################
### PCA of log2FC matrix
########################

log2FCmat<-cbind(ddsLSD$log2FoldChange,
       ddsAbeta1000$log2FoldChange,
       ddsAbeta500$log2FoldChange,
       ddsAbeta100$log2FoldChange,
       ddsAbeta1000LSD$log2FoldChange,
       ddsAbeta500LSD$log2FoldChange,
       ddsAbeta100LSD$log2FoldChange
)


pca<-PCA(t(log2FCmat))


df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               condition=c("LSD", "Abeta 100","Abeta 500","Abeta 1000", "Abeta 100+LSD", "Abeta 500+LSD", "Abeta 1000+LSD"))

pdf("results/Figures/invitro_PCA_log2FC.pdf", 6, 6)
ggplot(df, aes(x=PC1, y=PC2, label=condition))+geom_point(size=2)+geom_label()+theme_classic()
dev.off()


##########
## GO for Abeta and Abeta + LSD
#########

DEGs_seq<-c("ddsAbeta100", "ddsAbeta500", "ddsAbeta1000", "ddsAbeta100LSD", "ddsAbeta500LSD", "ddsAbeta1000LSD")

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down<-list()
ego_up<-list()
i<-0
for(DE in DEGs_seq){
  i<-i+1
   DEGs<-get(DE)
  genes_up<-rownames(DEGs)[which(DEGs$log2FoldChange>0 & DEGs$padj<0.05)]
  genes_dn<-rownames(DEGs)[which(DEGs$log2FoldChange<0 & DEGs$padj<0.05)]
  
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

save(ego_up, file="Results/ego_up_invitro.RData")
save(ego_down, file="Results/ego_down_invitro.RData")


ego_down_Abeta<- enrichGO(gene = Abeta_genes_dn,
                         keyType="SYMBOL",
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         universe=universe, OrgDb="org.Mm.eg.db")

ego_up_Abeta<- enrichGO(gene = Abeta_genes_up,
                       keyType="SYMBOL",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       universe=universe, OrgDb="org.Mm.eg.db")

save(ego_up_Abeta, file="Results/ego_up_Abeta.RData")
save(ego_down_Abeta, file="Results/ego_down_Abeta.RData")

