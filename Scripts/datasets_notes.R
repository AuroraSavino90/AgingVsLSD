colnames(GSE5392)<-gsub(".cel.gz", "", colnames(GSE5392))
colnames(GSE71620)<-unlist(strsplit(colnames(GSE71620),"_"))[seq(1,ncol(GSE71620)*2,2)]

colnames(GSE17612)<-gsub(".CEL.gz", "", colnames(GSE17612))##riguardare data e metadata
colnames(GSE22570)<-gsub(".CEL.gz", "", colnames(GSE22570))


################################
##########GSE17612
###############################
library(GEOquery)
gds<-getGEO("GSE17612",, AnnotGPL = TRUE)

GSE17612_meta<-pData(gds$GSE17612_series_matrix.txt.gz)

metadata<-metadata[-which(metadata$Dataset=="GSE17612"),]
metadata_GSE17612<-data.frame(Sample=colnames(GSE17612), 
                              Dataset=rep("GSE17612", ncol(GSE17612)),
                              Age=as.numeric(GSE17612_meta$`age:ch1`)*365,
                              Gender=GSE17612_meta$`gender:ch1`,
                              Race=rep(NA, ncol(GSE17612)),
                              Suicide=rep(NA, ncol(GSE17612)),
                              Diagnosis=rep("Healthy", ncol(GSE17612)),
                              Alcohol=rep(NA, ncol(GSE17612)),
                              Drugs=rep(NA, ncol(GSE17612)),
                              Smoke=rep(NA, ncol(GSE17612)),
                              Medication=rep(NA, ncol(GSE17612)),
                              Organism=rep("Homo sapiens", ncol(GSE17612)),
                              Region=rep("SFG", ncol(GSE17612)),
                              Platform=GSE17612_meta$platform_id,
                              Region_simpl=rep("SFG", ncol(GSE17612)))

metadata<-rbind.data.frame(metadata, metadata_GSE17612)
###Details on PFC anatomy
GSE161986 PFC -> SFG
GSE21138 PFC -> BA46
GSE44147 -> SFG
GSE4675 -> mouse
GSE51264 -> PFC not specified
GSE106669 -> PFC not specified
GSE49376 -> BA9
GSE5390 -> DLPFC
GSE85377 -> main sulci
GSE87082 -> rat

#####
## nel confronto con LSD considerare anche PFC generica (anche quando poi è in realtà SFG, si assume che abbiano preso la parte del SFG che fa parte anche della PFC)
# nel dataset integrato non considerare SFG come DLPFC, perché questa fa parte del middle gyrus
# (NON SEREVE) cambiare le annotazioni di Region_simpl GSE21138, GSE49376 e GSE5390 in DLPFC
