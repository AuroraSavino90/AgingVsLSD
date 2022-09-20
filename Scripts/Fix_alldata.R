#Fix RData
load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging/Aging_alldata27Apr.RData")

colnames(GSE161986)<-unlist(strsplit( colnames(GSE161986), "_"))[seq(1,ncol(GSE161986)*2,2)]
colnames(GSE113834)<-unlist(strsplit( colnames(GSE113834), "_"))[seq(1,ncol(GSE113834)*5,5)]
colnames(GSE11512)<-gsub("[.]CEL[.]gz", "", colnames(GSE11512))
colnames(GSE13564)<-gsub("[.]CEL[.]gz", "", colnames(GSE13564))
colnames(GSE5388)<-gsub("[.]cel[.]gz", "", colnames(GSE5388))
colnames(GSE59630)<-unlist(strsplit( colnames(GSE59630), "_"))[seq(1,ncol(GSE59630)*2,2)]
colnames(GSE60190)<-metadata$Sample[which(metadata$Dataset=="GSE60190")]#??????? CHECK!!!
colnames(GSE21138)<-gsub("[.]CEL[.]gz", "", colnames(GSE21138))
colnames(GSE49376)<-metadata$Sample[which(metadata$Dataset=="GSE49376")]#??????? CHECK!!!
colnames(GSE5390)<-gsub("[.]cel[.]gz", "", colnames(GSE5390))
colnames(GSE92538_GPL17027)<-GSE92538_GPL17027_meta[1,-1]

E_PFC<-log2(E_PFC+1)

rm(list=ls()[-which(ls() %in% c(unique(metadata$Dataset), "metadata"))])

save.image("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD/RData/Alldata_27May.RData")

#################
rm(list=ls())
load("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD_local/RData/Alldata_27May.RData")

load("Results/RData/GTEx_RPM_PFC_filt.RData")
metadata<-metadata[-which(metadata$Sample %in% setdiff(metadata$Sample[metadata$Dataset=="E_PFC"], colnames(GTEx_PFC))),]

metadata$Dataset[metadata$Dataset=="E_PFC"]<-"GTEx_PFC"
rm(E_PFC)


#####mi restringo a Homo sapiens
metadata<-metadata[metadata$Organism=="Homo sapiens",]
#remove GSE5392, superseries of GSE5388
metadata<-metadata[-which(metadata$Dataset=="GSE5392"),]
metadata$Region_simpl[which(metadata$Region=="STG")]<-"TC"
metadata$Diagnosis[which(metadata$Region=="STG")]<-"Healthy"
metadata$Gender[which(metadata$Gender=="female")]<-"F"
metadata$Gender[which(metadata$Gender=="male")]<-"M"



#numero di campioni con età superiore a 20
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens"]))
adult_num<-c()
for(dd in dat){
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens")]))
  adult_num<-c(adult_num, length(which(age>=7300)))
}

dataset_adult<-dat[which(adult_num>=10)]

metadata<-metadata[which(metadata$Dataset %in% dataset_adult),]

metadata$Smoke[metadata$Smoke %in% c("N")]<-"No"
metadata$Smoke[metadata$Smoke %in% c("Y")]<-"Yes"
metadata$Smoke[metadata$Smoke %in% c("U")]<-NA

##include GSE33000
################################
##########GSE33000 
###############################
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

library(GEOquery)
gds<-getGEO("GSE33000",destdir="Data/", AnnotGPL = TRUE)

GSE33000<-exprs(gds$GSE33000_series_matrix.txt.gz)
GSE33000_meta<-pData(gds$GSE33000_series_matrix.txt.gz)
GSE33000_anno<-fData(gds$GSE33000_series_matrix.txt.gz)

GSE33000<-changenames(data=GSE33000, anno=cbind(as.character(GSE33000_anno$ID), GSE33000_anno$`Gene symbol`))

metadata_GSE33000<-data.frame(Sample=colnames(GSE33000),
                              Dataset="GSE33000",
                              Age=as.numeric(gsub(" yrs", "", GSE33000_meta$`age:ch2`))*356,
                              Gender=GSE33000_meta$`gender:ch2`,
                              Race=NA,
                              Suicide=NA,
                              Diagnosis=GSE33000_meta$`disease status:ch2`,
                              Alcohol=NA,
                              Drugs=NA,
                              Smoke=NA,
                              Medication=NA,
                              Organism="Homo sapiens", 
                              Region=GSE33000_meta$`tissue:ch1`,
                              Platform=GSE33000_meta$platform_id,
                              Region_simpl="PFC")

metadata<-rbind.data.frame(metadata, metadata_GSE33000)

metadata$Gender[metadata$Gender=="female"]<-"F"
metadata$Gender[metadata$Gender=="male"]<-"M"
metadata$Diagnosis[metadata$Diagnosis=="Alzheimer's disease"]<-"AD"
metadata$Diagnosis[metadata$Diagnosis=="Huntington's disease"]<-"HD"
metadata$Diagnosis[metadata$Diagnosis=="non-demented"]<-"Healthy"

save.image("Results/RData/Alldata_20Sep.RData")

