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
load("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD/RData/Alldata_27May.RData")

load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging/GTEx_RPM_PFC_filt2.RData")
E_PFC<-E
E_PFC<-E_PFC[rowSums(E_PFC>0)>=10,]
metadata<-metadata[-which(metadata$Sample %in% setdiff(metadata$Sample[metadata$Dataset=="E_PFC"], colnames(E_PFC))),]

metadata$Dataset[metadata$Dataset=="E_PFC"]<-"GTEx_PFC"
GTEx_PFC<-E_PFC
rm(E)
rm(E_PFC)


#####mi restringo a Homo sapiens
metadata<-metadata[metadata$Organism=="Homo sapiens",]
#remove GSE5392, superseries of GSE5388
metadata<-metadata[-which(metadata$Dataset=="GSE5392"),]
metadata$Region_simpl[which(metadata$Region=="STG")]<-"TC"
metadata$Diagnosis[which(metadata$Region=="STG")]<-"Healthy"
metadata$Gender[which(metadata$Gender=="female")]<-"F"
metadata$Gender[which(metadata$Gender=="male")]<-"M"



#numero di campioni con età superiore a 18
dat<-na.omit(unique(metadata$Dataset[metadata$Organism=="Homo sapiens"]))
adult_num<-c()
for(dd in dat){
  age<-round(as.numeric(metadata$Age[which(metadata$Dataset==dd & metadata$Organism=="Homo sapiens")]))
  adult_num<-c(adult_num, length(which(age>6570)))
}

dataset_adult<-dat[which(adult_num>=10)]

metadata<-metadata[which(metadata$Dataset %in% dataset_adult),]

metadata$Smoke[metadata$Smoke %in% c("N")]<-"No"
metadata$Smoke[metadata$Smoke %in% c("Y")]<-"Yes"
metadata$Smoke[metadata$Smoke %in% c("U")]<-NA

metadata<-metadata[metadata$Dataset!="GSE30272",] #no raw data

save.image("~/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/AgingVSLSD/RData/Alldata_7Jul.RData")

