rm(GSE84422_1)#no PFC

GSE84422<-rbind(GSE84422_2, GSE84422_3)
rm(GSE84422_2)
rm(GSE84422_3)
rm(GSE64810)
#################################
### Metadata categories
##################################


categories<-c("Sample","Dataset","Age", "Gender","Race","Diagnosis",
              "Medication","Region", "Platform")
#age in days

#lista con tutti i dataset
datasetsNames<-setdiff(objects()[grep("GSE", objects())],objects()[grep("_", objects())])
#da aggiungere molti dataset con _

##per avere il numero totale di righe, riempire la colonna dataset
#DA FARE QUANDO SARANNO STATI AGGIUSTATI TUTTI I DATASET CON REPLICHE TECNICHE etc (da alldata_metadata)

datasetGEO<-c()
for(i in 1:length(datasetsNames)){
  datasetGEO<-c(datasetGEO, rep(datasetsNames[i], ncol(get(datasetsNames[i]))))
}

#creare una tabella metadata con tutti "NA" in cui poi inserire i metadati

metadata<-data.frame(matrix(nrow=length(datasetGEO), ncol=length(categories)))
colnames(metadata)<-categories

metadata$Dataset<-datasetGEO

#la colonna sample ? data dal nome delle colonne di ogni dataset
sampleGEO<-c()
for(i in 1:length(datasetsNames)){
  sampleGEO<-c(sampleGEO, colnames(get(datasetsNames[i])))
}

metadata$Sample<-sampleGEO
#in alcuni casi non sono corretti perch? non corrispondono ai GSE
#i dataset in cui questo avviene sono:
unique(metadata$Dataset[-grep("GSM", metadata$Sample)])
#per ognuno controllare il match tra le colonne e i metadati e cambiare sia Sample che colnames




################################
##########GSE84422_2
###############################
#neurofibrillary score
#plaque
#dementia rating
#subject id

metadata$Age[match(GSE84422_meta_2[,2], metadata$Sample)]<-as.numeric(GSE84422_meta_2$`age:ch1`)*356
metadata$Gender[match(GSE84422_meta_2[,2], metadata$Sample)]<-GSE84422_meta_2$`Sex:ch1`
metadata$Region[match(GSE84422_meta_2[,2], metadata$Sample)]<-GSE84422_meta_2$`brain region:ch1`
metadata$Race[match(GSE84422_meta_2[,2], metadata$Sample)]<-GSE84422_meta_2$`race:ch1`
metadata$Diagnosis[match(GSE84422_meta_2[,2], metadata$Sample)]<-GSE84422_meta_2$`neuropathological category:ch1`




################################
##########GSE3300
###############################

metadata$Age[match(GSE33000_meta[,2], metadata$Sample)]<-as.numeric(gsub(" yrs", "", GSE33000_meta$`age:ch2`))*356
metadata$Gender[match(GSE33000_meta[,2], metadata$Sample)]<-GSE33000_meta$`gender:ch2`
metadata$Region[match(GSE33000_meta[,2], metadata$Sample)]<-GSE33000_meta$`tissue:ch1`
metadata$Diagnosis[match(GSE33000_meta[,2], metadata$Sample)]<-GSE33000_meta$`disease status:ch2`

################################
##########GSE44770
###############################

metadata$Age[match(GSE44770_meta[,2], metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE44770_meta$characteristics_ch2.1))*356
metadata$Gender[match(GSE44770_meta[,2], metadata$Sample)]<-gsub("gender: ", "", GSE44770_meta$characteristics_ch2.2)
metadata$Region[match(GSE44770_meta[,2], metadata$Sample)]<-"DLPFC"
metadata$Diagnosis[match(GSE44770_meta[,2], metadata$Sample)]<-gsub("disease: ", "", GSE44770_meta$characteristics_ch2)


################################
##########GSE150696
###############################
#dementia with Lewy bodies (DLB) and Parkinson's disease dementia (PDD)

metadata$Age[match(GSE150696_meta[,2], metadata$Sample)]<-as.numeric(GSE150696_meta$`age:ch1`)*365
metadata$Gender[match(GSE150696_meta[,2], metadata$Sample)]<-GSE150696_meta$`Sex:ch1`
metadata$Region[match(GSE150696_meta[,2], metadata$Sample)]<-"prefrontal cortex BA9"
metadata$Diagnosis[match(GSE150696_meta[,2], metadata$Sample)]<-GSE150696_meta$description

################################
##########GSE193391 #only 700 genes
###############################
#Apoe genotype

metadata$Age[match(GSE193391_meta[,2], metadata$Sample)]<-as.numeric(GSE193391_meta$`age:ch1`)*365
metadata$Gender[match(GSE193391_meta[,2], metadata$Sample)]<-GSE193391_meta$`Sex:ch1`
metadata$Region[match(GSE193391_meta[,2], metadata$Sample)]<-GSE193391_meta$`tissue:ch1`
metadata$Diagnosis[match(GSE193391_meta[,2], metadata$Sample)]<-GSE193391_meta$`diagnosis:ch1`

metadata<-metadata[metadata$Dataset!="GSE193391",]

############
### Al controls named healthy
#############

metadata$Diagnosis[metadata$Diagnosis %in% c("Control", "CTRL", "N", "non-demented", "Normal")]<-"Healthy"
##############################
### Disease
############################
metadata$Diagnosis[metadata$Diagnosis %in% c("A", "AD", "Alzheimer's disease", "definite AD")]<-"Alzheimer"
metadata$Diagnosis[metadata$Diagnosis %in% c("Huntington's disease")]<-"Huntington"

##############################
### Gender
############################
metadata$Gender[metadata$Gender %in% c("F", "f", "Female")]<-"F"
metadata$Gender[metadata$Gender %in% c("M", "m", "Male")]<-"M"

##############################
### Region
############################
metadata$Region[metadata$Region %in% c( "prefrontal cortex brain", "Prefrontal Cortex")]<-"PFC"
metadata$Region[metadata$Region %in% c("prefrontal cortex BA9","Dorsolateral Prefrontal Cortex","DLPFC")]<-"DLPFC"


