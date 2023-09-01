load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging/alldata_13Apr.RData")
load("/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging/raw_data_aging.RData")

#da sistemare GSE30272
GSE113834<-changenames(data=as.matrix(GSE113834), anno=cbind(GSE113834_anno$ID, GSE113834_anno$`Gene Symbol`))

GSE13564<-changenames(data=exprs(GSE13564), anno=cbind(GSE13564_anno$ID, GSE13564_anno$`Gene symbol`))
GSE161986<-changenames(data=exprs(GSE161986), anno=cbind(GSE161986_anno$ID, GSE161986_anno$`Gene symbol`))
GSE17612<-changenames(data=exprs(GSE17612), anno=cbind(GSE161986_anno$ID, GSE161986_anno$`Gene symbol`))
GSE21138<-changenames(data=exprs(GSE21138), anno=cbind(GSE21138_anno$ID, GSE21138_anno$`Gene symbol`))
GSE21935<-changenames(data=exprs(GSE21935), anno=cbind(GSE21935_anno$ID, GSE21935_anno$`Gene symbol`))
GSE22570<-changenames(data=exprs(GSE22570), anno=cbind(GSE22570_anno$ID, GSE22570_anno$`Gene symbol`))
GSE49376<-GSE49376[rowSums(is.na(GSE49376))==0,]
GSE5388<-changenames(data=exprs(GSE5388), anno=cbind(GSE5388_anno$ID, GSE5388_anno$`Gene symbol`))
GSE5390<-changenames(data=exprs(GSE5390), anno=cbind(GSE5390_anno$ID, GSE5390_anno$`Gene symbol`))
GSE5392<-changenames(data=exprs(GSE5392), anno=cbind(GSE5392_anno$ID, GSE5392_anno$`Gene symbol`))

symbols<-gsub('([^/]*)//([^/]*)//.*', '\\2', GSE59630_anno$gene_assignment)
symbols<-gsub(' ', '', symbols)
GSE59630<-changenames(data=exprs(GSE59630), anno=cbind(GSE59630_anno$ID, symbols))

GSE60190<-changenames(data=GSE60190, anno=cbind(GSE60190_anno$ID, GSE60190_anno$`Gene symbol`))
GSE71620<-changenames(data=exprs(GSE71620), anno=cbind(GSE71620_anno$ID, GSE71620_anno$`Gene symbol`))
GSE11512<-changenames(data=exprs(GSE11512), anno=cbind(GSE5392_anno$ID, GSE5392_anno$`Gene symbol`))

GSE92538<-changenames(data=exprs(GSE92538), anno=cbind(GSE5392_anno$ID, GSE5392_anno$`Gene symbol`))
GSE92538_GPL17027<-GSE92538
rm(GSE92538)

GSE37981<-changenames(data=exprs(GSE37981), anno=cbind(GSE37981_anno$ID, GSE37981_anno$`Gene symbol`))

####
add_meta<-data.frame(matrix(nrow=ncol(GSE37981), ncol=length(categories)))
colnames(add_meta)<-categories

add_meta$Dataset<-"GSE37981"
add_meta$Sample<-GSE37981_meta$geo_accession
add_meta$Age<-as.numeric(GSE37981_meta$`age:ch1`)*356
add_meta$Gender<-GSE37981_meta$`gender:ch1`
add_meta$Region<-"STG"
add_meta$Platform<-GSE37981_meta$platform_id
add_meta$Organism<-GSE37981_meta$organism_ch1
add_meta$Region_simpl<-NA

metadata<-rbind.data.frame(metadata, add_meta )


