####Download CEL files
#FIX GSE30272 from 27Apr
library(GEOquery)
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
#Due dataset da sistemare
# expression
proj_path<-"/Users/aurora.savino/Library/CloudStorage/OneDrive-Htechnopole/Documents/Work/Projects/Metanalyses/Aging"

getGEOSuppFiles("GSE37981")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE37981/GSE37981_RAW/", sep=""))
GSE37981 <- affy::rma(Data)

getGEOSuppFiles("GSE113834")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE113834/GSE113834_RAW/", sep=""))
GSE113834 <- rma(Data)

getGEOSuppFiles("GSE11512")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE11512/GSE11512_RAW/", sep=""))
GSE11512 <- rma(Data)

getGEOSuppFiles("GSE125681")
library(illuminaio)
idat<-list.files(path = paste(proj_path, "/GSE125681_RAW/", sep=""), pattern = "idat")
Data <- read.idat(idat, path = paste(proj_path, "/GSE125681_RAW", sep=""), bgxfile =  "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz")
GSE125681<-neqc(Data)

getGEOSuppFiles("GSE13564")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE13564_RAW/", sep=""))
GSE13564 <- rma(Data)

getGEOSuppFiles("GSE161986")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE161986/GSE161986_RAW/", sep=""))
GSE161986 <- rma(Data)

getGEOSuppFiles("GSE17612")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE17612_RAW/", sep=""))
GSE17612 <- rma(Data)

getGEOSuppFiles("GSE21138")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE21138_RAW/", sep=""))
GSE21138 <- rma(Data)

getGEOSuppFiles("GSE21935")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE21935/GSE21935_RAW/", sep=""))
GSE21935 <- rma(Data)

getGEOSuppFiles("GSE22570")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE22570_RAW/", sep=""))
GSE22570 <- rma(Data)

getGEOSuppFiles("GSE30272")#txt files
files<-list.files(path = paste(proj_path, "/GSE30272_RAW/", sep=""), pattern = "txt")
aa<-read.csv(paste(proj_path, "/GSE30272_RAW/",idat[1], sep=""), sep="\t")
GSE30272<-matrix(nrow=nrow(aa), ncol=length(files))
rownames(GSE30272)<-aa$PlatePos
colnames(GSE30272)<-files
for(f in files){
  tmp<-read.csv(paste(proj_path, "/GSE30272_RAW/",f, sep=""), sep="\t")
  GSE30272[,f]<-tmp$calRatio
}
GSE30272<-oligo::rma(GSE30272)
GSE30272<-changenames(data=GSE30272, anno=cbind(GSE30272_anno$ID, GSE30272_anno$Gene_Symbol))
colnames(GSE30272)<-GSE30272_meta[match(colnames(GSE30272), GSE30272_meta[,1]),2]

# create an instance of ExpressionSet
ExpressionSet()
prova<-ExpressionSet(assayData=GSE30272_new[1:2,])
prova<-affy::rma(prova)
exprs(Data)<-GSE30272_new
Data$sample<-colnames(GSE30272_new)
prova<-affy::rma(Data)

getGEOSuppFiles("GSE49376")
Data<-read.csv(paste(proj_path, "/GSE49376/GSE49376_non-normalized.txt", sep=""), sep="\t", row.names = 1)
GSE49376<-neqc(Data[,seq(1,95,2)], detection.p = Data[,seq(2,96,2)])

getGEOSuppFiles("GSE5388")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE5388/GSE5388_RAW/", sep=""))
GSE5388 <- rma(Data)

getGEOSuppFiles("GSE5390")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE5390/GSE5390_RAW/", sep=""))
GSE5390 <- rma(Data)

getGEOSuppFiles("GSE5392")
Data <- ReadAffy(celfile.path = paste(proj_path, "/GSE5392_RAW/", sep=""))
GSE5392 <- rma(Data)

getGEOSuppFiles("GSE59630")
library(oligo)
celfiles <- list.files(path=paste(proj_path, "/GSE59630_RAW/", sep=""), full = TRUE)
rawData <- read.celfiles(celfiles)
GSE59630 <- rma(rawData)

getGEOSuppFiles("GSE60190")
Data<-read.csv(paste(proj_path, "/GSE60190/GSE60190_non-normalized.txt", sep=""), sep="\t", row.names = 1)
GSE60190<-neqc(Data[,seq(1,266,2)], detection.p = Data[,seq(2,266,2)])

getGEOSuppFiles("GSE71620")
celfiles <- list.files(path=paste(proj_path, "/GSE71620_RAW/", sep=""), full = TRUE)
rawData <- read.celfiles(celfiles)
GSE71620<- rma(rawData)

getGEOSuppFiles("GSE92538")
celfiles <- list.files(path=paste(proj_path, "/GSE92538_RAW/", sep=""), pattern="cel")
celfiles<-celfiles[grep("_133A_", celfiles)]
rawData <- read.celfiles(paste(proj_path, "/GSE92538_RAW/", celfiles, sep=""))
GSE92538<- rma(rawData)

