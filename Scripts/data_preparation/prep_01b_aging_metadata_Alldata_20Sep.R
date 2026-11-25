#####################################################################
## Data preparation, step 1b: sample metadata and final assembly of
## Data/Alldata_20Sep.RData
#####################################################################
##
## Inputs
##   EXPR_FILE   output of prep_01a_aging_expression.R
##   GEO_DIR     GEO series matrices downloaded in 2022 (sample annotations)
##   GTEX_COUNTS GTEx v8 gene read counts, frontal cortex BA9 (Data/); GTEx_PFC
##               is built from it as in Scripts/00_GTEx_filt.R, with the HGNC
##               symbols of Ensembl 107
##   GTEX_SUBJ   GTEx v8 subject phenotypes (Data/)
##   GTEX_SAMP   GTEx v8 sample attributes (GTEx portal)
##   GSE33000 and GPL4372 are read from DOWNLOAD_DIR (Data/raw_data_notclosed/dementia_GEO).
## Output: OUT_FILE with the expression matrices and `metadata`.
##
## Notes reproduced from the original file
##   - Ages are in days (some datasets were converted with 356 instead of
##     365 days per year).
##   - Column names of GSE17612, GSE21935, GSE22570, GSE37981, GSE71620 and
##     GSE125681 do not match metadata$Sample. None of these datasets has
##     DLPFC/PFC healthy samples, so they are never selected in the analyses.
#####################################################################

# Run from the repository root:
#   Rscript Scripts/data_preparation/prep_01b_aging_metadata_Alldata_20Sep.R
# Settings can be overridden with environment variables of the same name.
EXPR_FILE    <- Sys.getenv("EXPR_FILE", "Data/regenerated/aging_expression.RData")
GEO_DIR      <- Sys.getenv("GEO_DIR", "Data/raw_data_notclosed/aging_GEO")
# GTEx v8 (release 2017-06-05), open-access files of the GTEx portal:
# https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/counts-by-tissue/gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz
# https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
# https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
GTEX_COUNTS  <- Sys.getenv("GTEX_COUNTS", "Data/raw_data_notclosed/gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz")
GTEX_SUBJ    <- Sys.getenv("GTEX_SUBJ", "Data/raw_data_notclosed/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
GTEX_SAMP    <- Sys.getenv("GTEX_SAMP", "Data/raw_data_notclosed/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
GTEX_RELEASE <- as.integer(Sys.getenv("GTEX_RELEASE", "107"))   # Ensembl release for the GTEx gene symbols
ENSEMBL_CACHE <- Sys.getenv("ENSEMBL_CACHE", "Data/ensembl")
# Folder with GSE33000_series_matrix.txt.gz and GPL4372.annot.gz
DOWNLOAD_DIR <- Sys.getenv("DOWNLOAD_DIR", "Data/raw_data_notclosed/dementia_GEO")
OUT_FILE     <- Sys.getenv("OUT_FILE", "Data/regenerated/Alldata_20Sep.RData")

suppressMessages({
  library(GEOquery)
  library(data.table)
})
load(EXPR_FILE)

g <- function(...) file.path(GEO_DIR, ...)
series_meta <- function(file) pData(getGEO(filename = g(file), getGPL = FALSE))

#####################################################################
## 1. GEO sample annotations
#####################################################################
GSE71620_meta    <- series_meta("GSE71620_series_matrix.txt.gz")
GSE13564_meta    <- series_meta("GSE13564_series_matrix.txt.gz")
GSE102556_H_meta <- series_meta("GSE102556-GPL11154_series_matrix.txt.gz")
GSE30272_meta    <- series_meta("GSE30272_series_matrix.txt.gz")
GSE60190_meta    <- series_meta("GSE60190_series_matrix.txt.gz")
GSE22570_meta    <- series_meta("GSE22570_series_matrix.txt.gz")
GSE59630_meta    <- series_meta("GSE59630_series_matrix.txt.gz")
GSE101521_meta   <- series_meta("GSE101521-GPL16791_series_matrix.txt.gz")
GSE11512_meta    <- series_meta("GSE11512_series_matrix.txt.gz")
GSE59288_H_meta  <- series_meta("GSE59288-GPL11154_series_matrix.txt.gz")
GSE64810_meta    <- series_meta("GSE64810_series_matrix.txt.gz")
GSE5388_meta     <- series_meta("GSE5388_series_matrix.txt.gz")
GSE161986_meta   <- series_meta("GSE161986_series_matrix.txt.gz")
GSE21138_meta    <- series_meta("GSE21138_series_matrix.txt.gz")
GSE17612_meta    <- series_meta("GSE17757_series_matrix.txt.gz") 
GSE17757_meta    <- GSE17612_meta
GSE5392_meta     <- series_meta("GSE5392_series_matrix.txt.gz")
GSE49376_meta    <- series_meta("GSE49376_series_matrix.txt.gz")
GSE21935_meta    <- series_meta("GSE21935_series_matrix.txt.gz")
GSE125681_meta   <- series_meta("GSE125681_series_matrix.txt.gz")
GSE5390_meta     <- series_meta("GSE5390_series_matrix.txt.gz")
GSE102741_meta   <- series_meta("GSE102741_series_matrix.txt.gz")
GSE174409_meta   <- series_meta("GSE174409_series_matrix.txt.gz")
GSE37981_meta    <- series_meta("GSE37981_series_matrix.txt.gz")
GSE99349_meta    <- series_meta("GSE99349_series_matrix.txt.gz")
GSE99349_meta    <- GSE99349_meta[-which(GSE99349_meta$title == "Neun DLPFC_Total_RNAseq_28"),]
# GSE92538 (GPL17027): header lines of the series matrix as a table
GSE92538_GPL17027_meta <- read.csv(g("GSE92538-GPL17027_series_matrix_meta.txt"), sep = "\t", header = T)

#####################################################################
## 2. GTEx frontal cortex BA9 (00_GTEx_filt.R, Aging_PFC.R)
#####################################################################
# GTEx_PFC = Results/RData/GTEx_RPM_PFC_filt.RData of 2022: log2(RPM + 1),
# without the 10% of samples with the highest mitochondrial read counts and
# without genes expressed in fewer than 10 samples (same code as 00_GTEx_filt.R)
source("Scripts/utilities.R")                        
source("Scripts/data_preparation/ensembl_release_map.R")
GTEx_PFC <- fread(GTEX_COUNTS, data.table = FALSE)
rownames(GTEx_PFC) <- GTEx_PFC[, 2]
genes_conv <- ensembl_release_map(GTEX_RELEASE, "hsapiens", ENSEMBL_CACHE)
symbol <- genes_conv[match(gsub("\\..*", "", rownames(GTEx_PFC)), genes_conv[, 1]), 2]
GTEx_PFC <- changenames(GTEx_PFC[, -c(1:3)], anno = cbind(rownames(GTEx_PFC), symbol))
GTEx_PFC <- log2(t(t(GTEx_PFC)/colSums(GTEx_PFC))*1000000 + 1)
MT_genes <- colSums(2^GTEx_PFC[grep("MT-", rownames(GTEx_PFC)), ] - 1)
GTEx_PFC <- GTEx_PFC[, which(MT_genes < quantile(MT_genes, 0.9))]
GTEx_PFC <- GTEx_PFC[rowSums(GTEx_PFC > 0) >= 10, ]
rm(genes_conv, symbol, MT_genes)

subj <- fread(GTEX_SUBJ)
samp <- fread(GTEX_SAMP)
samp[, SUBJID := gsub('([^-]*)-([^-]*)-.*', '\\1-\\2', SAMPID)]
PFC_samp <- samp$SAMPID[samp$SMTSD == "Brain - Frontal Cortex (BA9)"]
PFC_subj <- samp$SUBJID[samp$SMTSD == "Brain - Frontal Cortex (BA9)"]
PFC_subj <- PFC_subj[which(PFC_samp %in% colnames(GTEx_PFC))]
PFC_samp <- PFC_samp[which(PFC_samp %in% colnames(GTEx_PFC))]
age_bin <- subj$AGE[match(PFC_subj, subj$SUBJID)]
age_bin[age_bin == "20-29"] <- 1
age_bin[age_bin == "30-39"] <- 2
age_bin[age_bin == "40-49"] <- 3
age_bin[age_bin == "50-59"] <- 4
age_bin[age_bin == "60-69"] <- 5
age_bin[age_bin == "70-79"] <- 6

#####################################################################
## 3. Metadata table
## One row per sample; Sample = sample names of each dataset as in the GEO
## series matrices (GSM IDs) or in the processed files.
#####################################################################
categories <- c("Sample", "Dataset", "Age", "Gender", "Race", "Suicide", "Diagnosis",
                "Alcohol", "Drugs", "Smoke", "Medication", "Organism", "Region", "Platform")
samples <- list(
  GSE101521 = colnames(GSE101521), GSE102741 = colnames(GSE102741),
  GSE11512 = GSE11512_meta[,2], GSE125681 = GSE125681_meta[,2], GSE13564 = GSE13564_meta[,2],
  GSE161986 = GSE161986_meta[,2], GSE174409 = colnames(GSE174409), GSE17612 = GSE17612_meta[,2],
  GSE21138 = GSE21138_meta[,2], GSE21935 = GSE21935_meta[,2], GSE22570 = GSE22570_meta[,2],
  GSE30272 = GSE30272_meta[,2], GSE49376 = GSE49376_meta[,2], GSE5388 = GSE5388_meta[,2],
  GSE5390 = GSE5390_meta[,2], GSE5392 = GSE5392_meta[,2], GSE59630 = GSE59630_meta[,2], GSE60190 = GSE60190_meta[,2],
  GSE64810 = colnames(GSE64810), GSE71620 = GSE71620_meta[,2], GSE99349 = colnames(GSE99349),
  GSE92538_GPL17027 = unlist(GSE92538_GPL17027_meta[1,-1]), GSE102556_H = colnames(GSE102556_H),
  E_PFC = PFC_samp, GSE59288_H = colnames(GSE59288_H))
metadata <- data.frame(matrix(nrow = length(unlist(samples)), ncol = length(categories)))
colnames(metadata) <- categories
metadata$Sample  <- unname(unlist(samples))
metadata$Dataset <- rep(names(samples), lengths(samples))

## GSE71620
metadata$Age[match(GSE71620_meta[,2], metadata$Sample)]<-as.numeric(GSE71620_meta$`age:ch1`)*356
metadata$Gender[match(GSE71620_meta[,2], metadata$Sample)]<-GSE71620_meta$`Sex:ch1`
metadata$Region[match(GSE71620_meta[,2], metadata$Sample)]<-GSE71620_meta$source_name_ch1
metadata$Race[match(GSE71620_meta[,2], metadata$Sample)]<-GSE71620_meta$`race:ch1`
metadata$Platform[match(GSE71620_meta[,2], metadata$Sample)]<-GSE71620_meta$platform_id
metadata$Organism[match(GSE71620_meta[,2], metadata$Sample)]<-GSE71620_meta$organism_ch1

## GSE92538_GPL17027
metadata$Age[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE92538_GPL17027_meta[17,-1]))*365
metadata$Region[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-(gsub("tissue: ", "", GSE92538_GPL17027_meta[20,-1]))
metadata$Race[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-(gsub("race: ", "", GSE92538_GPL17027_meta[16,-1]))
metadata$Gender[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-(gsub("gender: ", "", GSE92538_GPL17027_meta[15,-1]))
metadata$Diagnosis[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-(gsub("diagnosis: ", "", GSE92538_GPL17027_meta[11,-1]))
metadata$Suicide[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-(gsub("suicide (1=yes): ", "", GSE92538_GPL17027_meta[19,-1]))
metadata$Platform[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-( GSE92538_GPL17027_meta[34,-1])
metadata$Organism[match(GSE92538_GPL17027_meta[1,-1], metadata$Sample)]<-"Homo sapiens"

## GSE17757 (same samples as GSE17612)
metadata$Age[match(GSE17757_meta[,2], metadata$Sample)]<-gsub(" days", "", GSE17757_meta$`age:ch1`)
metadata$Region[match(GSE17757_meta[,2], metadata$Sample)]<-GSE17757_meta$`tissue:ch1`
metadata$Gender[match(GSE17757_meta[,2], metadata$Sample)]<-(GSE17757_meta$`sex:ch1`)
metadata$Platform[match(GSE17757_meta[,2], metadata$Sample)]<-(GSE17757_meta$platform_id)
metadata$Organism[match(GSE17757_meta[,2], metadata$Sample)]<-GSE17757_meta$organism_ch1

## GSE13564
metadata$Age[match(GSE13564_meta[,2], metadata$Sample)]<-as.numeric(gsub("DLPFC age in years | 1| 2", "", GSE13564_meta$source_name_ch1))*365
metadata$Region[match(GSE13564_meta[,2], metadata$Sample)]<-"DLPFC"
metadata$Organism[match(GSE13564_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Gender[match(GSE13564_meta[,2], metadata$Sample)]<-unlist(lapply(GSE13564_meta$characteristics_ch1, function(x){substr(x, nchar(x), nchar(x))}))

## GTEx (age = centre of the 10-year age bin)
metadata$Age[match(PFC_samp, metadata$Sample)]<-(as.numeric(age_bin)+1.5)*10*365
metadata$Gender[match(PFC_samp, metadata$Sample)]<-subj$SEX[match(PFC_subj, subj$SUBJID)]
metadata$Organism[match(PFC_samp, metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(PFC_samp, metadata$Sample)]<-"PFC"

## GSE102556_H
ind<-gsub(".AntIns|.BA11|.BA25|.BA8_9|.Nac|.Subic|X", "", colnames(GSE102556_H))
ind_meta<-gsub(": Anterior Insula \\(aINS\\)|: Nucleus Accumbens \\(Nac\\)|: Subiculum \\(Sub\\)|: Dorsolateral prefrontal cortex \\(dlPFC; BA8\\/9\\)|: Orbitofrontal \\(OFC; BA11\\)|: Cingulate gyrus 25 \\(Cg25\\)", "", GSE102556_H_meta[,1])
metadata$Age[match(colnames(GSE102556_H), metadata$Sample)]<-(as.numeric(GSE102556_H_meta$`age:ch1`)*365)[match(ind, ind_meta)]
metadata$Region[match(colnames(GSE102556_H), metadata$Sample)]<-c(rep("Anterior Insula",48), rep("OFC",48),rep("CG",46 ), rep("dlPFC",48), rep("Nucleus Accumbens",48), rep("Subiculum",43))
metadata$Gender[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("gender: ", "", GSE102556_H_meta$characteristics_ch1.1)[match(ind, ind_meta)]
metadata$Suicide[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("Cause of death: ", "", GSE102556_H_meta$characteristics_ch1.2)[match(ind, ind_meta)]
metadata$Diagnosis[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("phenotype: ", "", GSE102556_H_meta$characteristics_ch1.6)[match(ind, ind_meta)]
metadata$Alcohol[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("alcohol: ", "", GSE102556_H_meta$characteristics_ch1.7)[match(ind, ind_meta)]
metadata$Drugs[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("drug type: ", "", GSE102556_H_meta$characteristics_ch1.9)[match(ind, ind_meta)]
metadata$Medication[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("medication type: ", "", GSE102556_H_meta$characteristics_ch1.11)[match(ind, ind_meta)]
metadata$Smoke[match(colnames(GSE102556_H), metadata$Sample)]<-gsub("smoking: ", "", GSE102556_H_meta$characteristics_ch1.12)[match(ind, ind_meta)]
metadata$Platform[match(colnames(GSE102556_H), metadata$Sample)]<-GSE102556_H_meta$platform_id
metadata$Organism[match(colnames(GSE102556_H), metadata$Sample)]<-"Homo sapiens"

## GSE30272
metadata$Age[match(GSE30272_meta[,2], metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE30272_meta$characteristics_ch1.1))*365
metadata$Gender[match(GSE30272_meta[,2], metadata$Sample)]<-gsub("Sex: ", "", GSE30272_meta$characteristics_ch1.2)
metadata$Smoke[match(GSE30272_meta[,2], metadata$Sample)]<-gsub("smoke at death: ", "", GSE30272_meta$characteristics_ch1.7)
metadata$Race[match(GSE30272_meta[,2], metadata$Sample)]<-gsub("race: ", "", GSE30272_meta$characteristics_ch1.3)
metadata$Organism[match(GSE30272_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE30272_meta[,2], metadata$Sample)]<-"dorsolateral prefrontal cortex (DLPFC BA 46/9)"

## GSE60190
metadata$Age[match(GSE60190_meta[,2], metadata$Sample)]<-as.numeric(GSE60190_meta$`age:ch1`)*365
metadata$Gender[match(GSE60190_meta[,2], metadata$Sample)]<-gsub("Sex: ", "", GSE60190_meta$characteristics_ch1.7)
metadata$Race[match(GSE60190_meta[,2], metadata$Sample)]<-GSE60190_meta$`race:ch1`
metadata$Diagnosis[match(GSE60190_meta[,2], metadata$Sample)]<-GSE60190_meta$`dx:ch1`
metadata$Region[match(GSE60190_meta[,2], metadata$Sample)]<-"DLPFC"
metadata$Organism[match(GSE60190_meta[,2], metadata$Sample)]<-"Homo sapiens"

## GSE22570
metadata$Age[match(GSE22570_meta[,2], metadata$Sample)]<-as.numeric(gsub(" days|\\(prenatal individual; predicted to be ~30 days before birth\\)|\\(prenatal individual; predicted to be ~15 days before birth\\)", "", GSE22570_meta$`age:ch1`))
metadata$Gender[match(GSE22570_meta[,2], metadata$Sample)]<-GSE22570_meta$`gender:ch1`
metadata$Region[match(GSE22570_meta[,2], metadata$Sample)]<-GSE22570_meta$`tissue:ch1`
metadata$Organism[match(GSE22570_meta[,2], metadata$Sample)]<-GSE22570_meta$organism_ch1

## GSE59630
age_GSE59630<-GSE59630_meta$`age:ch1`
age_GSE59630[grep("yr", age_GSE59630)]<-as.numeric(gsub("yr", "", age_GSE59630[grep("yr", age_GSE59630)]))*365
age_GSE59630[grep("mo", age_GSE59630)]<-as.numeric(gsub("mo", "", age_GSE59630[grep("mo", age_GSE59630)]))*30
age_GSE59630[grep("wg", age_GSE59630)]<-(as.numeric(gsub("wg", "", age_GSE59630[grep("wg", age_GSE59630)]))-36)*7
metadata$Age[match(GSE59630_meta[,2], metadata$Sample)]<-as.numeric(age_GSE59630)
metadata$Gender[match(GSE59630_meta[,2], metadata$Sample)]<-gsub("Sex: ", "", GSE59630_meta[,13])
metadata$Region[match(GSE59630_meta[,2], metadata$Sample)]<-gsub("region: ", "", GSE59630_meta[,11])
metadata$Race[match(GSE59630_meta[,2], metadata$Sample)]<-GSE59630_meta$`race:ch1`
metadata$Diagnosis[match(GSE59630_meta[,2], metadata$Sample)]<-GSE59630_meta$`disease status:ch1`
metadata$Organism[match(GSE59630_meta[,2], metadata$Sample)]<-"Homo sapiens"

## GSE101521
metadata$Age[match(GSE101521_meta[,1], metadata$Sample)]<-as.numeric(GSE101521_meta$`age (yrs):ch1`)*365
metadata$Diagnosis[match(GSE101521_meta[,1], metadata$Sample)]<-GSE101521_meta$`diagnosis:ch1`
metadata$Gender[match(GSE101521_meta[,1], metadata$Sample)]<-gsub("Sex: ", "", GSE101521_meta$characteristics_ch1.2)
metadata$Organism[match(GSE101521_meta[,1], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE101521_meta[,1], metadata$Sample)]<-GSE101521_meta$source_name_ch1

## GSE5392 (superseries of GSE5388 and GSE5390)
metadata$Diagnosis[match(GSE5392_meta[,2], metadata$Sample)]<-GSE5392_meta$`Disease_status:ch1`
metadata$Age[match(GSE5392_meta[,2], metadata$Sample)]<-as.numeric(GSE5392_meta$`Age (years):ch1`)*365
metadata$Gender[match(GSE5392_meta[,2], metadata$Sample)]<-(GSE5392_meta$`Gender:ch1`)
metadata$Suicide[match(GSE5392_meta[,2], metadata$Sample)]<-(GSE5392_meta$`Suicide:ch1`)
metadata$Alcohol[match(GSE5392_meta[,2], metadata$Sample)]<-GSE5392_meta$`Alcohol abuse (ratings scale:ch1`
metadata$Organism[match(GSE5392_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE5392_meta[,2], metadata$Sample)]<-c(rep("DLPFC", 61), rep("OFC", 21))

## GSE11512
metadata$Region[match(GSE11512_meta[,2], metadata$Sample)]<-GSE11512_meta$`tissue:ch1`
metadata$Age[match(GSE11512_meta[,2], metadata$Sample)]<-as.numeric(gsub(" days", "", GSE11512_meta$`age:ch1`))
metadata$Gender[match(GSE11512_meta[,2], metadata$Sample)]<-GSE11512_meta$`sex:ch1`
metadata$Organism[match(GSE11512_meta[,2], metadata$Sample)]<-GSE11512_meta$organism_ch1

## GSE59288_H
metadata$Region[match(colnames(GSE59288_H), metadata$Sample)]<-GSE59288_H_meta$`tissue:ch1`
metadata$Age[match(colnames(GSE59288_H), metadata$Sample)]<-as.numeric(gsub(" days", "", GSE59288_H_meta$`age:ch1`))
metadata$Diagnosis[match(colnames(GSE59288_H), metadata$Sample)]<-GSE59288_H_meta$`sample group:ch1`
metadata$Gender[match(colnames(GSE59288_H), metadata$Sample)]<-GSE59288_H_meta$`gender:ch1`
metadata$Organism[match(colnames(GSE59288_H), metadata$Sample)]<-"Homo sapiens"

## GSE64810 (all samples male)
metadata$Region[match(GSE64810_meta[,1], metadata$Sample)]<-gsub("tissue: ", "", GSE64810_meta$characteristics_ch1)
metadata$Age[match(GSE64810_meta[,1], metadata$Sample)]<-as.numeric(gsub("age of death: ", "", GSE64810_meta$characteristics_ch1.3))*365
metadata$Diagnosis[match(GSE64810_meta[,1], metadata$Sample)]<-gsub("diagnosis: ", "", GSE64810_meta$characteristics_ch1.1)
metadata$Organism[match(GSE64810_meta[,1], metadata$Sample)]<-"Homo sapiens"
metadata$Gender[match(GSE64810_meta[,1], metadata$Sample)]<-"M"

## GSE5388
metadata$Diagnosis[match(GSE5388_meta[,2], metadata$Sample)]<-GSE5388_meta$`Disease_status:ch1`
metadata$Gender[match(GSE5388_meta[,2], metadata$Sample)]<-gsub("Gender: ","", GSE5388_meta$characteristics_ch1.2)
metadata$Age[match(GSE5388_meta[,2], metadata$Sample)]<-as.numeric(GSE5388_meta$`Age (years):ch1`)*365
metadata$Alcohol[match(GSE5388_meta[,2], metadata$Sample)]<-GSE5388_meta$`Alcohol abuse (ratings scale:ch1`
metadata$Organism[match(GSE5388_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE5388_meta[,2], metadata$Sample)]<-"DLPFC"

## GSE161986
metadata$Region[match(GSE161986_meta[,2], metadata$Sample)]<-GSE161986_meta$`tissue:ch1`
metadata$Age[match(GSE161986_meta[,2], metadata$Sample)]<-as.numeric(GSE161986_meta$`age:ch1`)*365
metadata$Diagnosis[match(GSE161986_meta[,2], metadata$Sample)]<-GSE161986_meta$`diagnosis:ch1`
metadata$Smoke[match(GSE161986_meta[,2], metadata$Sample)]<-GSE161986_meta$`smoking:ch1`
metadata$Gender[match(GSE161986_meta[,2], metadata$Sample)]<-GSE161986_meta$`Sex:ch1`
metadata$Organism[match(GSE161986_meta[,2], metadata$Sample)]<-"Homo sapiens"

## GSE21138
metadata$Diagnosis[match(GSE21138_meta[,2], metadata$Sample)]<-c(rep("CTRL", 29), rep("Schizophrenia", 30))
metadata$Age[match(GSE21138_meta[,2], metadata$Sample)]<-as.numeric(GSE21138_meta$`age:ch1`)*365
metadata$Organism[match(GSE21138_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE21138_meta[,2], metadata$Sample)]<-"PFC"
metadata$Gender[match(GSE21138_meta[,2], metadata$Sample)]<-GSE21138_meta$`Sex:ch1`

## GSE17612
metadata$Age[match(GSE17612_meta[,2], metadata$Sample)]<-as.numeric(gsub(" days", "", GSE17612_meta$`age:ch1`))
metadata$Organism[match(GSE17612_meta[,2], metadata$Sample)]<-GSE17612_meta$organism_ch1
metadata$Region[match(GSE17612_meta[,2], metadata$Sample)]<-GSE17612_meta$`tissue:ch1`
metadata$Gender[match(GSE17612_meta[,2], metadata$Sample)]<-GSE17612_meta$`sex:ch1`

## GSE49376
metadata$Diagnosis[match(GSE49376_meta[,2], metadata$Sample)]<-GSE49376_meta$`aud status:ch1`
metadata$Age[match(GSE49376_meta[,2], metadata$Sample)]<-as.numeric(GSE49376_meta$`age:ch1`)*365
metadata$Organism[match(GSE49376_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE49376_meta[,2], metadata$Sample)]<-"PFC"
metadata$Gender[match(GSE49376_meta[,2], metadata$Sample)]<-GSE49376_meta$`Sex:ch1`

## GSE21935
metadata$Diagnosis[match(GSE21935_meta[,2], metadata$Sample)]<-GSE21935_meta$`disease state:ch1`
metadata$Age[match(GSE21935_meta[,2], metadata$Sample)]<-as.numeric(GSE21935_meta$`age:ch1`)*365
metadata$Organism[match(GSE21935_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE21935_meta[,2], metadata$Sample)]<-GSE21935_meta$`tissue:ch1`
metadata$Gender[match(GSE21935_meta[,2], metadata$Sample)]<-GSE21935_meta$`gender:ch1`

## GSE125681
metadata$Diagnosis[match(GSE125681_meta[,2], metadata$Sample)]<-"Suicide"
metadata$Age[match(GSE125681_meta[,2], metadata$Sample)]<-as.numeric(GSE125681_meta$`age:ch1`)*365
metadata$Organism[match(GSE125681_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE125681_meta[,2], metadata$Sample)]<-"DLPFC"
metadata$Gender[match(GSE125681_meta[,2], metadata$Sample)]<-GSE125681_meta$`gender:ch1`
metadata$Suicide[match(GSE125681_meta[,2], metadata$Sample)]<-"Yes"

## GSE5390
metadata$Age[match(GSE5390_meta[,2], metadata$Sample)]<-as.numeric(GSE5390_meta$`Age (years):ch1`)*365
metadata$Organism[match(GSE5390_meta[,2], metadata$Sample)]<-"Homo sapiens"
metadata$Region[match(GSE5390_meta[,2], metadata$Sample)]<-"PFC"
metadata$Gender[match(GSE5390_meta[,2], metadata$Sample)]<-GSE5390_meta$`Gender:ch1`
metadata$Diagnosis[match(GSE5390_meta[,2], metadata$Sample)]<-GSE5390_meta$`Disease_status:ch1`

## GSE102741
metadata$Age[match(GSE102741_meta$description, metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE102741_meta$characteristics_ch1.4))*365
metadata$Organism[match(GSE102741_meta$description, metadata$Sample)]<-GSE102741_meta$organism_ch1
metadata$Region[match(GSE102741_meta$description, metadata$Sample)]<-"DLPFC"
metadata$Gender[match(GSE102741_meta$description, metadata$Sample)]<-GSE102741_meta$`Sex:ch1`
metadata$Diagnosis[match(GSE102741_meta$description, metadata$Sample)]<-GSE102741_meta$`disease status:ch1`

## GSE174409
metadata$Age[match(GSE174409_meta$title, metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE174409_meta$`age:ch1`))*365
metadata$Organism[match(GSE174409_meta$title, metadata$Sample)]<-GSE174409_meta$organism_ch1
metadata$Smoke[match(GSE174409_meta$title, metadata$Sample)]<-GSE174409_meta$`tobaccoatod:ch1`
metadata$Gender[match(GSE174409_meta$title, metadata$Sample)]<-GSE174409_meta$`Sex:ch1`
metadata$Diagnosis[match(GSE174409_meta$title, metadata$Sample)]<-GSE174409_meta$`diagnosis:ch1`
metadata$Region[match(GSE174409_meta$title, metadata$Sample)]<-GSE174409_meta$source_name_ch1

## GSE99349
metadata$Age[match(GSE99349_meta[,2], metadata$Sample)]<-as.numeric(gsub("age: ", "", GSE99349_meta$`age:ch1`))*365
metadata$Organism[match(GSE99349_meta[,2], metadata$Sample)]<-GSE99349_meta$organism_ch1
metadata$Race[match(GSE99349_meta[,2], metadata$Sample)]<-gsub("race: ", "", GSE99349_meta$characteristics_ch1.3)
metadata$Region[match(GSE99349_meta[,2], metadata$Sample)]<-"DLPFC"
metadata$Diagnosis[match(GSE99349_meta[,2], metadata$Sample)]<-GSE99349_meta$characteristics_ch1.1

#####################################################################
## 4. Harmonisation of the annotations
#####################################################################
## Controls -> "Healthy"; samples without a diagnosis are controls
metadata$Diagnosis[metadata$Diagnosis %in% c("case/control: Unaffected Control", "Control", "CTRL", "CTL", "Healthy control", "Neurologically normal", "non-psychiatric controls (CON)", "control", "healthy", "CONT")]<-"Healthy"
metadata$Diagnosis[is.na(metadata$Diagnosis)]<-"Healthy"

## Diagnoses
metadata$Diagnosis[metadata$Diagnosis %in% c("Major Depressive Disorder", "MDD", "DSM-IV major depressive disorder non-suicides (MDD)", "DSM-IV major depressive disorder suicides (MDD-S)")]<-"MDD"
metadata$Diagnosis[metadata$Diagnosis %in% c("Bipolar disorder", "Bipolar Disorder", "Bipolar")]<-"Bipolar"
metadata$Diagnosis[metadata$Diagnosis %in% c("DS", "Down syndrome")]<-"Down"
metadata$Diagnosis[metadata$Diagnosis %in% c("Alcohol", "Alcohol Abuse", "Alcohol dependence" )]<-"Alcohol"
metadata$Diagnosis[metadata$Diagnosis %in% c("autistic patient", "ASD", "Autism spectrum disorder (ASD)")]<-"ASD"
metadata$Diagnosis[metadata$Diagnosis %in% c("Schizophrenia", "schizophrenic")]<-"Schizophrenia"
metadata$Diagnosis[metadata$Diagnosis %in% c("case/control: Cocaine Use Disorder", "case/control: Cocaine Dependence")]<-"CUD"

## Regions
metadata$Region[metadata$Region %in% c("prefrontal cortex (PFC)", "prefrontal cortex of the brain", "prefrontal cortex", "Prefrontal Cortex")]<-"PFC"
metadata$Region[metadata$Region %in% c("Dorsolateral Prefrontal Cortex", "dorsolateral prefrontal cortex of the brain", "dlPFC")]<-"DLPFC"
metadata$Region[metadata$Region %in% c("Nucleus Accumbens")]<-"NAC"

## Simplified region (BA8 and BA9 are DLPFC)
metadata$Region_simpl<-NA
metadata$Region_simpl[metadata$Region %in% c("dorsolateral prefrontal cortex (DLPFC BA 46/9)","Prefrontal Cortex BA8/9", "Prefrontal Cortex BA9", "Prefrontal Cortex BA8", "DLPFC", "BA9 prefrontal cortex",
                                             "dorsal lateral prefrontal cortex (Brodmann Area 9)", "Postmortem dorsolateral prefrontal cortex (DLPFC BA 46/9)")]<-"DLPFC"
metadata$Region_simpl[metadata$Region %in% c("PFC", "superior frontal gyrus of the brain", "Brain prefrontal cortex, broadman areas 47", "brain (medial prefrontal cortex)",
                                             "medial prefrontal cortex (mPFC)", "superior frontal gyrus", "Brain prefrontal cortex, broadman areas 11", "Brain prefrontal cortex, broadman areas 47")]<-"PFC"
metadata$Region_simpl[metadata$Region %in% c("superior temporal gyrus (STG)", "superior temporal cortex (Brodmann Area 22, BA22)", "temporal cortex")]<-"TC"
metadata$Region_simpl[metadata$Region %in% c("cerebellar cortex", "cerebellar cortex of the brain", "CBC")]<-"CBC"
metadata$Region_simpl[metadata$Region=="DFC"]<-"DLPFC"
metadata$Region_simpl[metadata$Region %in% c("MFC")]<-"PFC"
metadata$Region_simpl[metadata$Region %in% c("ITC", "STC")]<-"TC"
metadata$Region_simpl[metadata$Region %in% c("OFC")]<-"OPFC"
metadata$Region_simpl[metadata$Dataset=="E_PFC"]<-"DLPFC"
metadata$Region[metadata$Dataset=="E_PFC"]<-"DLPFC"
metadata$Region_simpl[metadata$Region %in% c("Brain prefrontal cortex, broadman areas 47",
                                             "Brain prefrontal cortex, broadman areas 11", "Brain prefrontal cortex, broadman areas 47")]<-"OPFC"
metadata$Region_simpl[metadata$Region %in% c("superior frontal gyrus of the brain", "superior frontal gyrus")]<-"SFG"
metadata$Region_simpl[metadata$Region %in% c("VFC Left", "VFC Right", "VFC")]<-"VFC"
metadata$Region_simpl[metadata$Region %in% c("STC Left", "STC Right", "STC")]<-"TC"
metadata$Region_simpl[metadata$Region %in% c("M1C Right", "M1C Left")]<-"M1C"

## Alcohol
metadata$Alcohol[metadata$Alcohol %in% c("alcool: no", "alcool: No", "0 (no use) to 6 (heavy use)): 0")]<- "No"
metadata$Alcohol[metadata$Alcohol %in% c("alcool: yes", "0 (no use) to 6 (heavy use)): 6")]<- "Yes"
metadata$Alcohol[metadata$Alcohol %in% c("alcool: NA", "0 (no use) to 6 (heavy use)): Unknown")]<- NA

## Smoke
metadata$Smoke[metadata$Smoke %in% c("No", "smoking history: No", "no")]<- "No"
metadata$Smoke[metadata$Smoke %in% c("smoking history: Yes", "Yes", "yes (heavy)", "yes(moderate)")]<- "Yes"
metadata$Smoke[metadata$Smoke %in% c("smoking history: Unknown", "Unknown", "NA")]<- NA

## Gender
metadata$Gender[metadata$Gender %in% c("f", "F", "female", "Female")]<- "F"
metadata$Gender[metadata$Gender %in% c("m", "M", "male", "Male")]<- "M"

## Suicide, alcohol, race
metadata$Suicide[metadata$Suicide %in% c("Suicide", "suicide (1=yes): 1")]<-"Yes"
metadata$Suicide[metadata$Suicide %in% c("Accident", "Natural","suicide (1=yes): 0")]<-"No"
metadata$Alcohol[metadata$Alcohol %in% c("0 (no use) to 6 (heavy use)): 1", "0 (no use) to 6 (heavy use)): 2", "0 (no use) to 6 (heavy use)): 3")]<-"No"
metadata$Alcohol[metadata$Alcohol %in% c("0 (no use) to 6 (heavy use)): 4", "0 (no use) to 6 (heavy use)): 5")]<-"Yes"
metadata$Alcohol[metadata$Diagnosis=="Alcohol"]<-"Yes"
metadata$Diagnosis[metadata$Alcohol=="Yes" & metadata$Diagnosis=="Healthy"]<-"Alcohol"
metadata$Race[metadata$Race=="AA"]<-"African American"
metadata$Race[metadata$Race=="B"]<-"African American"
metadata$Race[metadata$Race=="W"]<-"Caucasian"
metadata$Race[metadata$Race=="CAUC"]<-"Caucasian"
metadata$Race[metadata$Race=="HISP"]<-"Hispanic"
metadata$Race[metadata$Race=="H"]<-"Hispanic"   # assumed
metadata$Race[metadata$Race=="AS"]<-"Asian"     # assumed
metadata$Race[metadata$Race=="CC"]<-"Caucasian" # assumed
metadata$Diagnosis[metadata$Suicide=="Yes" & metadata$Diagnosis=="Healthy"]<-"Suicide"
metadata$Race[match(GSE102741_meta$description, metadata$Sample)]<-GSE102741_meta$`race:ch1`

metadata$Gender[metadata$Gender %in% c("f", "F", "female", "Female", "2")]<- "F"
metadata$Gender[metadata$Gender %in% c("m", "M", "male", "Male", "1")]<- "M"
metadata$Gender[metadata$Gender=="5"]<-NA
metadata$Diagnosis[metadata$Diagnosis=="Tics"]<-"OCD"
metadata$Race[metadata$Race %in% c("AA")]<- "African American"
metadata$Race[metadata$Race %in% c("CAUC")]<- "Caucasian"
metadata$Diagnosis[metadata$Diagnosis=="Alcohol"]<-"AUD"
metadata$Diagnosis[metadata$Diagnosis=="Bipolar"]<-"BD"
metadata$Diagnosis[metadata$Diagnosis=="Down"]<-"DS"
metadata$Diagnosis[metadata$Diagnosis=="Huntington's Disease"]<-"HD"
metadata$Diagnosis[metadata$Diagnosis=="Schizophrenia"]<-"SCZ"
metadata$Region_simpl[metadata$Region_simpl=="OPFC"]<-"OFC"
metadata$Region_simpl[is.na(metadata$Region_simpl)]<-"Other"

#####################################################################
## 5. GSE37981, added after the harmonisation
#####################################################################
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
metadata<-rbind.data.frame(metadata, add_meta)

#####################################################################
## 6. Sample names of the expression matrices
#####################################################################
colnames(GSE161986)<-unlist(strsplit( colnames(GSE161986), "_"))[seq(1,ncol(GSE161986)*2,2)]
colnames(GSE11512)<-gsub("[.]CEL[.]gz", "", colnames(GSE11512))
colnames(GSE13564)<-gsub("[.]CEL[.]gz", "", colnames(GSE13564))
colnames(GSE5388)<-gsub("[.]cel[.]gz", "", colnames(GSE5388))
colnames(GSE59630)<-unlist(strsplit( colnames(GSE59630), "_"))[seq(1,ncol(GSE59630)*2,2)]
colnames(GSE60190)<-metadata$Sample[which(metadata$Dataset=="GSE60190")]
colnames(GSE21138)<-gsub("[.]CEL[.]gz", "", colnames(GSE21138))
colnames(GSE49376)<-metadata$Sample[which(metadata$Dataset=="GSE49376")]
colnames(GSE5390)<-gsub("[.]cel[.]gz", "", colnames(GSE5390))
colnames(GSE92538_GPL17027)<-GSE92538_GPL17027_meta[1,-1]

#####################################################################
## 7. Filters and recoding
#####################################################################
## GTEx: keep the samples present in GTEx_PFC
metadata<-metadata[!metadata$Sample %in% setdiff(metadata$Sample[metadata$Dataset=="E_PFC"], colnames(GTEx_PFC)),]
metadata$Dataset[metadata$Dataset=="E_PFC"]<-"GTEx_PFC"

## Homo sapiens only
metadata<-metadata[metadata$Organism=="Homo sapiens",]
## GSE5392 is the superseries of GSE5388 and GSE5390
metadata<-metadata[metadata$Dataset!="GSE5392",]
metadata$Region_simpl[which(metadata$Region=="STG")]<-"TC"
metadata$Diagnosis[which(metadata$Region=="STG")]<-"Healthy"
metadata$Gender[which(metadata$Gender=="female")]<-"F"
metadata$Gender[which(metadata$Gender=="male")]<-"M"

## Datasets with at least 10 samples aged >= 20 years (7300 days)
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

#####################################################################
## 8. GSE33000
#####################################################################
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
dir.create(DOWNLOAD_DIR, showWarnings = FALSE, recursive = TRUE)
f33000 <- file.path(DOWNLOAD_DIR, "GSE33000_series_matrix.txt.gz")
gds <- if(file.exists(f33000)) getGEO(filename = f33000, destdir = DOWNLOAD_DIR, AnnotGPL = TRUE) else
  getGEO("GSE33000", destdir = DOWNLOAD_DIR, AnnotGPL = TRUE)$GSE33000_series_matrix.txt.gz
GSE33000<-exprs(gds)
GSE33000_meta<-pData(gds)
GSE33000_anno<-fData(gds)
GSE33000<-changenames(data=GSE33000, anno=cbind(as.character(GSE33000_anno$ID), GSE33000_anno$`Gene symbol`))

metadata_GSE33000<-data.frame(Sample=colnames(GSE33000),
                              Dataset="GSE33000",
                              Age=as.numeric(gsub(" yrs", "", GSE33000_meta$`age:ch2`))*356,
                              Gender=GSE33000_meta$`gender:ch2`,
                              Race=NA, Suicide=NA,
                              Diagnosis=GSE33000_meta$`disease status:ch2`,
                              Alcohol=NA, Drugs=NA, Smoke=NA, Medication=NA,
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

#####################################################################
## Save
#####################################################################
# metadata keeps all datasets (it is the source of Table S2/S3 and Fig. S1), but
# only the expression matrices that the analysis scripts use are saved.
# GSE125681, GSE17612, GSE21935, GSE22570, GSE37981, GSE71620 and GSE59288_H are
# processed above only to build metadata; no analysis script reads their matrices.
datasets <- c("GSE101521", "GSE102556_H", "GSE102741", "GSE11512", "GSE13564",
              "GSE161986", "GSE174409", "GSE21138", "GSE30272", "GSE33000",
              "GSE49376", "GSE5388", "GSE5390", "GSE59630", "GSE60190",
              "GSE64810", "GSE92538_GPL17027", "GSE99349", "GTEx_PFC")
stopifnot(all(datasets %in% metadata$Dataset))
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
save(list = c(datasets, "metadata"), file = OUT_FILE)
