# =============================================================================
# prep_02_Psychedelics_PFC.R  -- data preparation, step 2
#
# Regenerates the differential-expression (DE) tables of psychedelics and other
# fast-acting antidepressants in prefrontal cortex that downstream scripts load
# from "Data/original_RData_notclosed/Psychedelics_PFC.RData" 
# plus the LSD count matrix GSE179379 used by 00_LSD_DEGs.R.
#
#
# Input data (RAW_DIR, default Data/raw_data_notclosed/psychedelics_GEO)
# Non-public input data (RAW_CLOSED_DIR, default Data/raw_data_closed/psychedelics):
#   MDMA data/PREFRONTAL CORTEX/PF1-9.CEL  MDMA, mouse PFC, Mouse 430 2.0 arrays (non-GEO, RMA)
#   readcount copy 2.xlsx  DMT / pharmahuasca / harmaline, rat PFC, RNA-seq (unpublished)
# Series matrices are read from RAW_DIR when present (copied to a temp cache so
# nothing is written into RAW_DIR), otherwise downloaded from GEO.
#
# Output (OUT_FILE), only the objects the analysis scripts use: DE_MDMA,
#   DE_GSE23728, DE_GSE26364, DE_GSE179380, DE_GSE161626, DE_GSE161626_48h,
#   DE_GSE161626_7d, DE_GSE81672, DE_DMT, DE_pharm, DE_harm,
#   DE_GSE209859_FC3hours, DE_GSE209859_FC4weeks, GSE179379 (LSD counts).
#   Without the non-public data (RAW_CLOSED_DIR missing or incomplete) DE_MDMA and
#   DE_DMT, DE_pharm, DE_harm are left out (the analysis always reads the deposited
#   input RData, which contain them: run_all.R).
#
# Run from the project root:  Rscript Scripts/data_preparation/prep_02_Psychedelics_PFC.R
# Paths can be overridden with the environment variables RAW_DIR, RAW_CLOSED_DIR,
# OUT_FILE, ENSEMBL_RELEASE and ENSEMBL_CACHE.
# Required packages: GEOquery, affy, limma, DESeq2, openxlsx,
#   mouse4302.db, mouse4302cdf, mogene10stv1cdf (Bioconductor).
# Annotation versions: mouse4302.db takes gene symbols from org.Mm.eg.db, which
# changes every Bioconductor release; org.Mm.eg.db 3.14.0 (Bioconductor 3.14)
# reproduces the original MDMA gene set exactly (later versions rename/drop a few
# hundred probes' symbols).
# =============================================================================

RAW_DIR  <- Sys.getenv("RAW_DIR", "Data/raw_data_notclosed/psychedelics_GEO")
RAW_CLOSED_DIR <- Sys.getenv("RAW_CLOSED_DIR", "Data/raw_data_closed/psychedelics")
OUT_FILE <- Sys.getenv("OUT_FILE", "Data/regenerated/Psychedelics_PFC.RData")
# Ensembl release for the Ensembl -> gene symbol mappings (GSE161626,
# GSE81672, GSE209859, DMT); see the note at the Ensembl section.
ENSEMBL_RELEASE <- as.integer(Sys.getenv("ENSEMBL_RELEASE", "107"))
ENSEMBL_CACHE   <- Sys.getenv("ENSEMBL_CACHE", "Data/ensembl")

stopifnot(!normalizePath(OUT_FILE, mustWork = FALSE) %in%   # never overwrite the deposited files
          normalizePath(file.path(c("Data/original_RData", "Data/original_RData_notclosed"), "Psychedelics_PFC.RData"), mustWork = FALSE))
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
GEO_CACHE <- file.path(tempdir(), "geo_cache")
dir.create(GEO_CACHE, showWarnings = FALSE)
options(timeout = 3600)

suppressPackageStartupMessages({
  library(GEOquery)
  library(affy)
  library(limma)
  library(DESeq2)
  library(openxlsx)
  library(mouse4302.db)   # MDMA probe annotation; CDF packages mouse4302cdf and
                          # mogene10stv1cdf are loaded automatically by affy
})

raw <- function(...) file.path(RAW_DIR, ...)
raw_closed <- function(...) file.path(RAW_CLOSED_DIR, ...)

# Non-public data (not deposited; available from the authors of the original
# studies): MDMA arrays (Fernandez-Castillo et al. 2012) and DMT / pharmahuasca /
# harmaline counts (Kelley et al. 2022). When they are not in RAW_CLOSED_DIR, their DE
# tables (DE_MDMA; DE_DMT, DE_pharm, DE_harm) are not built and not saved.
HAS_MDMA <- dir.exists(raw_closed("MDMA data", "PREFRONTAL CORTEX"))
HAS_DMT  <- file.exists(raw_closed("readcount copy 2.xlsx"))
if (!HAS_MDMA) message("MDMA data not found in RAW_CLOSED_DIR: DE_MDMA is not built")
if (!HAS_DMT)  message("readcount copy 2.xlsx not found in RAW_CLOSED_DIR: DE_DMT, DE_pharm, DE_harm are not built")

# Series matrix via getGEO, using the local copy in RAW_DIR if available
get_series <- function(gse, getGPL = FALSE) {
  f <- paste0(gse, "_series_matrix.txt.gz")
  if (file.exists(raw(f))) {   # local copy: read it directly (no GEO query)
    file.copy(raw(f), file.path(GEO_CACHE, f), overwrite = TRUE)
    return(getGEO(filename = file.path(GEO_CACHE, f), destdir = GEO_CACHE, getGPL = getGPL))
  }
  gds <- getGEO(gse, destdir = GEO_CACHE, getGPL = getGPL)
  gds[[1]]
}

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

# Probe/gene ID -> symbol; for duplicated symbols keep the row with highest mean
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

# limma DE, cond1 vs cond2 (columns = samples)
DExpr<-function(cond1, cond2){
  cond1<-data.frame(cond1)
  cond2<-data.frame(cond2)
  gs <- factor(c(rep("Cond1",ncol(cond1)),rep("Cond2",ncol(cond2))))
  design <- model.matrix(~gs + 0, cbind.data.frame(cond1, cond2))
  fit <- lmFit(cbind.data.frame(cond1, cond2), design)
  cts<-"gsCond1-gsCond2"
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
  return(tT)
}

# -----------------------------------------------------------------------------
# Ensembl ID -> symbol tables. 
# -----------------------------------------------------------------------------
source("Scripts/data_preparation/ensembl_release_map.R")
map.mouse <- ensembl_release_map(ENSEMBL_RELEASE, "mmusculus",   ENSEMBL_CACHE)
map.rat   <- ensembl_release_map(ENSEMBL_RELEASE, "rnorvegicus", ENSEMBL_CACHE)

# =============================================================================
# Preprocessing
# =============================================================================

## MDMA (mouse 430 2.0 CEL files, RMA) ----------------------------------------
if (HAS_MDMA) {
x <- mouse4302SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
xx<-cbind(names(unlist(xx)), unlist(xx))

Data <- ReadAffy(celfile.path = raw_closed("MDMA data", "PREFRONTAL CORTEX"))
MDMA <- rma(Data)
MDMA<-changenames(exprs(MDMA), xx)
}

## GSE161626 (DOI, RNA-seq; one featureCounts file per sample) ----------------
files<-list.files(path = raw("GSE161626"), pattern=".txt")

countdata <- read.delim(raw("GSE161626", files[1]), skip=1, row.names=1)
data<-countdata[,6]
for(i in 2:length(files)){
  countdata <- read.delim(raw("GSE161626", files[i]), skip=1, row.names=1)
  data<-cbind(data, countdata[,6])
}
rownames(data)<-rownames(countdata)
colnames(data)<-files
GSE161626<-data

# all mouse genes
geni_coord<-map.mouse
annotation_sel=geni_coord[match((rownames(GSE161626)), geni_coord[,1]),2]
GSE161626<-changenames(GSE161626, anno=cbind(rownames(GSE161626), annotation_sel))

GSE161626_conditions<-rep(c("Ctrl", "24h", "48h", "7d"), each=8)

## GSE14720 / GSE23728 (DOI, TIFR-Rat-8X15K) ----------------------------------
# Annotation for platform TIFR-Rat-8X15K v1.0 downloaded from
# https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=913 on 16/09/2021
GPL8160_anno<-read.delim(raw("GPL8160_noParents.an.txt"), sep="\t", skip=8)

gds<-get_series("GSE14720")
GSE14720<-exprs(gds)
GSE14720_meta<-pData(gds)
GSE14720<-changenames(data=GSE14720, anno=cbind(GPL8160_anno$ProbeName, GPL8160_anno$GeneSymbols))
GSE14720<-log2(GSE14720)
GSE14720_treatment<-GSE14720_meta$`treatment:ch1`
GSE14720_group<-GSE14720_meta$`group:ch1`

gds<-get_series("GSE23728")
GSE23728<-exprs(gds)
GSE23728_meta<-pData(gds)
GSE23728<-changenames(data=GSE23728, anno=cbind(GPL8160_anno$ProbeName, GPL8160_anno$GeneSymbols))
GSE23728<-log2(GSE23728)
GSE23728_treatment<-GSE23728_meta$`treatment:ch1`
GSE23728_group<-GSE23728_meta$`group:ch1`

## GSE81672 (ketamine, RNA-seq) -----------------------------------------------
gds<-get_series("GSE81672")
GSE81672<-read.xlsx(raw("GSE81672_HTSeq_counts_matrix.xlsx"),1, rowNames = T)
GSE81672_meta<-pData(gds)

ids<-map.mouse[map.mouse$ensembl_gene_id %in% rownames(GSE81672), ]
GSE81672<-changenames(data=GSE81672, anno=ids)

## GSE26364 (ketamine, MoGene 1.0 ST CEL files, RMA) ---------------------------
gpl <- raw("GPL6246.soft.gz")   # platform annotation for fData (gene_assignment)
if (file.exists(gpl)) file.copy(gpl, file.path(GEO_CACHE, basename(gpl)), overwrite = TRUE)
gds<-get_series("GSE26364", getGPL = TRUE)
GSE26364_anno<-fData(gds)

Data <- ReadAffy(celfile.path = raw("GSE26364", "GSE26364_RAW"))
GSE26364 <- rma(Data)

symbols<-gsub('([^/]*)//([^/]*)//.*', '\\2', GSE26364_anno$gene_assignment)
symbols<-gsub(' ', '', symbols)

GSE26364<-exprs(GSE26364)
GSE26364<-changenames(data=GSE26364, anno=cbind(GSE26364_anno$ID, symbols))
GSE26364_treatment<-c("Saline", "Saline", "Ketamine", "Ketamine")

## GSE179379 (LSD, mouse, RNA-seq counts) --------------------------------------
GSE179379<-read.csv(raw("GSE179379", "Raw_gene_counts_matrix.txt"), sep=" ", row.names=NULL)
nomi<-GSE179379[,1]
GSE179379<-GSE179379[,-c(1)]
GSE179379<-changenames(data=GSE179379, anno=cbind(rownames(GSE179379), nomi))

## GSE179380 (LSD, rat RG_U34A, RMA matrix from the authors) -------------------
GSE179380<-read.xlsx(raw("GSE179380", "Matrix_rma.xlsx"), 1, rowNames=T)
annotation<-read.csv(raw("GPL85-10790_data.txt"), sep="\t")
GSE179380<-changenames(data=as.matrix(GSE179380), anno=annotation[,c(1,11)])

## GSE209859 (psilocybin, RNA-seq counts) --------------------------------------
GSE209859<-read.csv(raw("GSE209859_counts.csv.gz"))
names_psilo<-GSE209859[,1]
GSE209859<-GSE209859[,-c(1:6)]
rownames(GSE209859)<-names_psilo
rownames(GSE209859)<-gsub("\\..*","",rownames(GSE209859))

ids<-map.mouse[map.mouse$ensembl_gene_id %in% rownames(GSE209859), ]
GSE209859<-changenames(data=GSE209859, anno=ids)

## DMT / pharmahuasca / harmaline (rat, RNA-seq read counts) --------------------
if (HAS_DMT) {
DMT<-read.xlsx(raw_closed("readcount copy 2.xlsx"),1, rowNames = T)
# all rat genes (getBM without filter, as in the original)
geni_coord<-map.rat
annotation_sel=geni_coord[match((rownames(DMT)), geni_coord[,1]),2]
DMT<-changenames(DMT, anno=cbind(rownames(DMT), annotation_sel))
}

# =============================================================================
# Differential expression (DEGs_PFC.R)
# =============================================================================

## Arrays: limma (DExpr) -------------------------------------------------------
if (HAS_MDMA) {
cond2<-MDMA[,7:9]   # ctrl
cond1<-MDMA[,4:6]   # passive administration
DE_MDMA<-DExpr(cond1, cond2)
}

# GSE14720 and GSE23728: control (non maternally stressed) animals only
cond2<-GSE14720[,GSE14720_group=="control" & GSE14720_treatment == "saline"]
cond1<-GSE14720[,GSE14720_group=="control" & GSE14720_treatment == "DOI"]
DE_GSE14720<-DExpr(cond1, cond2)

cond2<-GSE23728[,GSE23728_group=="control" & GSE23728_treatment == "saline"]
cond1<-GSE23728[,GSE23728_group=="control" & GSE23728_treatment == "DOI"]
DE_GSE23728<-DExpr(cond1, cond2)

cond2<-GSE26364[,GSE26364_treatment=="Saline"]
cond1<-GSE26364[,GSE26364_treatment=="Ketamine"]
DE_GSE26364<-DExpr(cond1, cond2)

cond2<-GSE179380[,1:3]
cond1<-GSE179380[,4:6]
DE_GSE179380<-DExpr(cond1, cond2)

## RNA-seq: DESeq2 -------------------------------------------------------------
# GSE161626: each time point vs Ctrl, two technical libraries per replicate
dds <- DESeqDataSetFromMatrix(countData = GSE161626[,GSE161626_conditions %in% c("24h", "Ctrl")],
                              colData = data.frame(treatment=GSE161626_conditions[GSE161626_conditions %in% c("Ctrl", "24h")], replicate=factor(c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4))),
                              design= ~ treatment+replicate)
dds <- DESeq(dds)
DE_GSE161626 <- results(dds, contrast = c("treatment","24h", "Ctrl"))

dds <- DESeqDataSetFromMatrix(countData = GSE161626[,GSE161626_conditions %in% c("48h", "Ctrl")],
                              colData = data.frame(treatment=GSE161626_conditions[GSE161626_conditions %in% c("Ctrl", "48h")], replicate=factor(c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4))),
                              design= ~ treatment+replicate)
dds <- DESeq(dds)
DE_GSE161626_48h <- results(dds, contrast = c("treatment","48h", "Ctrl"))

dds <- DESeqDataSetFromMatrix(countData = GSE161626[,GSE161626_conditions %in% c("7d", "Ctrl")],
                              colData = data.frame(treatment=GSE161626_conditions[GSE161626_conditions %in% c("Ctrl", "7d")], replicate=factor(c(1,1,2,2,3,3,4,4,1,1,2,2,3,3,4,4))),
                              design= ~ treatment+replicate)
dds <- DESeq(dds)
DE_GSE161626_7d <- results(dds, contrast = c("treatment","7d", "Ctrl"))

# GSE179379 (same call as Scripts/00_LSD_DEGs.R)
dds <- DESeqDataSetFromMatrix(countData = GSE179379,
                              colData = data.frame(treatment=rep(c("LSD", "Saline"), each=10)),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_GSE179379 <- results(dds, contrast = c("treatment","LSD", "Saline"))

# GSE81672: PFC, susceptible ketamine responders vs susceptible saline
GSE81672_PFC<-GSE81672[,GSE81672_meta$`brain region:ch1`=="Prefrontal Cortex"]
GSE81672_meta_PFC<-GSE81672_meta[GSE81672_meta$`brain region:ch1`=="Prefrontal Cortex",]
dds <- DESeqDataSetFromMatrix(countData = GSE81672_PFC,
                              colData = data.frame(group=GSE81672_meta_PFC$`group:ch1`),
                              design= ~ group)
dds <- DESeq(dds)
DE_GSE81672 <- results(dds, contrast = c("group","Susceptible, ketamine-treated, responder", "Susceptible, saline-treated"))

# DMT, pharmahuasca, harmaline vs CTRL (3 samples each)
if (HAS_DMT) {
dds <-DESeqDataSetFromMatrix(countData = DMT[,c(2,7,15,4,12,16,
                                                 3,10,14, 5,9,13)],
                              colData = data.frame(treatment=c("CTRL","CTRL","CTRL",
                                                               "DMT","DMT","DMT",
                                                               "Pharmahuasca", "Pharmahuasca", "Pharmahuasca",
                                                               "Harmaline", "Harmaline", "Harmaline")),
                              design= ~ treatment)
dds <- DESeq(dds)
DE_pharm <- results(dds, contrast = c("treatment","Pharmahuasca", "CTRL"))
DE_DMT <- results(dds, contrast = c("treatment","DMT", "CTRL"))
DE_harm <- results(dds, contrast = c("treatment","Harmaline", "CTRL"))
}

# GSE209859: frontal cortex, psilocybin vs saline at 3 hours and 4 weeks
region<-rep(NA, ncol(GSE209859))
region[grep("Hypothalamus", colnames(GSE209859))]<-"Hypothalamus"
region[grep("Frontal", colnames(GSE209859))]<-"Cortex"

treatment<-rep(NA, ncol(GSE209859))
treatment[grep("Psilocybin", colnames(GSE209859))]<-"Psilocybin"
treatment[grep("Saline", colnames(GSE209859))]<-"Saline"

time<-rep(NA, ncol(GSE209859))
time[grep("weeks", colnames(GSE209859))]<-"4weeks"
time[grep("hours", colnames(GSE209859))]<-"3hours"

anno<-data.frame(region=region,treatment=treatment, time=time)

dds_FC3hours <- DESeqDataSetFromMatrix(countData = GSE209859[,anno$time=="3hours" & anno$region=="Cortex"],
                                       colData = anno[anno$time=="3hours" & anno$region=="Cortex",],
                                       design= ~ treatment)
dds_FC3hours <- DESeq(dds_FC3hours)
DE_GSE209859_FC3hours <- results(dds_FC3hours, contrast = c("treatment","Psilocybin", "Saline"))

dds_FC4weeks <- DESeqDataSetFromMatrix(countData = GSE209859[,anno$time=="4weeks" & anno$region=="Cortex"],
                                       colData = anno[anno$time=="4weeks" & anno$region=="Cortex",],
                                       design= ~ treatment)
dds_FC4weeks <- DESeq(dds_FC4weeks)
DE_GSE209859_FC4weeks <- results(dds_FC4weeks, contrast = c("treatment","Psilocybin", "Saline"))

# =============================================================================
out_objects <- c("DE_MDMA", "DE_GSE23728", "DE_GSE26364", "DE_GSE179380",
     "DE_GSE161626", "DE_GSE161626_48h", "DE_GSE161626_7d", "DE_GSE81672",
     "DE_DMT", "DE_pharm", "DE_harm",
     "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks", "GSE179379")
if (!HAS_MDMA) out_objects <- setdiff(out_objects, "DE_MDMA")
if (!HAS_DMT)  out_objects <- setdiff(out_objects, c("DE_DMT", "DE_pharm", "DE_harm"))
save(list = out_objects, file = OUT_FILE)
message("Saved: ", OUT_FILE)
