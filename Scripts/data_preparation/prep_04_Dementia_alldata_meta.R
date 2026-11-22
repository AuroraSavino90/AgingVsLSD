################################################################################
# prep_04_Dementia_alldata_meta.R
#
# Data preparation, step 4: human dementia expression datasets + sample metadata
#
# PRODUCES (OUT_FILE, default Data/regenerated/Dementia_alldata_meta.RData):
#   GSE150696  gene x sample matrix, DLPFC (BA9), Control / AD / DLB / PDD
#   GSE33000   gene x sample matrix, PFC, Control / AD / Huntington
#   GSE44770   gene x sample matrix, DLPFC, Control / AD
#   GSE84422   gene x sample matrix, MSBB multi-region (HG-U133A + HG-U133B)
#   metadata   one row per sample of the four datasets above
#              (Sample, Dataset, Age [days], Gender, Race, Diagnosis,
#               Medication, Region, Platform)
#   These are the only objects of Data/Dementia_alldata_meta.RData read by the
#   downstream scripts (12_dementia_cohen.R, 17_AD_DEGs.R,
#   19_GO_revert_dementia.R, save_mean_logFC_selected_genes.R), which select
#   datasets via metadata$Region in {DLPFC, PFC} and metadata$Diagnosis.
#   Other datasets in the original file (GSE110731, GSE53697, GSE64810,
#   GSE193391, GSE129308, GSE159940, GSE71191, GSE110732) have no Region
#   annotation, are never selected downstream, and are not regenerated.
#
# PROVENANCE: condensed from the original scripts
#   Dementia_alldata.R  (download + probe->gene)
#   Dementia_meta.R     (metadata table)
#   (not distributed)
#   Processing is unchanged: the GEO series-matrix values are used as deposited
#   (no re-normalization); probes are mapped to gene symbols with the GEO
#   annotation (GPL .annot for Affymetrix/Agilent, GPL SOFT for GPL17585), and
#   multiple probes per gene are collapsed keeping the probe with the highest
#   mean expression (changenames()).
#
# DATA SOURCES (all public, NCBI GEO; no controlled-access data):
#   GSE150696 (GPL17585, HTA 2.0; RMA by submitters)
#   GSE33000  (GPL4372, Rosetta/Merck 44K; HBTRC)
#   GSE44770  (GPL4372; HBTRC)
#   GSE84422  (GPL96 + GPL97; Mount Sinai Brain Bank, RMA + covariate-corrected
#              by submitters). The GPL570 sub-series is not used (no PFC).
#   Annotation files: GPL4372.annot.gz, GPL96.annot.gz, GPL97.annot.gz,
#   GPL17585.soft.gz. Files are read from RAW_DIR if present, otherwise
#   downloaded there by GEOquery.
#
# RUN ORDER: step 4 of data preparation; run from the project root:
#   Rscript Scripts/data_preparation/prep_04_Dementia_alldata_meta.R [RAW_DIR] [OUT_FILE]
#   Runtime ~85 min (the original changenames() loop is slow on the ~925k
#   GSE150696 probesets); ~3 GB RAM.
# VERIFIED (R 4.5.1, GEOquery 2.76.0, GEO files of Nov 2024): matrices match the
#   original to <= 2e-15 (floating-point parsing), metadata identical.
################################################################################

## ---- configuration ----------------------------------------------------------
RAW_DIR  <- "Data/raw_data_notclosed/dementia_GEO" # GEO files of the original run (missing ones are downloaded here)
OUT_FILE <- "Data/regenerated/Dementia_alldata_meta.RData"

args <- commandArgs(trailingOnly = TRUE)                 # optional overrides
if (length(args) >= 1) RAW_DIR  <- args[1]
if (length(args) >= 2) OUT_FILE <- args[2]

stopifnot(!normalizePath(OUT_FILE, mustWork = FALSE) %in%   # never overwrite the deposited files
          normalizePath(file.path(c("Data/original_RData", "Data/original_RData_notclosed"), "Dementia_alldata_meta.RData"), mustWork = FALSE))
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(bit64)  # GPL4372 probe IDs are read as integer64; needed for as.character()
})
options(timeout = 3600)

## ---- helpers ----------------------------------------------------------------
# Load a GEO series matrix (ExpressionSet with GPL annotation in fData).
# Uses the local file in RAW_DIR when available, otherwise downloads it.
get_series <- function(gse, file) {
  f <- file.path(RAW_DIR, file)
  if (file.exists(f)) {
    getGEO(filename = f, destdir = RAW_DIR, AnnotGPL = TRUE)
  } else {
    getGEO(gse, destdir = RAW_DIR, AnnotGPL = TRUE)[[file]]
  }
}

# Probe -> gene symbol; for duplicated symbols keep the row with the highest
# mean expression; drop unannotated probes. Verbatim from Dementia_alldata.R.
changenames <- function(data, anno){
  annotation_sel = anno[match(rownames(data), anno[,1]), 2]

  if(length(which(annotation_sel == "")) > 0){
    data <- data[-which(annotation_sel == ""),]
    annotation_sel <- annotation_sel[-which(annotation_sel == "")]
  }

  a <- which(duplicated(annotation_sel))
  while(length(a) > 0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel == unique(annotation_sel)[i])) > 1){
        m = which.max(rowMeans(data[which(annotation_sel == unique(annotation_sel)[i]),], na.rm = T))
        data = data[-which(annotation_sel == unique(annotation_sel)[i])[-m],]
        annotation_sel = annotation_sel[-which(annotation_sel == unique(annotation_sel)[i])[-m]]
      }
    }

    data = data[which(is.na(annotation_sel) == F),]
    annotation_sel = na.omit(annotation_sel)
    a <- which(duplicated(annotation_sel))
  }

  rownames(data) = annotation_sel
  return(data)
}

## ---- GSE84422 (MSBB; GPL96 = HG-U133A, GPL97 = HG-U133B) --------------------
es <- get_series("GSE84422", "GSE84422-GPL96_series_matrix.txt.gz")
GSE84422_meta_2 <- pData(es)
anno <- fData(es)
GSE84422_2 <- changenames(exprs(es), cbind(anno$ID, anno$`Gene symbol`))

es <- get_series("GSE84422", "GSE84422-GPL97_series_matrix.txt.gz")
anno <- fData(es)
GSE84422_3 <- changenames(exprs(es), cbind(anno$ID, anno$`Gene symbol`))

# Same samples on both arrays: stack genes, then keep one row per gene
# (highest mean) for genes measured on both arrays.
GSE84422 <- rbind(GSE84422_2, GSE84422_3)
GSE84422 <- changenames(GSE84422, cbind(rownames(GSE84422), rownames(GSE84422)))
rm(GSE84422_2, GSE84422_3)

## ---- GSE33000 (HBTRC, PFC) ---------------------------------------------------
es <- get_series("GSE33000", "GSE33000_series_matrix.txt.gz")
GSE33000_meta <- pData(es)
anno <- fData(es)
GSE33000 <- changenames(exprs(es), cbind(as.character(anno$ID), anno$`Gene symbol`))

## ---- GSE44770 (HBTRC, DLPFC) -------------------------------------------------
es <- get_series("GSE44770", "GSE44770_series_matrix.txt.gz")
GSE44770_meta <- pData(es)
anno <- fData(es)
GSE44770 <- changenames(exprs(es), cbind(as.character(anno$ID), anno$`Gene symbol`))

## ---- GSE150696 (DLPFC BA9; AD, DLB, PDD) -------------------------------------
es <- get_series("GSE150696", "GSE150696_series_matrix.txt.gz")
GSE150696_meta <- pData(es)
anno <- fData(es)
GSE150696 <- changenames(exprs(es), cbind(anno$ID, anno$gene_symbols))
rm(es, anno); invisible(gc())

## ---- metadata (Dementia_meta.R) -----------------------------------------------
# Dataset order as in the original (alphabetical object names).
categories    <- c("Sample", "Dataset", "Age", "Gender", "Race", "Diagnosis",
                   "Medication", "Region", "Platform")
datasetsNames <- c("GSE150696", "GSE33000", "GSE44770", "GSE84422")

metadata <- data.frame(matrix(nrow = sum(sapply(datasetsNames, function(d) ncol(get(d)))),
                              ncol = length(categories)))
colnames(metadata) <- categories
metadata$Dataset <- unlist(lapply(datasetsNames, function(d) rep(d, ncol(get(d)))))
metadata$Sample  <- unlist(lapply(datasetsNames, function(d) colnames(get(d))))

# Age in days. Factors 356 and 365 are kept exactly as in the original script.
# GSE84422 (GPL96 sample annotation)
idx <- match(GSE84422_meta_2[,2], metadata$Sample)
metadata$Age[idx]       <- as.numeric(GSE84422_meta_2$`age:ch1`)*356
metadata$Gender[idx]    <- GSE84422_meta_2$`Sex:ch1`
metadata$Region[idx]    <- GSE84422_meta_2$`brain region:ch1`
metadata$Race[idx]      <- GSE84422_meta_2$`race:ch1`
metadata$Diagnosis[idx] <- GSE84422_meta_2$`neuropathological category:ch1`

# GSE33000
idx <- match(GSE33000_meta[,2], metadata$Sample)
metadata$Age[idx]       <- as.numeric(gsub(" yrs", "", GSE33000_meta$`age:ch2`))*356
metadata$Gender[idx]    <- GSE33000_meta$`gender:ch2`
metadata$Region[idx]    <- GSE33000_meta$`tissue:ch1`
metadata$Diagnosis[idx] <- GSE33000_meta$`disease status:ch2`

# GSE44770
idx <- match(GSE44770_meta[,2], metadata$Sample)
metadata$Age[idx]       <- as.numeric(gsub("age: ", "", GSE44770_meta$characteristics_ch2.1))*356
metadata$Gender[idx]    <- gsub("gender: ", "", GSE44770_meta$characteristics_ch2.2)
metadata$Region[idx]    <- "DLPFC"
metadata$Diagnosis[idx] <- gsub("disease: ", "", GSE44770_meta$characteristics_ch2)

# GSE150696
idx <- match(GSE150696_meta[,2], metadata$Sample)
metadata$Age[idx]       <- as.numeric(GSE150696_meta$`age:ch1`)*365
metadata$Gender[idx]    <- GSE150696_meta$`Sex:ch1`
metadata$Region[idx]    <- "prefrontal cortex BA9"
metadata$Diagnosis[idx] <- GSE150696_meta$description

# Harmonize labels
metadata$Diagnosis[metadata$Diagnosis %in% c("Control", "CTRL", "N", "non-demented", "Normal")] <- "Healthy"
metadata$Diagnosis[metadata$Diagnosis %in% c("A", "AD", "Alzheimer's disease", "definite AD")] <- "Alzheimer"
metadata$Diagnosis[metadata$Diagnosis %in% c("Huntington's disease")] <- "Huntington"
metadata$Gender[metadata$Gender %in% c("F", "f", "Female")] <- "F"
metadata$Gender[metadata$Gender %in% c("M", "m", "Male")]   <- "M"
metadata$Region[metadata$Region %in% c("prefrontal cortex brain", "Prefrontal Cortex")] <- "PFC"
metadata$Region[metadata$Region %in% c("prefrontal cortex BA9", "Dorsolateral Prefrontal Cortex", "DLPFC")] <- "DLPFC"

## ---- save ---------------------------------------------------------------------
save(GSE150696, GSE33000, GSE44770, GSE84422, metadata, file = OUT_FILE)
message("Saved: ", OUT_FILE)
print(table(metadata$Dataset, metadata$Region %in% c("DLPFC", "PFC")))
