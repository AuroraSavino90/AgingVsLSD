# =============================================================================
# prep_03_ExerciseAndAlcohol.R  -- data preparation, step 3
#
# Regenerates the differential-expression (DE) tables of the non-psychedelic
# interventions (exercise / enriched environment / cognitive performance and
# alcohol) that downstream scripts load from "Data/original_RData_notclosed/ExerciseAndAlcohol.RData"
# (09_age_reversal_otherpsychedelics.R, 10_GSEA_integrated_otherpsychedelics.R,
#  12_dementia_cohen.R, 19_GO_revert_dementia.R, save_mean_logFC_selected_genes.R).
#
# Provenance: the original Data/ExerciseAndAlcohol.RData is an RStudio workspace
# dump (no save() call exists). Its DE objects were produced interactively by the
# old scripts in Discarded_scripts/Download data/:
#   Exercise.R  - download/preprocessing (getGEO, changenames, DExpr), DE_GSE64607
#   Alcohol.R / Controls.R - DE_GSE72507_0, DE_GSE60676_0, DE_GSE28515_C57
#   Controls.R  - DE_GSE111273_Y/O and DE_GSE164798_chronic (DESeq2)
# (The dump also held DE_GSE105453_Y/O, cognitive performance; that dataset is
#  not part of the paper and is no longer computed.)
# The dump also holds aging datasets and plotting leftovers that downstream
# scripts do not take from this file; they are not regenerated here.
#
# GEO accessions (Mus musculus):
#   GSE64607  voluntary running, lateral entorhinal cortex, day 28 (array, limma)
#   GSE72507  chronic intermittent ethanol, PFC, 0 h (array, limma)
#   GSE60676  ethanol (array, limma)
#   GSE28515  ethanol vs saline, C57 (array, limma)
#   GSE111273 enriched environment vs none, 5 / 24 months (counts, DESeq2)
#   GSE164798 30-day wheel running vs untreated, C57BL/6 (counts, DESeq2)
#
# Output (OUT_FILE): DE_GSE111273_O, DE_GSE111273_Y, DE_GSE164798_chronic, DE_GSE28515_C57, DE_GSE60676_0,
#   DE_GSE64607, DE_GSE72507_0
#
# Run from the project root:  Rscript Scripts/data_preparation/prep_03_ExerciseAndAlcohol.R
# Paths can be overridden with the environment variables RAW_DIR and OUT_FILE.
# =============================================================================

RAW_DIR  <- Sys.getenv("RAW_DIR",  "Data/raw_data_notclosed/exercise_alcohol_GEO")  # GEO files
OUT_FILE <- Sys.getenv("OUT_FILE", "Data/regenerated/ExerciseAndAlcohol.RData")
# Ensembl release used for the Ensembl -> MGI symbol mapping of GSE111273
# (the original used the then-current Ensembl via useMart; see note below).
ENSEMBL_RELEASE <- as.integer(Sys.getenv("ENSEMBL_RELEASE", "107"))
ENSEMBL_CACHE   <- Sys.getenv("ENSEMBL_CACHE", "Data/ensembl")

stopifnot(!normalizePath(OUT_FILE, mustWork = FALSE) %in%   # never overwrite the deposited files
          normalizePath(file.path(c("Data/original_RData", "Data/original_RData_notclosed"), "ExerciseAndAlcohol.RData"), mustWork = FALSE))
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(DESeq2)
})

# -----------------------------------------------------------------------------
# Functions from the original scripts (verbatim except the eBayes legacy flag)
# -----------------------------------------------------------------------------

# Probe/gene ID -> symbol; for duplicated symbols keep the row with highest mean
changenames <- function(data, anno){
  annotation_sel = anno[match(rownames(data), anno[,1]), 2]
  if(length(which(annotation_sel == "")) > 0){
    data <- data[-which(annotation_sel == ""), ]
    annotation_sel <- annotation_sel[-which(annotation_sel == "")]
  }
  a <- which(duplicated(annotation_sel))
  while(length(a) > 0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel == unique(annotation_sel)[i])) > 1){
        m = which.max(rowMeans(data[which(annotation_sel == unique(annotation_sel)[i]), ], na.rm = T))
        data = data[-which(annotation_sel == unique(annotation_sel)[i])[-m], ]
        annotation_sel = annotation_sel[-which(annotation_sel == unique(annotation_sel)[i])[-m]]
      }
    }
    data = data[which(is.na(annotation_sel) == F), ]
    annotation_sel = na.omit(annotation_sel)
    a <- which(duplicated(annotation_sel))
  }
  rownames(data) = annotation_sel
  return(data)
}

# limma two-group comparison; positive logFC = higher in cond1 (treated)
DExpr <- function(cond1, cond2){
  cond1 <- data.frame(cond1)
  cond2 <- data.frame(cond2)
  gs <- factor(c(rep("Cond1", ncol(cond1)), rep("Cond2", ncol(cond2))))
  design <- model.matrix(~gs + 0, cbind.data.frame(cond1, cond2))
  fit <- lmFit(cbind.data.frame(cond1, cond2), design)
  cts <- "gsCond1-gsCond2"
  cont.matrix <- makeContrasts(contrasts = cts, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  # limma >= 3.61 changed the eBayes default when residual df differ between
  # genes (missing values, e.g. GSE60676); legacy = TRUE keeps the original
  # (pre-3.61) behaviour. No effect on datasets without missing values.
  fit2 <- if("legacy" %in% names(formals(eBayes))) eBayes(fit2, 0.01, legacy = TRUE) else eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = nrow(fit2))
  return(tT)
}

# Series matrix (with GPL annotation) from GEO, cached in RAW_DIR. A local copy
# is read directly (getGEO(GSE) would query the GEO file listing every time).
get_series <- function(gse){
  f <- file.path(RAW_DIR, paste0(gse, "_series_matrix.txt.gz"))
  if(file.exists(f)) return(getGEO(filename = f, destdir = RAW_DIR, AnnotGPL = TRUE))
  gds <- getGEO(gse, destdir = RAW_DIR, AnnotGPL = TRUE)
  gds[[paste0(gse, "_series_matrix.txt.gz")]]
}

# Supplementary file from GEO, cached in RAW_DIR
get_suppl <- function(gse, fname){
  f <- file.path(RAW_DIR, fname)
  if(!file.exists(f)) getGEOSuppFiles(gse, makeDirectory = FALSE, baseDir = RAW_DIR,
                                      filter_regex = fname)
  f
}

# -----------------------------------------------------------------------------
# Microarray datasets (limma, Exercise.R / Alcohol.R)
# -----------------------------------------------------------------------------

## GSE64607 - voluntary running vs saline vehicle, lateral entorhinal cortex, day 28
es <- get_series("GSE64607")
GSE64607_meta <- pData(es); GSE64607_anno <- fData(es)
GSE64607 <- changenames(data = exprs(es), anno = cbind(GSE64607_anno$ID, GSE64607_anno$`Gene symbol`))
ctrl <- which(GSE64607_meta$`tissue:ch1` == "lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1` == "saline vehicle" & GSE64607_meta$`timepoint:ch1` == "day 28")
trt  <- which(GSE64607_meta$`tissue:ch1` == "lateral entorhinal cortex" & GSE64607_meta$`intervention:ch1` == "voluntary running" & GSE64607_meta$`timepoint:ch1` == "day 28")
DE_GSE64607 <- DExpr(GSE64607[, trt], GSE64607[, ctrl])

## GSE72507 - chronic intermittent ethanol (CIE) vs control, 0 h after last session
es <- get_series("GSE72507")
GSE72507_meta <- pData(es); GSE72507_anno <- fData(es)
GSE72507 <- changenames(data = exprs(es), anno = cbind(GSE72507_anno$ID, GSE72507_anno$`Gene symbol`))
ctrl <- intersect(grep("control", GSE72507_meta[, 1]), which(GSE72507_meta$`sacrifice time:ch1` == "0 hour post final vapor chamber session"))
trt  <- intersect(grep("CIE",     GSE72507_meta[, 1]), which(GSE72507_meta$`sacrifice time:ch1` == "0 hour post final vapor chamber session"))
DE_GSE72507_0 <- DExpr(GSE72507[, trt], GSE72507[, ctrl])

## GSE60676 - ethanol
es <- get_series("GSE60676")
GSE60676_meta <- pData(es); GSE60676_anno <- fData(es)
GSE60676 <- changenames(data = exprs(es), anno = cbind(GSE60676_anno$ID, GSE60676_anno$`Gene symbol`))
DE_GSE60676_0 <- DExpr(GSE60676[, trt], GSE60676[, ctrl])

## GSE28515 - ethanol vs saline, C57 strain
es <- get_series("GSE28515")
GSE28515_meta <- pData(es); GSE28515_anno <- fData(es)
GSE28515 <- changenames(data = exprs(es), anno = cbind(GSE28515_anno$ID, GSE28515_anno$`Gene symbol`))
ctrl <- which(GSE28515_meta$source_name_ch1 == "saline_C57")
trt  <- which(GSE28515_meta$source_name_ch1 == "etoh_C57")
DE_GSE28515_C57 <- DExpr(GSE28515[, trt], GSE28515[, ctrl])

# -----------------------------------------------------------------------------
# GSE111273 - enriched environment vs none, 5 and 24 months (DESeq2)
# Ensembl gene IDs mapped to MGI symbols. The original used
#   getBM(c("ensembl_gene_id", "mgi_symbol"), filters = "ensembl_gene_id",
#         values = rownames(GSE111273), mart = useMart("ensembl", ...))
# on the live Ensembl; the table left in the original workspace (43,926 pairs)
# is identical to Ensembl 107 (Jul 2022). ensembl_release_map() rebuilds that
# release's table from the Ensembl FTP (see ensembl_release_map.R).
# -----------------------------------------------------------------------------
GSE111273_meta <- pData(get_series("GSE111273"))
GSE111273 <- read.csv(get_suppl("GSE111273", "GSE111273_sample_counts.txt.gz"), sep = "\t")
rownames(GSE111273) <- GSE111273[, 1]

source("Scripts/data_preparation/ensembl_release_map.R")
map.mouse <- ensembl_release_map(ENSEMBL_RELEASE, "mmusculus", ENSEMBL_CACHE)
ids <- map.mouse[map.mouse$ensembl_gene_id %in% rownames(GSE111273), ]
GSE111273 <- changenames(data = GSE111273[, -c(1)], anno = ids)

GSE111273_age <- GSE111273_meta$`age:ch1`
GSE111273_trt <- GSE111273_meta$`perturbation:ch1`

dds <- DESeqDataSetFromMatrix(countData = GSE111273[, GSE111273_age %in% c("5 months")],
                              colData = data.frame(treatment = GSE111273_trt[GSE111273_age %in% c("5 months")]),
                              design = ~ treatment)
dds <- DESeq(dds)
DE_GSE111273_Y <- results(dds, contrast = c("treatment", "enriched environment", "none"))

dds <- DESeqDataSetFromMatrix(countData = GSE111273[, GSE111273_age %in% c("24 months")],
                              colData = data.frame(treatment = GSE111273_trt[GSE111273_age %in% c("24 months")]),
                              design = ~ treatment)
dds <- DESeq(dds)
DE_GSE111273_O <- results(dds, contrast = c("treatment", "enriched environment", "none"))

# -----------------------------------------------------------------------------
# GSE164798 - 30 days voluntary wheel running vs untreated, C57BL/6 (DESeq2)
# -----------------------------------------------------------------------------
GSE164798_meta <- pData(get_series("GSE164798"))
GSE164798 <- read.csv(get_suppl("GSE164798", "GSE164798_Raw_gene_counts_matrix.txt.gz"), sep = "\t")
rownames(GSE164798) <- GSE164798[, 1]
GSE164798 <- GSE164798[, -c(1:6)]

GSE164798_trt    <- GSE164798_meta$`treatment:ch1`
GSE164798_strain <- GSE164798_meta$`strain:ch1`

dds <- DESeqDataSetFromMatrix(countData = GSE164798[, GSE164798_strain %in% c("C57BL/6")],
                              colData = data.frame(treatment = GSE164798_trt[GSE164798_strain %in% c("C57BL/6")]),
                              design = ~ treatment)
dds <- DESeq(dds)
DE_GSE164798_chronic <- results(dds, contrast = c("treatment", "30 days voluntary saucer wheel running", "untreated"))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
save(DE_GSE111273_O, DE_GSE111273_Y,
     DE_GSE164798_chronic, DE_GSE28515_C57, DE_GSE60676_0, DE_GSE64607,
     DE_GSE72507_0, file = OUT_FILE)
message("Saved: ", OUT_FILE)
sessionInfo()
