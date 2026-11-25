#####################################################################
## Data preparation, step 1a: expression matrices of the human aging
## datasets contained in Data/Alldata_20Sep.RData
#####################################################################
## Provenance
## Data/Alldata_20Sep.RData was built in 2022 

## This script and prep_01b_* generate Data/Alldata_20Sep.RData, originally built in 2022,
##  with the same processing. Checked against the
## original objects: identical (max abs difference < 1e-13) for all datasets.
##
## Input: GEO_DIR, a folder with the files downloaded from GEO
##   - <GSE>_RAW folders with CEL files (or <GSE>/<GSE>_RAW)
##   - series matrices (<GSE>_series_matrix.txt.gz) and GPL annotation files
##     (GPLxxx.annot.gz / GPLxxx.soft) as downloaded in 2022: the gene symbols
##     come from these annotation files, newer versions give different symbols
##   - supplementary processed files (listed below)
## Output: OUT_FILE with one expression matrix per dataset (used by prep_01b).
##
## Required packages: GEOquery, affy, oligo, limma, openxlsx,
## illuminaio (GSE125681 only), and the CDF / pd.* annotation packages that
## affy and oligo install automatically on first use.
#####################################################################

# Run from the repository root:
#   Rscript Scripts/data_preparation/prep_01a_aging_expression.R
# Settings can be overridden with environment variables of the same name.
GEO_DIR  <- Sys.getenv("GEO_DIR", "Data/raw_data_notclosed/aging_GEO")
OUT_FILE <- Sys.getenv("OUT_FILE", "Data/regenerated/aging_expression.RData")
# Ensembl release used for Ensembl ID -> HGNC symbol (section 1). 
ENSEMBL_RELEASE <- as.integer(Sys.getenv("ENSEMBL_RELEASE", "105"))
ENSEMBL_CACHE   <- Sys.getenv("ENSEMBL_CACHE", "Data/ensembl")

suppressMessages({
  library(GEOquery)
  library(affy)
  library(limma)
  library(openxlsx)
})

# Probe -> gene symbol; for duplicated symbols keeps the probe with the
# highest mean expression 
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

g <- function(...) file.path(GEO_DIR, ...)

# Probe annotation of a GEO series (fData with GPL annotation, restricted to
# the probes of the series), read from the local files in GEO_DIR
series_anno <- function(file, annot = TRUE){
  fData(getGEO(filename = g(file), AnnotGPL = annot, destdir = GEO_DIR))
}

#####################################################################
## 1. Datasets with Ensembl gene IDs: HGNC symbols of Ensembl 105
#####################################################################
source("Scripts/data_preparation/ensembl_release_map.R")
ens_map <- ensembl_release_map(ENSEMBL_RELEASE, "hsapiens", ENSEMBL_CACHE)
ens_symbols <- function(ids) ens_map[ens_map$ensembl_gene_id %in% ids, ]

# GSE101521: counts
GSE101521 <- read.csv(g("GSE101521_totalRNA_counts.csv.gz"), row.names = 1)
ids <- ens_symbols(rownames(GSE101521))
GSE101521 <- changenames(data = GSE101521, anno = cbind(rownames(GSE101521), ids[match(rownames(GSE101521), ids[,1]), 2]))
GSE101521 <- log2(GSE101521 + 1)

# GSE59288: human samples (first 34 columns)
GSE59288_H <- read.csv(g("GSE59288_exp_mRNA.txt.gz"), sep = "\t", row.names = 1)[, 1:34]
ids <- ens_symbols(rownames(GSE59288_H))
GSE59288_H <- changenames(data = GSE59288_H, anno = cbind(rownames(GSE59288_H), ids[match(rownames(GSE59288_H), ids[,1]), 2]))
GSE59288_H <- log2(GSE59288_H + 1)
# dataset (34 ASD samples, superior frontal gyrus) is never selected by the
# analyses (healthy samples only).
GSE59288_H <- log2(GSE59288_H + 1)

# GSE64810: DESeq2-normalised counts, versioned Ensembl IDs
GSE64810 <- read.csv(g("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz"), sep = "\t", row.names = 1)
ids <- ens_symbols(gsub('\\.[0-9]', '', rownames(GSE64810)))
GSE64810 <- changenames(data = GSE64810, anno = cbind(rownames(GSE64810), ids[match(gsub('\\.[0-9]', '', rownames(GSE64810)), ids[,1]), 2]))
GSE64810 <- log2(GSE64810 + 1)

# GSE102741: log2 RPKM
GSE102741 <- read.xlsx(g("GSE102741_log2RPKMcounts.xlsx"), 1, rowNames = T)
ids <- ens_symbols(rownames(GSE102741))
GSE102741 <- changenames(data = GSE102741, anno = cbind(rownames(GSE102741), ids[match(rownames(GSE102741), ids[,1]), 2]))

# GSE174409: raw counts -> log2 CPM
GSE174409 <- read.csv(g("GSE174409_raw_counts_02102020.csv.gz"), row.names = 1)
GSE174409 <- log2(t(t(GSE174409) / colSums(GSE174409)) * 1000000 + 1)
ids <- ens_symbols(rownames(GSE174409))
GSE174409 <- changenames(data = GSE174409, anno = cbind(rownames(GSE174409), ids[match(rownames(GSE174409), ids[,1]), 2]))

#####################################################################
## 2. Affymetrix 3' arrays: raw CEL files, affy::rma,
##    gene symbols from the annotation of a GEO series
#####################################################################
affy_sets <- list(
  GSE37981  = c(cel = "GSE37981/GSE37981_RAW",   anno = "GSE37981_series_matrix.txt.gz"),
  GSE11512  = c(cel = "GSE11512/GSE11512_RAW",   anno = "GSE5392_series_matrix.txt.gz"),
  GSE13564  = c(cel = "GSE13564_RAW",            anno = "GSE13564_series_matrix.txt.gz"),
  GSE161986 = c(cel = "GSE161986/GSE161986_RAW", anno = "GSE161986_series_matrix.txt.gz"),
  GSE17612  = c(cel = "GSE17612_RAW",            anno = "GSE161986_series_matrix.txt.gz"),
  GSE21138  = c(cel = "GSE21138_RAW",            anno = "GSE21138_series_matrix.txt.gz"),
  GSE21935  = c(cel = "GSE21935/GSE21935_RAW",   anno = "GSE21935_series_matrix.txt.gz"),
  GSE22570  = c(cel = "GSE22570_RAW",            anno = "GSE22570_series_matrix.txt.gz"),
  GSE5388   = c(cel = "GSE5388/GSE5388_RAW",     anno = "GSE5388_series_matrix.txt.gz"),
  GSE5390   = c(cel = "GSE5390/GSE5390_RAW",     anno = "GSE5390_series_matrix.txt.gz")
)
for(d in names(affy_sets)){
  raw  <- ReadAffy(celfile.path = paste0(g(affy_sets[[d]]["cel"]), "/"))
  anno <- series_anno(affy_sets[[d]]["anno"])
  assign(d, changenames(data = exprs(affy::rma(raw)),
                        anno = cbind(as.character(anno$ID), anno$`Gene symbol`)))
}

#####################################################################
## 3. Affymetrix gene/exon arrays and HG-U133A subset: oligo::rma
#####################################################################
suppressMessages(library(oligo))

# GSE59630 (Human Exon 1.0 ST): symbols from the gene_assignment field
raw  <- read.celfiles(list.files(g("GSE59630_RAW"), full.names = TRUE))
anno <- series_anno("GSE59630_series_matrix.txt.gz", annot = FALSE)
symbols <- gsub(' ', '', gsub('([^/]*)//([^/]*)//.*', '\\2', anno$gene_assignment))
GSE59630 <- changenames(data = exprs(oligo::rma(raw)), anno = cbind(as.character(anno$ID), symbols))

# GSE71620 (Human Gene 1.1 ST)
raw  <- read.celfiles(list.files(g("GSE71620_RAW"), full.names = TRUE))
anno <- series_anno("GSE71620_series_matrix.txt.gz")
GSE71620 <- changenames(data = exprs(oligo::rma(raw)), anno = cbind(as.character(anno$ID), anno$`Gene symbol`))

# GSE92538, HG-U133A arrays 
cel  <- list.files(g("GSE92538_RAW"), pattern = "cel")
cel  <- cel[grep("_133A_", cel)]
raw  <- read.celfiles(g("GSE92538_RAW", cel))
anno <- series_anno("GSE5392_series_matrix.txt.gz")
GSE92538_GPL17027 <- changenames(data = exprs(oligo::rma(raw)), anno = cbind(as.character(anno$ID), anno$`Gene symbol`))

#####################################################################
## 4. Illumina arrays: neqc normalisation
#####################################################################
# GSE60190: non-normalised intensities + detection p-values
raw <- read.csv(g("GSE60190", "GSE60190_non-normalized.txt"), sep = "\t", row.names = 1)
GSE60190 <- neqc(as.matrix(raw[, seq(1, 266, 2)]), detection.p = as.matrix(raw[, seq(2, 266, 2)]))
anno <- series_anno("GSE60190_series_matrix.txt.gz")
GSE60190 <- changenames(data = GSE60190, anno = cbind(as.character(anno$ID), anno$`Gene symbol`))

# GSE49376: rows are already gene symbols; rows with NA removed 
raw <- read.csv(g("GSE49376", "GSE49376_non-normalized.txt"), sep = "\t", row.names = 1)
GSE49376 <- neqc(as.matrix(raw[, seq(1, 95, 2)]), detection.p = as.matrix(raw[, seq(2, 96, 2)]))
GSE49376 <- GSE49376[rowSums(is.na(GSE49376)) == 0,]

# GSE125681: idat files, kept as a limma EList with probe IDs 
#(no sample of this dataset is "Healthy", so it is not used in the
# aging analyses). Column names are the idat file paths.
idat <- list.files(g("GSE125681_RAW"), pattern = "idat$")
GSE125681 <- neqc(read.idat(g("GSE125681_RAW", idat),
                            bgxfile = g("GSE125681_RAW", "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz")))

#####################################################################
## 5. Processed supplementary files
#####################################################################
# GSE30272: cleaned expression matrix, custom platform annotation
GSE30272 <- read.csv(g("GSE30272_ExprsMtxCleanedN269_31SVN.txt.gz"), sep = "\t", row.names = 1)
anno <- read.csv(g("GPL4611_GSE30272.txt"), sep = "\t", skip = 11)
GSE30272 <- changenames(data = GSE30272, anno = cbind(anno$ID, anno$Gene_Symbol))
meta <- pData(getGEO(filename = g("GSE30272_series_matrix.txt.gz"), getGPL = FALSE))
colnames(GSE30272) <- meta[match(colnames(GSE30272), meta[,1]), 2]

# GSE102556 (human samples): FPKM, gene symbol in column 282
GSE102556_H <- read.csv(g("GSE102556_HumanMDD_fpkmtab.txt.gz"), sep = "\t", row.names = 1)
GSE102556_H <- changenames(data = GSE102556_H[, c(1:281)], anno = cbind(rownames(GSE102556_H), GSE102556_H[, 282]))
GSE102556_H <- log2(GSE102556_H + 1)

# GSE99349: gene symbol in column 1; sample names from the series matrix
GSE99349 <- read.xlsx(g("GSE99349_Expression_Matrix.xlsx"), 1, rowNames = T)
GSE99349 <- changenames(data = GSE99349[, -c(1:6)], anno = cbind(rownames(GSE99349), GSE99349[, 1]))
meta <- pData(getGEO(filename = g("GSE99349_series_matrix.txt.gz"), getGPL = FALSE))
meta <- meta[-which(meta$title == "Neun DLPFC_Total_RNAseq_28"),]
colnames(GSE99349) <- meta[, 2]

#####################################################################
## Save
#####################################################################
datasets <- c("GSE101521", "GSE102556_H", "GSE102741", "GSE11512", "GSE125681", "GSE13564",
              "GSE161986", "GSE174409", "GSE17612", "GSE21138", "GSE21935", "GSE22570",
              "GSE30272", "GSE37981", "GSE49376", "GSE5388", "GSE5390", "GSE59288_H",
              "GSE59630", "GSE60190", "GSE64810", "GSE71620", "GSE92538_GPL17027", "GSE99349")
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
save(list = datasets, file = OUT_FILE)
