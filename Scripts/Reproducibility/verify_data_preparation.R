################################################################################
# verify_data_preparation.R - compares the four RData files rebuilt from the raw
# data (<run_dir>/rebuilt_RData/*.RData) with the deposited ones (<ref_dir>/*.RData:
# Data/original_RData/ or Data/original_RData_notclosed/), which are the files the analysis uses
#
# For every object of the regenerated file: same row names, same column names,
# largest absolute difference of the numeric values (matrices, data.frames,
# DESeq2 results, limma ELists); data.frames of metadata are compared exactly.
# The rebuilt files hold only DE tables and expression matrices: Data/original_RData_notclosed/
# lacks only the raw MDMA and DMT tables, which are never rebuilt, so both references
# give the same comparison. Without the non-public raw data, DE_MDMA, DE_DMT, DE_pharm,
# DE_harm are not rebuilt and are listed as such.
# Writes <run_dir>/comparison/data_preparation_verification.csv.
# Run from the folder root:
#   Rscript Scripts/Reproducibility/verify_data_preparation.R <run_dir> [ref_dir, default Data/original_RData]
################################################################################
suppressPackageStartupMessages({ library(limma); library(DESeq2) })
run_dir <- commandArgs(trailingOnly = TRUE)[1]
ref_dir <- commandArgs(trailingOnly = TRUE)[2]
if (is.na(ref_dir)) ref_dir <- "Data/original_RData"
stopifnot(!is.na(run_dir), dir.exists(file.path(run_dir, "rebuilt_RData")), dir.exists(ref_dir))
message("reference: ", ref_dir)
closed_de <- c("DE_MDMA", "DE_DMT", "DE_pharm", "DE_harm")
dir.create(file.path(run_dir, "comparison"), showWarnings = FALSE)

as_tab <- function(x) {
  if (is(x, "EList")) return(x$E)
  if (is(x, "DataFrame")) return(as.data.frame(x))
  x
}
max_num_diff <- function(a, b) {
  a <- as.data.frame(a); b <- as.data.frame(b)
  num <- intersect(names(a)[sapply(a, is.numeric)], names(b)[sapply(b, is.numeric)])
  if (!length(num)) return(NA)
  rn <- intersect(rownames(a), rownames(b))
  d <- unlist(lapply(num, function(k) abs(a[rn, k] - b[rn, k])))
  if (all(is.na(d))) return(0)
  max(d, na.rm = TRUE)
}

files <- c("Alldata_20Sep.RData", "Psychedelics_PFC.RData", "ExerciseAndAlcohol.RData", "Dementia_alldata_meta.RData")
res <- list()
for (f in files) {
  new <- file.path(run_dir, "rebuilt_RData", f)
  if (!file.exists(new)) { message("missing ", new); next }
  ref <- file.path(ref_dir, f)
  if (!file.exists(ref)) { message("no ", ref, ": not compared"); next }
  o <- new.env(); load(ref, o)
  n <- new.env(); load(new, n)
  for (k in setdiff(intersect(closed_de, ls(o)), ls(n)))
    res[[length(res) + 1]] <- data.frame(file = f, object = k,
      note = "not rebuilt: non-public raw data not available")
  for (k in ls(n)) {
    if (!exists(k, o, inherits = FALSE)) { res[[length(res) + 1]] <- data.frame(file = f, object = k, note = "not in original"); next }
    a <- as_tab(get(k, o)); b <- as_tab(get(k, n))
    meta <- is.data.frame(a) && k == "metadata"
    if (!meta && !(is.matrix(a) || is.data.frame(a))) {   # not a table: compare as a whole
      res[[length(res) + 1]] <- data.frame(file = f, object = k,
        identical_values = isTRUE(all.equal(a, b, check.attributes = FALSE)), note = "not a table")
      next
    }
    res[[length(res) + 1]] <- tryCatch(data.frame(
      file = f, object = k,
      dim_original = paste(dim(a), collapse = "x"), dim_new = paste(dim(b), collapse = "x"),
      same_rownames = identical(rownames(a), rownames(b)),
      same_colnames = identical(colnames(a), colnames(b)),
      max_abs_diff = if (meta) NA else signif(max_num_diff(a, b), 3),
      identical_values = if (meta) { ra <- a; rb <- b; rownames(ra) <- rownames(rb) <- NULL; isTRUE(all.equal(ra, rb)) } else NA,
      note = ""), error = function(e) data.frame(file = f, object = k, note = paste("not compared:", conditionMessage(e))))
  }
  rm(o, n); invisible(gc())
}
res <- do.call(rbind, lapply(res, function(r) { for (c in c("dim_original","dim_new","same_rownames","same_colnames","max_abs_diff","identical_values","note")) if (is.null(r[[c]])) r[[c]] <- NA; r[, c("file","object","dim_original","dim_new","same_rownames","same_colnames","max_abs_diff","identical_values","note")] }))
print(res, row.names = FALSE)
write.csv(res, file.path(run_dir, "comparison", "data_preparation_verification.csv"), row.names = FALSE)
