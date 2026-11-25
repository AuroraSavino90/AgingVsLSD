################################################################################
# make_notclosed_RData.R - builds Data/original_RData_notclosed/ from Data/original_RData/.
# Accessory script, NOT part of the pipeline: it was run once to make the input
# RData that can be distributed (Zenodo).
#
#   Rscript Scripts/Reproducibility/make_notclosed_RData.R
#
# - Alldata_20Sep.RData, ExerciseAndAlcohol.RData, Dementia_alldata_meta.RData:
#   no non-public data, copied unchanged (same bytes).
# - Psychedelics_PFC.RData: the original workspace without the two tables of raw
#   non-public data, MDMA (MDMA arrays, RMA expression matrix) and DMT (read counts of
#   DMT, pharmahuasca, harmaline). Every other object is kept as it is, including the
#   DE tables DE_MDMA, DE_DMT, DE_pharm, DE_harm and the random state .Random.seed.
# After saving, the new file is reloaded and every object is checked to have the same
# serialized content as the original one.
# Data/original_RData/ is only read. The four RData must not yet exist in Data/original_RData_notclosed/.
#
#   Rscript Scripts/Reproducibility/make_notclosed_RData.R corrected
#
# Corrected pipeline (Script_corrected/): copies the three *_corrected.RData written by
# Script_corrected/make_corrected_RData.R (Alldata_20Sep, ExerciseAndAlcohol,
# Dementia_alldata_meta; no non-public data) unchanged into Data/original_RData_notclosed/.
################################################################################
src <- "Data/original_RData"; dest <- "Data/original_RData_notclosed"
corrected <- identical(commandArgs(trailingOnly = TRUE), "corrected")
if (corrected) {
  files <- c("Alldata_20Sep_corrected.RData", "ExerciseAndAlcohol_corrected.RData",
             "Dementia_alldata_meta_corrected.RData")
  stopifnot(all(file.exists(file.path(src, files))))
  if (any(file.exists(file.path(dest, files)))) stop(dest, " already has the corrected RData")
  for (f in files) {
    stopifnot(file.copy(file.path(src, f), file.path(dest, f), copy.date = TRUE))
    stopifnot(unname(tools::md5sum(file.path(src, f)) == tools::md5sum(file.path(dest, f))))
    message(f, ": copied unchanged")
  }
  quit(save = "no")
}
files  <- c("Alldata_20Sep.RData", "ExerciseAndAlcohol.RData", "Dementia_alldata_meta.RData",
            "Psychedelics_PFC.RData")
closed <- c("MDMA", "DMT")   # raw data of the non-public datasets in Psychedelics_PFC.RData
stopifnot(all(file.exists(file.path(src, files))))
if (any(file.exists(file.path(dest, files))))
  stop(dest, " already has the input RData")
dir.create(dest, showWarnings = FALSE)

# open RData: copied unchanged
for (f in files[1:3]) {
  stopifnot(file.copy(file.path(src, f), file.path(dest, f), copy.date = TRUE))
  stopifnot(unname(tools::md5sum(file.path(src, f)) == tools::md5sum(file.path(dest, f))))
  message(f, ": copied unchanged")
}

# Psychedelics_PFC.RData without the raw non-public data
f <- "Psychedelics_PFC.RData"
o <- new.env(); load(file.path(src, f), envir = o)
stopifnot(all(closed %in% ls(o)),
          all(c("DE_MDMA", "DE_DMT", "DE_pharm", "DE_harm", ".Random.seed") %in% ls(o, all.names = TRUE)))
keep <- setdiff(ls(o, all.names = TRUE), closed)
save(list = keep, envir = o, file = file.path(dest, f))

n <- new.env(); load(file.path(dest, f), envir = n)
stopifnot(setequal(ls(n, all.names = TRUE), keep))
# serialized content: identical() would compare environments (functions, AffyBatch,
# DESeqDataSet, ...) by address, not by content
same <- vapply(keep, function(k) identical(serialize(get(k, o), NULL), serialize(get(k, n), NULL)), logical(1))
if (!all(same)) stop("objects not identical after saving: ", paste(keep[!same], collapse = ", "))
message(f, ": removed ", paste(closed, collapse = ", "), "; ", length(keep),
        " objects kept, all identical to the original (same serialized content, .Random.seed included)")
