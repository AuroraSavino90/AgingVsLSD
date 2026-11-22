################################################################################
# check_ensembl_releases.R - shows that the Ensembl ID -> symbol tables left in
# the original R workspaces are identical to Ensembl 105 (human, aging data)
# and Ensembl 107 (mouse, rat).
# Run from the main folder:  Rscript Scripts/Reproducibility/check_ensembl_releases.R
# The release tables are rebuilt from the BioMart dumps in Data/ensembl/
# (downloaded from https://ftp.ensembl.org/pub/release-<N>/mysql/ensembl_mart_<N>/).
################################################################################
source("Scripts/data_preparation/ensembl_release_map.R")
d <- "Data/ensembl/original_workspace_tables"
checks <- list(
  c("Aging_alldata27Apr_ids_hsapiens_hgnc.tsv",        "105", "hsapiens"),
  c("Aging_alldata27Apr_anno_hsapiens_hgnc.tsv",       "105", "hsapiens"),
  c("Psychedelics_PFC_ids_mmusculus_mgi.tsv",          "107", "mmusculus"),
  c("Psychedelics_PFC_geni_coord_rnorvegicus_rgd.tsv", "107", "rnorvegicus"),
  c("ExerciseAndAlcohol_ids_mmusculus_mgi.tsv",        "107", "mmusculus"))
for (ck in checks) {
  o <- read.delim(file.path(d, ck[1]), colClasses = "character", na.strings = character(0))
  m <- ensembl_release_map(as.integer(ck[2]), ck[3], "Data/ensembl")
  m <- m[m[, 1] %in% o[, 1], ]
  p <- function(x) paste(x[, 1], x[, 2])
  first <- function(x) x[match(unique(o[, 1]), x[, 1]), 2]
  cat(sprintf("%-48s Ensembl %s: %6d pairs | only in original %d | only in release %d | first symbol differs %d\n",
              ck[1], ck[2], nrow(o), sum(!p(o) %in% p(m)), sum(!p(m) %in% p(o)), sum(first(o) != first(m))))
}
