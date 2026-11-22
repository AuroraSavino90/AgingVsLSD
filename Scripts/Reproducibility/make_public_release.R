################################################################################
# make_public_release.R - copies the repository folder to <dest> without the data that cannot
# be distributed, for the public repository.
#
#   Rscript Scripts/Reproducibility/make_public_release.R <dest>          # copy
#   Rscript Scripts/Reproducibility/make_public_release.R <dest> --list   # only list what would be copied/excluded
#
# Excluded:
#   - non-public data (README_public.md, section 4): Data/raw_data_closed/ (MDMA arrays,
#     DMT/pharmahuasca/harmaline read counts, AD portal/ = MSBB)
#   - Data/original_RData/ and Zenodo/ (their Psychedelics_PFC.RData contains the raw MDMA
#     and DMT tables), and any leftover copy Data/<input>.RData. The distributed input RData
#     are in Data/original_RData_notclosed/ (Scripts/Reproducibility/make_notclosed_RData.R)
#   - the tables of Data/closed_data_intermediate/ that the pipeline no longer reads
#     (MSBB_gene_correlations.RData, gene-level correlations without sample data, is
#     copied: 13_AD_validation.R reads it when MSBB is missing)
#   - Data/raw_data_notclosed/GSE70138_*.gctx, *.gctx.gz: LINCS data from which
#     Data/original_RData_notclosed/cmap/ was prepared (08_1_cMAP_data.R), not read by the pipeline (31 GB)
#   - internal material: .git/, _internal/, _discarded/, original_figures/, README.md (internal),
#     this script
# README_public.md is copied as README.md.
# <dest> must not exist or be empty. This folder is only read.
################################################################################
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript Scripts/Reproducibility/make_public_release.R <dest> [--list]")
dest <- args[1]; only_list <- "--list" %in% args
inputs <- c("Alldata_20Sep.RData", "Psychedelics_PFC.RData", "ExerciseAndAlcohol.RData", "Dementia_alldata_meta.RData")
stopifnot(file.exists("README_public.md"), file.exists("Data/closed_data_intermediate/MSBB_gene_correlations.RData"),
          all(file.exists(file.path("Data/original_RData_notclosed", inputs))))

exclude_dirs  <- c("Data/raw_data_closed", "Data/original_RData", "Zenodo", ".git",
                   "_internal", "_discarded", "original_figures")
unused_closed <- setdiff(list.files("Data/closed_data_intermediate"), "MSBB_gene_correlations.RData")
gctx          <- list.files("Data/raw_data_notclosed", "^GSE70138_.*[.]gctx([.]gz)?$")
exclude_files <- c(file.path("Data", inputs), file.path("Data/raw_data_notclosed", gctx),
                   file.path("Data/closed_data_intermediate", unused_closed),
                   "README.md", "README_public.md", "Scripts/Reproducibility/make_public_release.R")

files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
drop  <- files %in% exclude_files |
         Reduce(`|`, lapply(exclude_dirs, function(d) startsWith(files, paste0(d, "/"))))
keep  <- files[!drop]

size <- function(f) sum(file.size(f), na.rm = TRUE) / 2^30
message(sprintf("copy: %d files, %.2f GB; excluded: %d files, %.2f GB",
                length(keep), size(keep), sum(drop), size(files[drop])))
if (only_list) {
  cat("EXCLUDED\n"); writeLines(paste0("  ", files[drop]))
  quit(save = "no")
}

if (dir.exists(dest) && length(list.files(dest, all.files = TRUE, no.. = TRUE)))
  stop(dest, " exists and is not empty")
for (d in unique(dirname(keep))) dir.create(file.path(dest, d), recursive = TRUE, showWarnings = FALSE)
ok <- file.copy(keep, file.path(dest, keep), copy.date = TRUE)
if (!all(ok)) stop("not copied: ", paste(keep[!ok], collapse = ", "))
file.copy("README_public.md", file.path(dest, "README.md"), copy.date = TRUE)

# check: none of the excluded data is in <dest>
left <- c(file.exists(file.path(dest, exclude_files[exclude_files != "README.md"])),
          dir.exists(file.path(dest, exclude_dirs)))
if (any(left)) stop("excluded files found in ", dest)
message("public release written to ", normalizePath(dest))
