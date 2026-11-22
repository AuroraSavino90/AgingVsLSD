################################################################################
# run_all.R - runs the analysis pipeline from the raw data, in one command.
#
#   Rscript run_all.R <run_dir> [step ...] [--rdata=full|notclosed]
#
# Run from this folder, after `Rscript Scripts/Reproducibility/install_environment.R`.
# <run_dir>: an empty or new folder OUTSIDE Dropbox (about 8 GB). Everything that
# is produced goes there; this folder is only read.
#
# Steps (default: all, in this order):
#   prepare   rebuild the four input RData files from the raw data into <run_dir>/rebuilt_RData/
#             (Scripts/data_preparation, steps 1a 1b 2 3 4; about 4 h)
#   verify    compare them with the deposited RData -> comparison/data_preparation_verification.csv
#   analysis  run the analysis scripts of the paper (Scripts/Reproducibility/run_pipeline.R; about 5 h)
#             on the deposited RData of Data/original_RData_notclosed/
#
# Raw data: Data/raw_data_notclosed/ (public data) and, if present, Data/raw_data_closed/
# (non-public MDMA and DMT/pharmahuasca/harmaline data). Without Data/raw_data_closed/,
# prepare builds Psychedelics_PFC.RData without DE_MDMA, DE_DMT, DE_pharm, DE_harm.
#
# Deposited input RData:
#   Data/original_RData/            the original files (Psychedelics_PFC.RData also has the raw
#                                   MDMA and DMT tables; not distributed)
#   Data/original_RData_notclosed/  the same files, Psychedelics_PFC.RData without the raw MDMA
#                                   and DMT tables (Scripts/Reproducibility/make_notclosed_RData.R), plus
#                                   cmap/, seeded_intermediate/, clusterProfiler_cache_2022.rds
# The analysis scripts read Data/original_RData_notclosed/ (no script uses the raw MDMA and
# DMT tables), never the rebuilt files. --rdata only chooses the files verify compares the
# rebuilt ones with (default full if Data/original_RData/ has the four files).
#
# Data/ layout (the scripts read these paths): original_RData/, original_RData_notclosed/,
# intermediate_RData/, raw_data_notclosed/ (public raw data and the other public inputs),
# raw_data_closed/ (non-public data, AD portal/ included), closed_data_intermediate/,
# homologs/, ensembl/. <run_dir>/Data/ links these folders (except original_RData/ and packages/).
#
# No internet connection is used: all R processes get a non-working proxy, so any
# attempt to download something fails visibly. Set ALLOW_NETWORK=1 to allow it.
# R libraries: REPRO_LIBS (default <USERPROFILE>/AgingVsLSD_Rlib).
################################################################################
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript run_all.R <run_dir> [prepare verify analysis] [--rdata=full|notclosed]")
stopifnot(dir.exists("Scripts/data_preparation"), dir.exists("Scripts/Reproducibility"))
run_dir <- args[1]
opt     <- grep("^--rdata=", args[-1], value = TRUE)
steps   <- setdiff(args[-1], opt)
if (!length(steps)) steps <- c("prepare", "verify", "analysis")
stopifnot(all(steps %in% c("prepare", "verify", "analysis")))
inputs  <- c("Alldata_20Sep.RData", "Psychedelics_PFC.RData", "ExerciseAndAlcohol.RData", "Dementia_alldata_meta.RData")
rdata   <- if (length(opt)) sub("^--rdata=", "", opt[1]) else
             if (all(file.exists(file.path("Data/original_RData", inputs)))) "full" else "notclosed"
stopifnot(rdata %in% c("full", "notclosed"))
inter_dir <- "Data/original_RData_notclosed"   # distributed intermediate files
rdata_dir <- c(full = "Data/original_RData", notclosed = inter_dir)[[rdata]]
if (!all(file.exists(file.path(rdata_dir, inputs)))) stop(rdata_dir, "/ does not have the four input RData files")
raw_open   <- "Data/raw_data_notclosed"
raw_closed <- "Data/raw_data_closed/psychedelics"
has_raw_closed <- dir.exists(file.path(raw_closed, "MDMA data")) && file.exists(file.path(raw_closed, "readcount copy 2.xlsx"))
if ("prepare" %in% steps) stopifnot(dir.exists(raw_open))
message("verify compares with: ", rdata_dir, "/ (--rdata=", rdata, "); non-public raw data: ",
        if (has_raw_closed) raw_closed else "not available")
if (grepl("dropbox", tolower(normalizePath(run_dir, mustWork = FALSE))))
  warning("the run folder is inside Dropbox: it will hold about 8 GB of results")

libroot <- Sys.getenv("REPRO_LIBS", file.path(Sys.getenv("USERPROFILE"), "AgingVsLSD_Rlib"))
lib <- function(...) paste(file.path(libroot, c(...)), collapse = .Platform$path.sep)
if (!dir.exists(file.path(libroot, "base"))) stop("R libraries not found in ", libroot, ": run Rscript Scripts/Reproducibility/install_environment.R")
if (Sys.getenv("ALLOW_NETWORK") != "1")
  Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
             HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", ftp_proxy = "http://127.0.0.1:9")
rscript <- file.path(R.home("bin"), "Rscript")

# ---- run folder: Data/ = links to the folders of Data/ ---------------------------------
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "comparison"), showWarnings = FALSE)
run_dir <- normalizePath(run_dir, winslash = "/")
dd <- file.path(run_dir, "Data"); dir.create(dd, showWarnings = FALSE)
for (f in setdiff(list.dirs("Data", recursive = FALSE), c("Data/original_RData", "Data/packages"))) {
  target <- file.path(dd, basename(f))
  if (file.exists(target)) next
  if (.Platform$OS.type == "windows")
    shell(sprintf('mklink /J "%s" "%s"', normalizePath(target, winslash = "\\", mustWork = FALSE),
                  normalizePath(f, winslash = "\\")), intern = TRUE)
  else file.symlink(normalizePath(f), target)
}

run <- function(script, args = character(), libs, log, env = character()) {
  message(format(Sys.time(), "%H:%M"), "  ", basename(script), " ", paste(args, collapse = " "))
  old <- Sys.getenv(c("R_LIBS_USER", names(env)), unset = NA)
  do.call(Sys.setenv, as.list(c(R_LIBS_USER = libs, env)))
  on.exit(for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, setNames(list(old[[n]]), n)))
  st <- system2(rscript, c(shQuote(script), shQuote(args)), stdout = log, stderr = log)
  if (st != 0) stop(basename(script), " failed: see ", log)
}

# ---- 1. data preparation --------------------------------------------------------------
if ("prepare" %in% steps) {
  prep <- function(s, env = character(), args = character(), log)
    run(file.path("Scripts/data_preparation", s), args, lib("prep", "base"),
        file.path(run_dir, "logs", log), env)
  inter <- file.path(run_dir, "prep_intermediate"); dir.create(inter, showWarnings = FALSE)
  rb <- file.path(run_dir, "rebuilt_RData"); dir.create(rb, showWarnings = FALSE)
  unlink(file.path(rb, inputs))
  prep("prep_01a_aging_expression.R", log = "prep_01a.log", env = c(
    GEO_DIR = file.path(raw_open, "aging_GEO"), OUT_FILE = file.path(inter, "aging_expression.RData"),
    ENSEMBL_CACHE = "Data/ensembl"))
  prep("prep_01b_aging_metadata_Alldata_20Sep.R", log = "prep_01b.log", env = c(
    EXPR_FILE = file.path(inter, "aging_expression.RData"), GEO_DIR = file.path(raw_open, "aging_GEO"),
    GTEX_COUNTS = file.path(raw_open, "gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz"),
    GTEX_SUBJ = file.path(raw_open, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
    GTEX_SAMP = file.path(raw_open, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
    ENSEMBL_CACHE = "Data/ensembl", DOWNLOAD_DIR = file.path(raw_open, "dementia_GEO"),
    OUT_FILE = file.path(rb, "Alldata_20Sep.RData")))
  prep("prep_02_Psychedelics_PFC.R", log = "prep_02.log", env = c(
    RAW_DIR = file.path(raw_open, "psychedelics_GEO"), RAW_CLOSED_DIR = raw_closed,
    OUT_FILE = file.path(rb, "Psychedelics_PFC.RData"),
    ENSEMBL_CACHE = "Data/ensembl"))
  prep("prep_03_ExerciseAndAlcohol.R", log = "prep_03.log", env = c(
    RAW_DIR = file.path(raw_open, "exercise_alcohol_GEO"), OUT_FILE = file.path(rb, "ExerciseAndAlcohol.RData"),
    ENSEMBL_CACHE = "Data/ensembl"))
  prep("prep_04_Dementia_alldata_meta.R", log = "prep_04.log",
       args = c(file.path(raw_open, "dementia_GEO"), file.path(rb, "Dementia_alldata_meta.RData")))
}

# ---- 2. comparison with the deposited RData ------------------------------------------
if ("verify" %in% steps)
  run("Scripts/Reproducibility/verify_data_preparation.R", c(run_dir, rdata_dir), lib("base"),
      file.path(run_dir, "logs", "verify_data_preparation.log"))

# ---- 3. analysis scripts of the paper --------------------------------------------------
if ("analysis" %in% steps) {
  # the deposited input RData (never the rebuilt ones)
  if (!all(file.exists(file.path(inter_dir, inputs)))) stop(inter_dir, "/ does not have the four input RData files")
  message("  analysis input RData: ", inter_dir, "/")
  Sys.setenv(REPRO_LIBS = libroot)
  run("Scripts/Reproducibility/run_pipeline.R", run_dir, lib("base"), file.path(run_dir, "logs", "run_pipeline.log"))
}

message(format(Sys.time(), "%H:%M"), "  done")
