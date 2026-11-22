# ============================================================
# Run the analysis scripts of the paper in <run_dir>, in dependency order, each
# in a fresh R session with its own R library (Data/envs/environments.csv).
# Normally called by run_all.R, which prepares <run_dir>/Data first.
#
# From the folder root:
#   Rscript Scripts/Reproducibility/run_pipeline.R <run_dir> [script ...]
#
# <run_dir>/Data must link the folders of Data/ (run_all.R does it); the scripts read the
# input RData from Data/original_RData_notclosed/.
# This script adds to <run_dir>:
#   Scripts/                      copy of Scripts/ (the scripts read the random states from
#                                 Data/original_RData_notclosed/rng_states.rds)
#   Results/RData/                CMAP_NPC_mean.RData, compounds_NPC.RData (from Data/original_RData_notclosed/cmap)
#   logs/                         .log, .status, .sessionInfo.txt, _objects.RData per script
# Libraries: REPRO_LIBS (default <USERPROFILE>/AgingVsLSD_Rlib), built by Scripts/Reproducibility/install_environment.R.
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1, dir.exists("Scripts"), file.exists("Data/envs/environments.csv"))
top     <- normalizePath(".", winslash = "/")
run_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
only    <- if (length(args) >= 2) args[-1] else character(0)
need <- file.path(run_dir, "Data", "original_RData_notclosed", c("Alldata_20Sep.RData", "Psychedelics_PFC.RData",
                                    "ExerciseAndAlcohol.RData", "Dementia_alldata_meta.RData"))
if (!all(file.exists(need))) stop("missing input RData in ", run_dir, "/Data/original_RData_notclosed: run run_all.R")

layers <- list(
  L0 = c("00_LSD_DEGs.R", "00_GTEx_filt.R", "16_invitro_RNAseq.R"),
  L1 = c("01_aging_metanalysis.R", "03_Age_reversal.R", "04_integrated_dataset.R",
         "05_DEvsAging_validation.R", "06_GO_aging.R", "09_age_reversal_otherpsychedelics.R",
         "11_validation_datasets.R", "12_dementia_cohen.R", "13_AD_validation.R",
         "save_mean_logFC_selected_genes.R"),
  L2 = c("08_2_random_drugs.R", "14_Signature_agreement.R", "17_AD_DEGs.R", "18_GO_allLSD.R",
         "19_GO_revert_dementia.R", "20_NES_figures.R", "21_cohen_figures.R", "Cumulative_logFC.R"))

# ---- non-public data -----------------------------------------------------------
# Psychedelics_PFC.RData must hold the DE tables of the non-public datasets (DE_MDMA,
# DE_DMT, DE_pharm, DE_harm); Data/original_RData_notclosed/Psychedelics_PFC.RData has them
# (README, "Non-public data").
# Without the MSBB data (Data/raw_data_closed/AD portal/, controlled access), 13_AD_validation.R reads
# the gene correlations from Data/closed_data_intermediate/MSBB_gene_correlations.RData.
psy <- new.env(); load(file.path(run_dir, "Data", "original_RData_notclosed", "Psychedelics_PFC.RData"), envir = psy)
miss <- setdiff(c("DE_MDMA", "DE_DMT", "DE_pharm", "DE_harm"), ls(psy)); rm(psy)
if (length(miss)) stop(run_dir, "/Data/original_RData_notclosed/Psychedelics_PFC.RData has no ", paste(miss, collapse = ", "))

# ---- working copies ------------------------------------------------------------
unlink(file.path(run_dir, "Scripts"), recursive = TRUE)
file.copy("Scripts", run_dir, recursive = TRUE, copy.date = TRUE)
dir.create(file.path(run_dir, "Results", "RData"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "Results", "Figures"), showWarnings = FALSE)
file.copy(list.files("Data/original_RData_notclosed/cmap", full.names = TRUE), file.path(run_dir, "Results", "RData"), overwrite = TRUE)
dir.create(file.path(run_dir, "logs"), showWarnings = FALSE)
writeLines(format(Sys.time()), file.path(run_dir, "logs", "pipeline_start.txt"))

# ---- libraries ---------------------------------------------------------------------
envs <- read.csv("Data/envs/environments.csv", colClasses = "character")
root <- Sys.getenv("REPRO_LIBS", file.path(Sys.getenv("USERPROFILE"), "AgingVsLSD_Rlib"))
libpath <- function(e) paste(c(if (e != "base") file.path(root, e), file.path(root, "base")),
                             collapse = .Platform$path.sep)
for (l in file.path(root, unique(c("base", envs$environment))))
  if (!dir.exists(l)) stop("library ", l, " not found: run Rscript Scripts/Reproducibility/install_environment.R first")
old_libs <- Sys.getenv("R_LIBS_USER")

rscript <- file.path(R.home("bin"), "Rscript")
runner  <- file.path(top, "Scripts", "Reproducibility", "run_script.R")
for (L in names(layers)) {
  todo <- layers[[L]]
  if (length(only)) todo <- intersect(todo, only)
  if (!length(todo)) next
  message(format(Sys.time(), "%H:%M"), "  ", L, ": ", paste(todo, collapse = ", "))
  st <- file.path(run_dir, "logs", sub("\\.R$", ".status", todo))
  unlink(st)
  for (s in todo) {
    e <- envs$environment[envs$script == s]
    if (length(e) != 1) stop("no environment for ", s, " in environments.csv")
    Sys.setenv(R_LIBS_USER = libpath(e))
    system2(rscript, c(shQuote(runner), shQuote(run_dir), shQuote(s)),
            stdout = file.path(run_dir, "logs", sub("\\.R$", ".log", s)),
            stderr = file.path(run_dir, "logs", sub("\\.R$", ".log", s)), wait = FALSE)
    Sys.setenv(R_LIBS_USER = old_libs)
  }
  repeat { Sys.sleep(20); if (all(file.exists(st))) break }
  for (i in seq_along(todo))
    message("   ", format(todo[i], width = 40), readLines(st[i])[1], "  (",
            sub("minutes: ", "", readLines(st[i])[2]), " min)")
}
message(format(Sys.time(), "%H:%M"), "  analysis finished")
