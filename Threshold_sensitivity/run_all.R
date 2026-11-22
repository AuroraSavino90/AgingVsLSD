# ============================================================
# Threshold sensitivity: APR (aging, dementia) and CRB (Fig. 1E-F)
# with the SAME gene-selection threshold for every intervention dataset.
#
# Run from the repository root (AgingVsLSD/), not from inside Threshold_sensitivity/.
# Inputs: Data/original_RData/ (the four complete input RData),
# Data/intermediate_RData/DE_GSE179379.RData (intermediate file of 00_LSD_DEGs.R),
# Data/raw_data_notclosed/ (GSE75772), Data/homologs/, Data/ensembl/, Scripts/utilities.R,
# Scripts/data_preparation/ensembl_release_map.R, Data/original_RData_notclosed/rng_states.rds.
# The main-pipeline Results/ is not needed.
#   Rscript Threshold_sensitivity/run_all.R            # both modes
#   Rscript Threshold_sensitivity/run_all.R pvalue     # nominal p < 0.05
#   Rscript Threshold_sensitivity/run_all.R fdr        # FDR < 0.05
#
# Outputs: Threshold_sensitivity/Results_pvalue, Threshold_sensitivity/Results_FDR
# (main figures and tables are collected in <Results_*>/Main_figures).
#
# R libraries: if REPRO_LIBS is set (the libraries built by install_environment.R),
# each script uses the same library as the main-pipeline script it is a copy of
# (Data/envs/environments.csv). Otherwise the current R library is used.
#
# Data: for now this analysis needs the complete RData, non-public data included
# (DE_MDMA, DE_DMT, DE_pharm, DE_harm in Data/original_RData/Psychedelics_PFC.RData). It stops
# if they are missing.
# ============================================================

args  <- commandArgs(trailingOnly = TRUE)
modes <- if (length(args)) tolower(args) else c("pvalue", "fdr")
stopifnot(all(modes %in% c("pvalue", "fdr")), file.exists("Data"),
          file.exists("Data/intermediate_RData/DE_GSE179379.RData"))

rscript <- file.path(R.home("bin"), "Rscript")
here    <- "Threshold_sensitivity"
scripts <- c(
  "03_Age_reversal.R"                   = "base",        # LSD (chronic) GSEA vs aging
  "09_age_reversal_otherpsychedelics.R" = "base",        # other interventions + rat aging; APR aging
  "12_dementia_cohen.R"                 = "ggplot2_35",  # APR dementia
  "save_mean_logFC_selected_genes.R"    = "ann322",      # input for CRB
  "Cumulative_logFC.R"                  = "base"         # CRB panels E-F
)

# ---- libraries -------------------------------------------------------------------------
libroot <- Sys.getenv("REPRO_LIBS")
libpath <- function(e) paste(c(if (e != "base") file.path(libroot, e), file.path(libroot, "base")),
                             collapse = .Platform$path.sep)

# ---- complete data required ------------------------------------------------------------
closed <- c("DE_MDMA", "DE_DMT", "DE_pharm", "DE_harm")
psy <- new.env(); load("Data/original_RData/Psychedelics_PFC.RData", envir = psy)
missing <- setdiff(closed, ls(psy)); rm(psy)
if (length(missing))
  stop("Threshold_sensitivity needs the complete RData; missing in Data/original_RData/Psychedelics_PFC.RData: ",
       paste(missing, collapse = ", "))
todo <- scripts

for (mode in modes) {
  label <- if (mode == "fdr") "FDR" else "pvalue"
  out   <- file.path(here, paste0("Results_", label))
  dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(SELECTION_MODE = mode, SENS_OUT = out)

  for (s in names(todo)) {
    log <- file.path(out, "logs", sub("\\.R$", ".log", s))
    message(format(Sys.time(), "%H:%M"), "  [", label, "] ", s)
    if (nzchar(libroot)) { old <- Sys.getenv("R_LIBS_USER"); Sys.setenv(R_LIBS_USER = libpath(todo[[s]])) }
    status <- system2(rscript, shQuote(file.path(here, "Scripts", s)), stdout = log, stderr = log)
    if (nzchar(libroot)) Sys.setenv(R_LIBS_USER = old)
    if (status != 0) stop(s, " failed (", label, "); see ", log)
  }

  # collect the main figures / tables under explicit names
  main <- file.path(out, "Main_figures")
  dir.create(main, showWarnings = FALSE)
  copy <- c(
    "Figures/All_psychedelicsAndCTRL_rataging.pdf" = paste0("APR_aging_", label, ".pdf"),
    "Figures/All_psychedelicsAndCTRL_dementia.pdf" = paste0("APR_dementia_", label, ".pdf"),
    "aging_opposition_scores.xlsx"                 = paste0("APR_aging_scores_", label, ".xlsx"),
    "dementia_opposition_scores.xlsx"              = paste0("APR_dementia_scores_", label, ".xlsx"),
    setNames(paste0(c("CRB_panel_E_aging_", "CRB_panel_F_dementia_"), label, ".pdf"),
             paste0(c("CRB_panel_E_aging_", "CRB_panel_F_dementia_"), label, ".pdf")),
    setNames(paste0("CRB_values_panels_E_F_", label, ".xlsx"), paste0("CRB_values_panels_E_F_", label, ".xlsx"))
  )
  ok <- file.copy(file.path(out, names(copy)), file.path(main, copy), overwrite = TRUE)
  if (!all(ok)) warning("missing outputs: ", paste(names(copy)[!ok], collapse = ", "))
  message("[", label, "] done -> ", main)
}
