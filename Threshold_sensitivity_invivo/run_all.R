# ============================================================
# Threshold sensitivity, in vivo LSD (acute, mouse PFC; Fig. 2B-D, Tables S15-S18):
# genes of each LSD dose selected with nominal p < 0.05 (as in the paper code)
# or with FDR < 0.05 (as stated in the paper text).
#
# Run from the repository root (AgingVsLSD/), not from inside Threshold_sensitivity_invivo/:
#   Rscript Threshold_sensitivity_invivo/run_all.R            # both modes
#   Rscript Threshold_sensitivity_invivo/run_all.R pvalue     # nominal p < 0.05
#   Rscript Threshold_sensitivity_invivo/run_all.R fdr        # FDR < 0.05
#
# R libraries: if REPRO_LIBS is set (the libraries built by install_environment.R),
# each script uses the same library as the main-pipeline script it is a copy of
# (Data/envs/environments.csv). Otherwise the current R library is used.
# Main-pipeline outputs come from Data/intermediate_RData/; Scripts/ is only read.
# ============================================================

args  <- commandArgs(trailingOnly = TRUE)
modes <- if (length(args)) tolower(args) else c("pvalue", "fdr")
stopifnot(all(modes %in% c("pvalue", "fdr")), file.exists("Data"))

rscript <- file.path(R.home("bin"), "Rscript")
scripts <- file.path("Threshold_sensitivity_invivo", "Scripts", c(
  "05_DEvsAging_validation.R",  # in vivo DEGs, GSEA vs aging, APR (Fig. 2B, S15-S16), GO per dose
  "18_GO_allLSD.R"              # GO shared by LSD and aging (Fig. 2C-D, S17-S18)
))
envs <- c("05_DEvsAging_validation.R" = "ggplot2_35", "18_GO_allLSD.R" = "base")
libroot <- Sys.getenv("REPRO_LIBS")
libpath <- function(e) paste(c(if (e != "base") file.path(libroot, e), file.path(libroot, "base")),
                             collapse = .Platform$path.sep)

for (mode in modes) {
  label <- if (mode == "fdr") "FDR" else "pvalue"
  out   <- file.path("Threshold_sensitivity_invivo", paste0("Results_", label))
  dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(SELECTION_MODE = mode, SENS_OUT = out)

  for (s in scripts) {
    log <- file.path(out, "logs", sub("[.]R$", ".log", basename(s)))
    message(format(Sys.time(), "%H:%M"), "  [", label, "] ", basename(s))
    if (nzchar(libroot)) { old <- Sys.getenv("R_LIBS_USER"); Sys.setenv(R_LIBS_USER = libpath(envs[[basename(s)]])) }
    status <- system2(rscript, shQuote(s), stdout = log, stderr = log)
    if (nzchar(libroot)) Sys.setenv(R_LIBS_USER = old)
    if (status != 0) stop(basename(s), " failed (", label, "); see ", log)
  }

  # collect the main figures / tables under explicit names
  main <- file.path(out, "Main_figures")
  dir.create(main, showWarnings = FALSE)
  copy <- c(
    "Figures/All_DE_aging.pdf"              = "Fig2B_APR_invivo",
    "Figures/GOup_top_allLSD_drugonly.pdf"  = "Fig2C_left_GOup_LSD",
    "Figures/GO_shared_up_top.pdf"          = "Fig2C_right_GOup_aging",
    "Figures/GOdn_top_allLSD_drugonly.pdf"  = "Fig2D_left_GOdn_LSD",
    "Figures/GO_shared_dn_top.pdf"          = "Fig2D_right_GOdn_aging",
    "all_GSEA_invivo.xlsx"                  = "TableS15_GSEA_invivo",
    "invivo_metap_all.xlsx"                 = "TableS16_APR_invivo",
    "GO_all_stats_up.xlsx"                  = "TableS17_GO_up",
    "GO_all_stats_dn.xlsx"                  = "TableS18_GO_dn"
  )
  dest <- paste0(copy, "_", label, ".", tools::file_ext(names(copy)))
  ok <- file.copy(file.path(out, names(copy)), file.path(main, dest), overwrite = TRUE)
  if (!all(ok)) warning("missing outputs: ", paste(names(copy)[!ok], collapse = ", "))
  message("[", label, "] done -> ", main)
}
