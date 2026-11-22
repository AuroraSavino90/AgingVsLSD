################################################################################
# install_environment.R - installs, OFFLINE, every R package needed to reproduce
# the paper, from the files in packages/. The user's R library is not changed.
#
#   Rscript Scripts/Reproducibility/install_environment.R            # into <USERPROFILE>/AgingVsLSD_Rlib
#   Rscript Scripts/Reproducibility/install_environment.R <root>     # into another folder (keep it outside Dropbox)
#
# R 4.5 on Windows. Four libraries are built in <root> (about 2 GB):
#   base             Bioconductor 3.21 software and annotations
#                    (packages/offline_repo_R4.5_windows + packages/bioc3.21)
#   prep             org.Mm.eg.db 3.14.0, mouse4302.db 3.13.0, array-reading packages:
#                    put in front of base for the data preparation (DE_MDMA needs them)
#   ann322           AnnotationDbi 1.72.0 + org.Hs/Mm/Rn.eg.db and GO.db 3.22.0
#   ann322_orgMm321  as ann322 but org.Mm.eg.db 3.21.0
#   ggplot2_35       ggplot2 3.5.2 (packages/cran_archive) + ggrepel 0.9.6 built against
#                    ggplot2 3.5 (packages/cran_archive/built_for_ggplot2_3.5): the figures of
#                    01, 05 and 12 were drawn with ggplot2 3.5 (theme: grey axis text; Fig S1A:
#                    histogram bins). The ggrepel 0.9.6 binary of the offline repo was built
#                    against ggplot2 4 and fails with ggplot2 3.5 ("could not find gg_par").
#   factominer212    FactoMineR 2.12 (packages/cran_archive): FactoMineR 2.18 draws random
#                    numbers in PCA(), which shifts the random state before the GSEA of Fig 1A (04)
# Which script uses which library: Data/envs/environments.csv (each overlay is put in
# front of base).
################################################################################
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else file.path(Sys.getenv("USERPROFILE"), "AgingVsLSD_Rlib")
lib_names <- c("base", "prep", "ann322", "ann322_orgMm321", "ggplot2_35", "factominer212")
libs <- setNames(file.path(root, lib_names), lib_names)
for (l in libs) dir.create(l, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type != "windows" || R.version$major != "4" || !startsWith(R.version$minor, "5"))
  warning("packages/ holds Windows binaries for R 4.5; this is ", R.version.string)
has <- function(p, lib) nzchar(system.file(package = p, lib.loc = lib))
src <- function(f, lib) install.packages(f, lib = lib, repos = NULL, type = "source")

# 1. base -------------------------------------------------------------------
.libPaths(c(libs[["base"]], .Library))
repo <- paste0("file:///", normalizePath("packages/offline_repo_R4.5_windows", winslash = "/"))
pkgs <- read.csv("packages/package_versions.csv")$Package
pkgs <- setdiff(pkgs, rownames(installed.packages(lib.loc = c(libs[["base"]], .Library))))
if (length(pkgs)) install.packages(pkgs, lib = libs[["base"]], repos = repo, type = "binary")
for (f in list.files("packages/bioc3.21", "tar.gz$", full.names = TRUE))
  if (!has(sub("_.*", "", basename(f)), libs[["base"]])) src(f, libs[["base"]])

# 2. prep (data preparation) --------------------------------------------------
lib <- libs[["prep"]]
.libPaths(c(lib, libs[["base"]], .Library))
if (!has("org.Mm.eg.db", lib)) src("packages/org.Mm.eg.db_3.14.0.tar.gz", lib)
if (!has("mouse4302.db", lib)) src("packages/mouse4302.db_3.13.0.tar.gz", lib)
for (f in list.files("packages/array_software/source", full.names = TRUE)) {
  p <- sub("_.*$", "", basename(f))
  if (has(p, c(lib, libs[["base"]]))) next
  # without Rtools (no make): Windows binary build if available
  bin <- list.files("packages/array_software/windows_binary", paste0("^", p, "_"), full.names = TRUE)
  if (.Platform$OS.type == "windows" && !nzchar(Sys.which("make")) && length(bin)) {
    install.packages(bin, lib = lib, repos = NULL, type = "binary")
  } else src(f, lib)
}
for (f in list.files("packages/array_annotation", full.names = TRUE))
  if (!has(sub("_.*$", "", basename(f)), c(lib, libs[["base"]]))) src(f, lib)

# 3. annotation overlays for the analysis ---------------------------------------
# The 3.22 annotation packages need AnnotationDbi >= 1.71.1 (Bioconductor 3.22).
for (ov in c("ann322", "ann322_orgMm321")) {
  lib <- libs[[ov]]
  .libPaths(c(lib, libs[["base"]], .Library))
  if (!has("AnnotationDbi", lib))
    install.packages("packages/bioc3.22/AnnotationDbi_1.72.0.zip", lib = lib, repos = NULL, type = "binary")
  ann <- list.files("packages/bioc3.22", "tar.gz$", full.names = TRUE)
  if (ov == "ann322_orgMm321")
    ann <- c(ann[!grepl("^org.Mm", basename(ann))], "packages/bioc3.21/org.Mm.eg.db_3.21.0.tar.gz")
  for (f in ann) if (!has(sub("_.*", "", basename(f)), lib)) src(f, lib)
}

# 4. ggplot2 3.5.2 (pure R: no Rtools needed) + ggrepel built against it, for 01, 05, 12 ---
lib <- libs[["ggplot2_35"]]
.libPaths(c(lib, libs[["base"]], .Library))
if (!has("ggplot2", lib)) src("packages/cran_archive/ggplot2_3.5.2.tar.gz", lib)
if (!has("ggrepel", lib))
  install.packages("packages/cran_archive/built_for_ggplot2_3.5/ggrepel_0.9.6.zip", lib = lib, repos = NULL, type = "binary")

# 5. FactoMineR 2.12 for 04_integrated_dataset.R ---------------------------------------
lib <- libs[["factominer212"]]
.libPaths(c(lib, libs[["base"]], .Library))
if (!has("FactoMineR", lib))
  install.packages("packages/cran_archive/FactoMineR_2.12.zip", lib = lib, repos = NULL, type = "binary")

for (n in names(libs)) {
  ipm <- installed.packages(lib.loc = libs[[n]]); ip <- setNames(ipm[, "Version"], ipm[, "Package"])
  key <- grepl("^(AnnotationDbi|GO.db|org\\..*|mouse4302.db|fgsea|clusterProfiler|ggplot2|ggrepel|FactoMineR|oligo|affy)$", names(ip))
  cat(sprintf("%-16s %3d packages | %s\n", n, length(ip), paste(names(ip)[key], ip[key], collapse = ", ")))
}
cat("Libraries in", root, "\n")
