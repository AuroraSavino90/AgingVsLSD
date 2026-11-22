# ============================================================
# Common setup for the threshold-sensitivity copies.
# Sourced by every script in Threshold_sensitivity/Scripts in place of
# setwd("workdir"). Working directory must be the repository root.
#
# Environment variables (set by run_all.R):
#   SELECTION_MODE  "pvalue" -> nominal p < 0.05 for ALL intervention datasets
#                   "fdr"    -> BH-adjusted p < 0.05 for ALL intervention datasets
#   SENS_OUT        output folder (default Threshold_sensitivity/Results_<mode>)
# ============================================================

selection_mode <- tolower(Sys.getenv("SELECTION_MODE", "pvalue"))
stopifnot(selection_mode %in% c("pvalue", "fdr"))

# column used for gene selection: limma (arrays) / DESeq2 (RNA-seq)
P_ARRAY <- if (selection_mode == "fdr") "adj.P.Val" else "P.Value"
P_SEQ   <- if (selection_mode == "fdr") "padj"      else "pvalue"

OUT <- Sys.getenv("SENS_OUT",
                  file.path("Threshold_sensitivity",
                            paste0("Results_", if (selection_mode == "fdr") "FDR" else "pvalue")))
dir.create(file.path(OUT, "RData"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "Figures"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists("Data") || !file.exists("Data/intermediate_RData/DE_GSE179379.RData"))
  stop("Run from the repository root (AgingVsLSD/): needs Data/intermediate_RData/DE_GSE179379.RData.")

set.seed(42)
suppressMessages(library(DESeq2))  # DESeqResults objects need DESeq2/S4Vectors loaded (rownames() returns NULL otherwise), as in Scripts/

# ------------------------------------------------------------
# Shims. With FDR some interventions have no (or one-direction-only)
# selected genes. The original code assumes non-empty UP and DN sets;
# these wrappers return NA for an empty set instead of stopping.
# When both sets are non-empty they return exactly what the original
# functions return.
# ------------------------------------------------------------

# fgsea: always return one row per set, ordered DN, UP (as fgsea does)
fgsea <- function(pathways, stats, ...) {
  pw <- pathways[vapply(pathways, function(g) length(intersect(g, names(stats))) > 0, logical(1))]
  res <- if (length(pw)) fgsea::fgsea(pw, stats, ...) else NULL
  miss <- setdiff(c("DN", "UP"), if (is.null(res)) character(0) else res$pathway)
  if (length(miss)) {
    add <- data.table::data.table(pathway = miss, pval = NA_real_, padj = NA_real_, log2err = NA_real_,
                                  ES = NA_real_, NES = NA_real_, size = 0L,
                                  leadingEdge = rep(list(character(0)), length(miss)))
    res <- if (is.null(res)) add else rbind(res, add, fill = TRUE)
  }
  res[order(res$pathway)]
}

plotEnrichment <- function(pathway, stats, ...) {
  if (length(intersect(pathway, names(stats))) == 0)
    return(list(data = data.frame(rank = numeric(0), ES = numeric(0))))
  fgsea::plotEnrichment(pathway, stats, ...)
}

# metap: drop NA p-values (empty sets) before combining
two2one <- function(p, two = NULL, invert = NULL) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p) & !is.na(invert)
  if (any(ok)) out[ok] <- metap::two2one(p[ok], two = two[ok], invert = invert[ok])
  out
}
sumlog <- function(p, ...) {
  p <- p[!is.na(p)]
  if (!length(p)) return(list(chisq = NA_real_, df = NA_real_, p = NA_real_, validp = p))
  metap::sumlog(p, ...)
}

# GSVA ssGSEA: NA scores for empty gene sets
ssgseaParam <- function(exprData, geneSets, ...) list(expr = exprData, geneSets = geneSets)
gsva <- function(param, ...) {
  gs <- param$geneSets
  out <- matrix(NA_real_, nrow = length(gs), ncol = ncol(param$expr),
                dimnames = list(names(gs), colnames(param$expr)))
  ok <- names(gs)[lengths(gs) > 0]
  if (length(ok)) {
    r <- GSVA::gsva(GSVA::ssgseaParam(param$expr, gs[ok]), verbose = FALSE)
    out[rownames(r), ] <- r
  }
  out
}
t.test <- function(x, y, ...) {
  if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2)
    return(list(statistic = NA, parameter = NA, p.value = NA_real_))
  stats::t.test(x, y, ...)
}

# per-intervention diagnostic GSEA plots cannot be drawn for empty sets: skip them
print <- function(x, ...) {
  if (inherits(x, "ggplot")) invisible(try(base::print(x, ...), silent = TRUE)) else base::print(x, ...)
}

# 12_dementia_cohen.R calls data.table() before loading the package (original bug)
suppressPackageStartupMessages(library(data.table))
