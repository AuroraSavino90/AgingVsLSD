# ============================================================
# Common setup for the in vivo threshold-sensitivity copies of
# Scripts/05_DEvsAging_validation.R and Scripts/18_GO_allLSD.R.
# Sourced in place of setwd("workdir"). Working directory must be the
# repository root.
#
# Environment variables (set by run_all.R):
#   SELECTION_MODE  "pvalue" -> nominal p < 0.05 for the in vivo LSD doses (as in the paper code)
#                   "fdr"    -> BH-adjusted p < 0.05 for the in vivo LSD doses
#   SENS_OUT        output folder (default Threshold_sensitivity_invivo/Results_<mode>)
# Chronic LSD (GSE179379) is selected with FDR < 0.05 in both modes, as in the paper.
# ============================================================

selection_mode <- tolower(Sys.getenv("SELECTION_MODE", "pvalue"))
stopifnot(selection_mode %in% c("pvalue", "fdr"))

# DESeq2 column used for gene selection
P_SEQ <- if (selection_mode == "fdr") "padj" else "pvalue"

OUT <- Sys.getenv("SENS_OUT",
                  file.path("Threshold_sensitivity_invivo",
                            paste0("Results_", if (selection_mode == "fdr") "FDR" else "pvalue")))
dir.create(file.path(OUT, "RData"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "Figures"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists("Data") || !file.exists("Data/intermediate_RData/DE_GSE179379.RData"))
  stop("Run from the repository root (AgingVsLSD/): needs Data/intermediate_RData/DE_GSE179379.RData.")

# ------------------------------------------------------------
# Shims. With FDR some doses have no (or very few) selected genes.
# The original code assumes non-empty UP and DN sets; these wrappers
# return NA / empty results for an empty set instead of stopping.
# When the sets are non-empty they return exactly what the original
# functions return (and draw the same random numbers).
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

# metap: drop NA p-values (empty sets) before combining.
# Without NAs the original function is called with the original arguments
# (18 passes `two`/`invert` of length 5 with 3 p-values in the in vitro part).
two2one <- function(p, two = NULL, invert = NULL) {
  if (!anyNA(p) && !anyNA(invert)) return(metap::two2one(p, two = two, invert = invert))
  two <- rep_len(two, length(p)); invert <- rep_len(invert, length(p))
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p) & !is.na(invert)
  if (any(ok)) out[ok] <- metap::two2one(p[ok], two = two[ok], invert = invert[ok])
  out
}
sumlog <- function(p, ...) {
  if (!anyNA(p) && length(p)) return(metap::sumlog(p, ...))
  p <- p[!is.na(p)]
  if (!length(p)) return(list(chisq = NA_real_, df = NA_real_, p = NA_real_, validp = p))
  metap::sumlog(p, ...)
}

# enrichGO: an empty table instead of NULL when no gene can be tested
# (a NULL would drop the dose from the list of GO results)
enrichGO <- function(gene, ...) {
  r <- if (length(gene)) clusterProfiler::enrichGO(gene = gene, ...) else NULL
  if (is.null(r)) data.frame(Description = character(0), pvalue = numeric(0)) else r
}

# GSEA plots of an empty dose cannot be drawn: skip them
print <- function(x, ...) {
  if (inherits(x, "ggplot")) invisible(try(base::print(x, ...), silent = TRUE)) else base::print(x, ...)
}
