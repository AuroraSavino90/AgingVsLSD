# ============================================================
# THRESHOLD-SENSITIVITY COPY of Scripts/save_mean_logFC_selected_genes.R
# Changes vs the original: (1) gene selection threshold set by
# SELECTION_MODE for ALL datasets (P_ARRAY / P_SEQ columns);
# (2) outputs written to OUT (Threshold_sensitivity/Results_<mode>);
# (3) see sensitivity_setup.R. Run via Threshold_sensitivity/run_all.R
# ============================================================
# ============================================================
# Mean intervention logFC of the selected UP/DN genes (Table S8)
# and of the genes changing in the same / opposite direction
# relative to human AGING (Table S10) and DEMENTIA (Table S14)
# signatures. The two output files are the input of
# Cumulative_logFC.R.
# ============================================================

rm(list = ls())
source("Threshold_sensitivity/Scripts/sensitivity_setup.R")  # SELECTION_MODE, OUT, P_ARRAY/P_SEQ, shims

# Gene selection threshold, applied identically to every intervention dataset:
#   "pvalue" = nominal p < 0.05
#   "fdr"    = BH-adjusted p < 0.05 (limma adj.P.Val / DESeq2 padj)
# Set with Sys.setenv(SELECTION_MODE = "fdr") before sourcing.
selection_mode <- tolower(Sys.getenv("SELECTION_MODE", "pvalue"))
stopifnot(selection_mode %in% c("pvalue", "fdr"))
out_suffix <- if (selection_mode == "fdr") "_FDR" else "_pvalue"

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(limma)
  library(org.Rn.eg.db)
  library(AnnotationDbi)
})


# ------------------------------------------------------------
# 1. LOAD DATA IN SEPARATE ENVIRONMENTS
#    This avoids metadata objects overwriting one another.
# ------------------------------------------------------------

interv_env <- new.env(parent = globalenv())
aging_env  <- new.env(parent = globalenv())
dem_env    <- new.env(parent = globalenv())

load("Data/original_RData/Psychedelics_PFC.RData", envir = interv_env)
load("Data/original_RData/ExerciseAndAlcohol.RData", envir = interv_env)
load("Data/intermediate_RData/DE_GSE179379.RData", envir = interv_env)

load("Data/original_RData/Alldata_20Sep.RData", envir = aging_env)
load("Data/original_RData/Dementia_alldata_meta.RData", envir = dem_env)

rat_homologs   <- read.csv("Data/homologs/Human rat homologs.txt")  # Ensembl 105: BioMart export, human genes (GRCh38.p13) with rat homologs, March 2022
mouse_homologs <- read.csv("Data/homologs/Human mouse homologs.txt", sep = "\t")  # Ensembl 107: BioMart export, human genes (GRCh38.p13) with mouse homologs, August 2022

# ------------------------------------------------------------
# 2. INTERVENTION DATASETS: SAME OBJECTS / ORDER AS PROJECT
# ------------------------------------------------------------

DEGs_array <- c(
  "DE_MDMA", "DE_GSE26364", "DE_GSE23728", "DE_GSE179380",
  "DE_GSE64607", "DE_GSE72507_0", "DE_GSE60676_0", "DE_GSE28515_C57"
)

DEGs_seq <- c(
  "DE_DMT", "DE_GSE161626", "DE_GSE161626_48h", "DE_GSE161626_7d",
  "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks", "DE_GSE81672",
  "DE_harm", "DE_pharm", "DE_GSE164798_chronic", "DE_GSE111273_Y", "DE_GSE111273_O"
)

mouse_data <- c(
  "DE_GSE64607", "DE_GSE164798_chronic", "DE_GSE72507_0",
  "DE_GSE60676_0", "DE_GSE28515_C57", "DE_MDMA", "DE_GSE161626",
  "DE_GSE161626_48h", "DE_GSE161626_7d", "DE_GSE81672", "DE_GSE26364",
  "DE_GSE209859_FC3hours", "DE_GSE209859_FC4weeks",
  "DE_GSE111273_Y", "DE_GSE111273_O"
)

intervention_info <- data.table(
  DE = c(DEGs_array, DEGs_seq, "DE_GSE179379", "DE_GSE75772"),
  intervention = c(
    "MDMA", "Ketamine (long term)", "DOI (2h)", "LSD (single)",
    "Exercise (28d)", "Alcohol (vapour chamber)", "Alcohol (vapour chamber 2)", "Alcohol (IP injection)",
    "DMT", "DOI (24h)", "DOI (48h)", "DOI (7days)",
    "Psilocybin (3h)", "Psilocybin (4weeks)", "Ketamine (single)",
    "Harmaline", "Pharmahuasca", "Exercise (30d)", "EE (Young)", "EE (Old)",
    "LSD (chronic)", "Rat aging"
  ),
  type = c(
    "Psychoplastogen", "Psychoplastogen", "Psychoplastogen", "Psychoplastogen",
    "Positive Control",
    "Negative Control", "Negative Control", "Negative Control",
    "Psychoplastogen", "Psychoplastogen", "Psychoplastogen", "Psychoplastogen",
    "Psychoplastogen", "Psychoplastogen", "Psychoplastogen", "Psychoplastogen",
    "Psychoplastogen",
    "Positive Control", "Positive Control", "Positive Control",
    "Psychoplastogen", "Negative Control"
  )
)

intervention_info[, species := fifelse(DE %in% mouse_data, "mouse", "rat")]
intervention_info[, data_class := fifelse(DE %in% c(DEGs_array, "DE_GSE75772"), "array", "RNA-seq")]

# ------------------------------------------------------------
# 3. RAT AGING CONTROL (GSE75772): Aged vs Young
# ------------------------------------------------------------

GSE75772 <- read.csv("Data/raw_data_notclosed/GSE75772_DEseq-normalized.counts_mPFC.csv.gz",
                     row.names = 1, check.names = FALSE)

gene_symbols <- as.character(AnnotationDbi::mapIds(
  org.Rn.eg.db, keys = gsub("\\..*$", "", rownames(GSE75772)),
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"
))
keep <- !is.na(gene_symbols) & gene_symbols != ""
GSE75772 <- GSE75772[keep, , drop = FALSE]
gene_symbols <- gene_symbols[keep]

# duplicated symbols: keep the most highly expressed row
keep_rows <- unlist(lapply(split(seq_along(gene_symbols), gene_symbols), function(ii) {
  if (length(ii) == 1) return(ii)
  ii[which.max(rowMeans(GSE75772[ii, , drop = FALSE], na.rm = TRUE))]
}), use.names = FALSE)
GSE75772 <- GSE75772[keep_rows, , drop = FALSE]
rownames(GSE75772) <- gene_symbols[keep_rows]
GSE75772 <- log2(GSE75772 + 1)

age_lab <- unlist(read.csv("Data/raw_data_notclosed/GSE75772_series_matrix.txt", sep = "\t")[10, 2:32])
aged  <- data.frame(GSE75772[, age_lab == "age: Aged",  drop = FALSE])
young <- data.frame(GSE75772[, age_lab == "age: Young", drop = FALSE])

gs <- factor(c(rep("Cond1", ncol(aged)), rep("Cond2", ncol(young))))
design <- model.matrix(~ gs + 0)
fit <- lmFit(cbind.data.frame(aged, young), design)
fit <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = "gsCond1-gsCond2", levels = design)), 0.01)
assign("DE_GSE75772", topTable(fit, adjust = "fdr", sort.by = "B", number = nrow(fit)), envir = interv_env)

# ------------------------------------------------------------
# 4. SELECTED INTERVENTION GENES AND THEIR LOGFC
# ------------------------------------------------------------

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

extract_selected_genes <- function(DE_name, species, data_class) {
  x <- get(DE_name, envir = interv_env, inherits = FALSE)

  if (data_class == "array") {
    logfc_col <- "logFC"
    p_col <- if (selection_mode == "fdr") "adj.P.Val" else "P.Value"
  } else {
    logfc_col <- "log2FoldChange"
    p_col <- if (selection_mode == "fdr") "padj" else "pvalue"
  }

  dt <- data.table(
    original_gene = rownames(x),
    original_logFC = as.numeric(x[[logfc_col]]),
    selection_p = as.numeric(x[[p_col]])
  )
  dt <- dt[!is.na(original_logFC) & !is.na(selection_p) &
           selection_p < 0.05 & original_logFC != 0]
  dt[, direction := fifelse(original_logFC > 0, "UP", "DN")]

  # map to human genes; if several rodent orthologues map to the same
  # human gene in the same direction, use their mean logFC
  h <- if (species == "mouse") mouse_homologs else rat_homologs
  map <- unique(data.table(original_gene = as.character(h[, 2]), human_gene = as.character(h[, 3])))
  map <- map[!is.na(original_gene) & original_gene != "" & !is.na(human_gene) & human_gene != ""]

  merge(dt, map, by = "original_gene", allow.cartesian = TRUE)[, .(
    intervention_logFC = mean(original_logFC, na.rm = TRUE)
  ), by = .(human_gene, direction)]
}

selected_genes <- rbindlist(lapply(seq_len(nrow(intervention_info)), function(i) {
  info <- intervention_info[i]
  z <- extract_selected_genes(info$DE, info$species, info$data_class)
  z[, `:=`(DE = info$DE, intervention = info$intervention, type = info$type)]
  z
}))

# mean logFC of all selected genes (Table S8)
mean_logFC_selected <- selected_genes[, .(
  mean_logFC_selected = safe_mean(intervention_logFC)
), by = .(DE, intervention, type, direction)]

# ------------------------------------------------------------
# 5. HUMAN AGING SIGNATURES
#    Gene-wise correlation of expression with age (subjects >= 20 years).
# ------------------------------------------------------------

build_aging_signatures <- function(env) {
  metadata <- get("metadata", envir = env, inherits = FALSE)
  out <- list()

  for (region in c("DLPFC", "PFC")) {
    dats <- na.omit(unique(metadata$Dataset[
      metadata$Organism == "Homo sapiens" & metadata$Diagnosis == "Healthy" &
      metadata$Region_simpl %in% region
    ]))
    dats <- setdiff(dats, c("GSE102741", "GSE5388"))

    for (dd in dats) {
      if (!exists(dd, envir = env, inherits = FALSE)) next
      expr <- get(dd, envir = env, inherits = FALSE)
      md <- metadata[which(metadata$Dataset == dd & metadata$Organism == "Homo sapiens" &
                           metadata$Diagnosis == "Healthy" & metadata$Region_simpl %in% region), , drop = FALSE]

      sam <- intersect(md$Sample, colnames(expr))
      if (length(sam) == 0) next
      md <- md[match(sam, md$Sample), , drop = FALSE]

      age <- suppressWarnings(as.numeric(md$Age))
      keep <- !is.na(age) & age >= 7300
      if (sum(keep) < 10) next

      expr2 <- expr[, sam[keep], drop = FALSE]
      age_cor <- as.numeric(suppressWarnings(cor(t(expr2), age[keep], use = "pairwise.complete.obs")))

      out[[length(out) + 1]] <- data.table(
        reference_dataset = dd,
        reference_region = region,
        reference_condition = "Healthy aging",
        human_gene = rownames(expr2),
        reference_effect = age_cor
      )
    }
  }

  ans <- rbindlist(out)
  ans[!is.na(reference_effect) & reference_effect != 0]
}

# ------------------------------------------------------------
# 6. HUMAN DEMENTIA SIGNATURES
#    Gene-wise effect = disease - Healthy (limma).
# ------------------------------------------------------------

dementia_contrasts <- data.table(
  reference_region = c("DLPFC", "PFC", "PFC", "DLPFC", "DLPFC"),
  ctrl = "Healthy",
  case = c("Alzheimer", "Alzheimer", "Huntington", "PDD", "DLB")
)

build_dementia_signatures <- function(env) {
  metadata <- get("metadata", envir = env, inherits = FALSE)
  out <- list()

  for (cc in seq_len(nrow(dementia_contrasts))) {
    reg  <- dementia_contrasts$reference_region[cc]
    ctrl <- dementia_contrasts$ctrl[cc]
    case <- dementia_contrasts$case[cc]

    dats <- na.omit(unique(metadata$Dataset[
      metadata$Diagnosis %in% c(ctrl, case) & metadata$Region %in% reg
    ]))

    for (dd in dats) {
      if (!exists(dd, envir = env, inherits = FALSE)) next
      expr <- get(dd, envir = env, inherits = FALSE)
      # discard genes with missing expression values, as in the dementia scripts
      expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]

      md <- metadata[which(metadata$Dataset == dd & metadata$Diagnosis %in% c(ctrl, case) &
                           metadata$Region %in% reg), , drop = FALSE]
      sam <- intersect(md$Sample, colnames(expr))
      if (length(sam) == 0) next
      md <- md[match(sam, md$Sample), , drop = FALSE]

      grp <- factor(md$Diagnosis, levels = c(ctrl, case))
      ok <- !is.na(grp)
      sam <- sam[ok]
      grp <- droplevels(grp[ok])
      if (sum(grp == ctrl) < 5 || sum(grp == case) < 5) next

      design <- model.matrix(~ 0 + grp)
      colnames(design) <- c(ctrl, case)
      fit <- lmFit(expr[, sam, drop = FALSE], design)
      fit <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = paste0(case, "-", ctrl), levels = design)))
      tt <- topTable(fit, number = Inf, sort.by = "none")

      out[[length(out) + 1]] <- data.table(
        reference_dataset = dd,
        reference_region = reg,
        reference_condition = paste0(case, " vs ", ctrl),
        human_gene = rownames(tt),
        reference_effect = as.numeric(tt$logFC)
      )
    }
  }

  ans <- rbindlist(out)
  ans[!is.na(reference_effect) & reference_effect != 0]
}

# ------------------------------------------------------------
# 7. COMPARE INTERVENTION DIRECTION WITH REFERENCE SIGNATURE
#    (Tables S10 and S14)
# ------------------------------------------------------------

compare_with_reference <- function(intervention_dt, reference_dt) {
  x <- merge(intervention_dt, reference_dt, by = "human_gene", allow.cartesian = TRUE)
  x[, relation := fifelse(sign(intervention_logFC) == sign(reference_effect),
                          "same_direction", "opposite_direction")]

  ans <- x[, .(
    n_genes_tested = uniqueN(human_gene),
    n_opposite_direction = uniqueN(human_gene[relation == "opposite_direction"]),
    mean_logFC_opposite_direction = safe_mean(intervention_logFC[relation == "opposite_direction"]),
    n_same_direction = uniqueN(human_gene[relation == "same_direction"]),
    mean_logFC_same_direction = safe_mean(intervention_logFC[relation == "same_direction"])
  ), by = .(DE, intervention, type, direction,
            reference_dataset, reference_region, reference_condition)]
  setorder(ans, DE, intervention, type, direction)
  ans
}

aging_summary <- compare_with_reference(selected_genes, build_aging_signatures(aging_env))

# rat aging is an aging-only control
dementia_summary <- compare_with_reference(selected_genes[DE != "DE_GSE75772"],
                                           build_dementia_signatures(dem_env))

# ------------------------------------------------------------
# 8. SAVE
# ------------------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb, "Mean_logFC_selected")
writeData(wb, "Mean_logFC_selected", mean_logFC_selected)
addWorksheet(wb, "Detailed_summary")
writeData(wb, "Detailed_summary", aging_summary)
saveWorkbook(wb, paste0(OUT, "/logFC_vs_aging_signatures", out_suffix, ".xlsx"), overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "Detailed_summary")
writeData(wb, "Detailed_summary", dementia_summary)
saveWorkbook(wb, paste0(OUT, "/logFC_vs_dementia_signatures", out_suffix, ".xlsx"), overwrite = TRUE)
