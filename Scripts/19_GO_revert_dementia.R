rm(list=ls())
setwd("workdir") 


dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/RData", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# LOAD
# -----------------------------
load("Data/Psychedelics_PFC.RData")
load("Data/Alldata_20Sep.RData")
load("Data/ExerciseAndAlcohol.RData")
load("Data/Dementia_alldata_meta.RData")
load("Results/RData/DE_GSE179379.RData")
homologs <- read.csv("Data/Human rat homologs.txt")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(metap)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------
# Helpers
# -----------------------------
quietly <- function(expr) {
  withCallingHandlers(
    suppressWarnings(suppressMessages(expr)),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

# t-stat gene-wise per (case - ctrl), expr genes x samples
fast_tstat <- function(expr, grp, ctrl_level, case_level) {
  grp <- factor(grp, levels = c(ctrl_level, case_level))
  if (any(is.na(grp))) stop("grp contains NA after factoring.")
  if (sum(grp == ctrl_level) < 3 || sum(grp == case_level) < 3) return(rep(NA_real_, nrow(expr)))
  
  x0 <- expr[, grp == ctrl_level, drop = FALSE]
  x1 <- expr[, grp == case_level, drop = FALSE]
  
  n0 <- ncol(x0); n1 <- ncol(x1)
  m0 <- rowMeans(x0, na.rm = TRUE)
  m1 <- rowMeans(x1, na.rm = TRUE)
  v0 <- apply(x0, 1, var, na.rm = TRUE)
  v1 <- apply(x1, 1, var, na.rm = TRUE)
  
  se <- sqrt(v0 / n0 + v1 / n1)
  t  <- (m1 - m0) / se
  t[!is.finite(t)] <- NA_real_
  t
}

make_div_breaks <- function(mat, paletteLength = 50) {
  myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
  rng <- range(as.numeric(mat), na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) rng <- c(-1, 1)
  myBreaks <- c(
    seq(rng[1], 0, length.out = ceiling(paletteLength/2) + 1),
    seq(rng[2] / paletteLength, rng[2], length.out = floor(paletteLength/2))
  )
  list(color=myColor, breaks=myBreaks)
}

meta_dir_sumlog <- function(p_vec, nes_vec, expect_negative = TRUE) {
  istwo <- rep(TRUE, length(p_vec))
  if (expect_negative) {
    toinvert <- ifelse(sign(nes_vec) == 1, TRUE, FALSE)   # want NES < 0
  } else {
    toinvert <- ifelse(sign(nes_vec) == -1, TRUE, FALSE)  # want NES > 0
  }
  sumlog(two2one(p_vec, two = istwo, invert = toinvert))
}

load("Results/ego_up_valid.RData")     
load("Results/ego_down_valid.RData")   

DEG_down<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange<0 & DE_GSE179379$padj<0.05)]
DEG_up<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

library(clusterProfiler)
library(org.Rn.eg.db)   
ego_down[[5]]<- enrichGO(gene = DEG_down,
                         keyType="SYMBOL",
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")

ego_up[[5]]<- enrichGO(gene = DEG_up,
                       keyType="SYMBOL",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       universe=rownames(DE_GSE179379), OrgDb="org.Rn.eg.db")


build_lsd_meta_p <- function(ego_list) {
  allpaths <- c()
  for (i in seq_along(ego_list)) allpaths <- union(allpaths, ego_list[[i]]$Description)
  
  allp_mat <- matrix(NA_real_, nrow = length(allpaths), ncol = length(ego_list))
  rownames(allp_mat) <- allpaths
  for (i in seq_along(ego_list)) {
    allp_mat[ego_list[[i]]$Description, i] <- ego_list[[i]]$pvalue
  }
  
  pmeta <- vector("list", length(allpaths))
  names(pmeta) <- allpaths
  
  for (k in seq_along(allpaths)) {
    path <- allpaths[k]
    missing <- which(is.na(allp_mat[path, ]))
    istwo <- rep(FALSE, ncol(allp_mat))
    toinvert <- rep(FALSE, ncol(allp_mat)) 
    
    if (length(missing) == 0) {
      pmeta[[k]] <- sumlog(two2one(allp_mat[path, ], two = istwo, invert = toinvert))[[3]]
    } else if (length(missing) < ncol(allp_mat)) {
      pmeta[[k]] <- sumlog(two2one(allp_mat[path, ][-missing], two = istwo[-missing], invert = toinvert[-missing]))[[3]]
    } else {
      pmeta[[k]] <- NA_real_
    }
  }
  
  pmeta <- unlist(pmeta)
  pmeta <- p.adjust(pmeta, method = "fdr")
  list(allp_mat = allp_mat, pmeta = pmeta)
}

lsd_up <- build_lsd_meta_p(ego_up)
lsd_dn <- build_lsd_meta_p(ego_down)

LSD_up_terms <- names(lsd_up$pmeta)[which(lsd_up$pmeta < 0.05)]
LSD_dn_terms <- names(lsd_dn$pmeta)[which(lsd_dn$pmeta < 0.05)]


req_cols <- c("Dataset","Sample","Diagnosis","Region")
miss_cols <- setdiff(req_cols, colnames(metadata))
if (length(miss_cols) > 0) stop(paste("metadata missing columns:", paste(miss_cols, collapse=", ")))

get_dat <- function(region, diagnosis_vec) {
  na.omit(unique(metadata$Dataset[metadata$Diagnosis %in% diagnosis_vec & metadata$Region %in% region]))
}

contrast_list <- list(
  list(region="DLPFC", ctrl="Healthy", case="Alzheimer",  datasets=get_dat("DLPFC", c("Healthy","Alzheimer"))),
  list(region="PFC",   ctrl="Healthy", case="Alzheimer",  datasets=get_dat("PFC",   c("Healthy","Alzheimer"))),
  list(region="PFC",   ctrl="Healthy", case="Huntington", datasets=get_dat("PFC",   c("Healthy","Huntington"))),
  list(region="DLPFC", ctrl="Healthy", case="PDD",        datasets=get_dat("DLPFC", c("Healthy","PDD"))),
  list(region="DLPFC", ctrl="Healthy", case="DLB",        datasets=get_dat("DLPFC", c("Healthy","DLB")))
)

gsea_dem_list <- list()

for (cc in contrast_list) {
  reg  <- cc$region
  ctrl <- cc$ctrl
  case <- cc$case
  dats <- cc$datasets
  if (length(dats) == 0) next
  
  for (dd in dats) {
    if (!exists(dd, envir = .GlobalEnv)) next
    expr <- get(dd)
    
    if (sum(is.na(expr)) > 0) {
      expr <- expr[-which(is.na(expr), arr.ind = TRUE)[,1], , drop = FALSE]
    }
    
    idx <- which(metadata$Dataset == dd & metadata$Diagnosis %in% c(ctrl, case) & metadata$Region %in% reg)
    if (length(idx) < 10) next
    
    sam <- metadata$Sample[idx]
    dis <- metadata$Diagnosis[idx]
    
    sam <- intersect(sam, colnames(expr))
    if (length(sam) < 10) next
    
    expr2 <- expr[, sam, drop = FALSE]
    
    
    dis2 <- metadata$Diagnosis[match(sam, metadata$Sample)]
    if (any(is.na(dis2))) next
    
    if (sum(dis2 == ctrl) < 5 || sum(dis2 == case) < 5) next
    
    tstat <- fast_tstat(expr2, dis2, ctrl_level = ctrl, case_level = case)
    names(tstat) <- rownames(expr2)
    tstat <- tstat[!is.na(tstat)]
    if (length(tstat) < 200) next
    tstat <- sort(tstat, decreasing = TRUE)
    
    gsea_dem <- quietly(gseGO(
      geneList = tstat,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = FALSE,
      by = "fgsea"
    ))
    
    name_ds <- paste(dd, reg, case, sep="__")
    gsea_dem_list[[name_ds]] <- gsea_dem
  }
}

save(gsea_dem_list, file = "Results/RData/gsea_GO_dementia_exactContrasts.RData")

all_terms_dem <- c()
for (nm in names(gsea_dem_list)) {
  all_terms_dem <- union(all_terms_dem, gsea_dem_list[[nm]]$Description)
}

p_mat_dem <- matrix(NA_real_, nrow=length(all_terms_dem), ncol=length(gsea_dem_list),
                    dimnames=list(all_terms_dem, names(gsea_dem_list)))
nes_mat_dem <- p_mat_dem

for (nm in names(gsea_dem_list)) {
  p_mat_dem[gsea_dem_list[[nm]]$Description, nm]   <- gsea_dem_list[[nm]]$pvalue
  nes_mat_dem[gsea_dem_list[[nm]]$Description, nm] <- gsea_dem_list[[nm]]$NES
}

dir_meta_for_terms <- function(terms, expect_negative = TRUE) {
  out_p <- setNames(rep(NA_real_, length(terms)), terms)
  
  for (tt in terms) {
    if (!tt %in% rownames(p_mat_dem)) next
    pv <- p_mat_dem[tt, ]
    nv <- nes_mat_dem[tt, ]
    ok <- which(!is.na(pv) & !is.na(nv))
    if (length(ok) < 2) next
    
    res <- meta_dir_sumlog(pv[ok], nv[ok], expect_negative = expect_negative)
    out_p[tt] <- res$p
  }
  
  
  keep <- which(!is.na(out_p))
  out_q <- out_p
  out_q[keep] <- p.adjust(out_p[keep], method="fdr")
  out_q
}

LSD_up_terms_both<-intersect(LSD_up_terms, rownames(p_mat_dem))
q_dem_for_LSDup <- dir_meta_for_terms(LSD_up_terms_both, expect_negative = TRUE)

LSD_dn_terms_both<-intersect(LSD_dn_terms, rownames(p_mat_dem))
q_dem_for_LSDdn <- dir_meta_for_terms(LSD_dn_terms_both, expect_negative = FALSE)

library(openxlsx)
colnames(nes_mat_dem)<-paste(colnames(nes_mat_dem), "NES")
colnames(p_mat_dem)<-paste(colnames(p_mat_dem), "p-value")
all_stats_up<-cbind(pathway=LSD_up_terms_both, nes_mat_dem[LSD_up_terms_both, ], p_mat_dem[LSD_up_terms_both, ], data.frame(lsd_up)[LSD_up_terms_both,-6], 
                    logFDRLSD=-log10(unlist(data.frame(lsd_up)[LSD_up_terms_both,6])),logFDRAging=-log10(q_dem_for_LSDup[LSD_up_terms_both]))
write.xlsx(all_stats_up, "results/GO_all_stats_up_dementia.xlsx")

all_stats_dn<-cbind(pathway=LSD_dn_terms_both, nes_mat_dem[LSD_dn_terms_both, ], p_mat_dem[LSD_dn_terms_both, ], data.frame(lsd_dn)[LSD_dn_terms_both,-6], 
                    logFDRLSD=-log10(unlist(data.frame(lsd_dn)[LSD_dn_terms_both,6])),logFDRAging=-log10(q_dem_for_LSDup[LSD_dn_terms_both]))
write.xlsx(all_stats_dn, "results/GO_all_stats_dn_dementia.xlsx")

# ============================================================
# IN VITRO validation 
# ============================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(openxlsx)
  library(sva)
  library(GSVA)
  library(limma)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

thr_invitro <- 0.05
outdir_inv <- "Results/OppositeGO_invitro_tests_DEMENTIA"   
dir.create(outdir_inv, showWarnings = FALSE, recursive = TRUE)

xlsx_metadata  <- "Data/In vitro RNA-seq_all/Metadata.xlsx"
tsv_counts     <- "Data/In vitro RNA-seq_all/salmon.merged.gene_counts.tsv"

# -----------------------------
# Helper functions 
# -----------------------------
changenames <- function(data, anno){
  annotation_sel <- anno[match(rownames(data), anno[,1]), 2]
  
  if(length(which(annotation_sel == "")) > 0){
    data <- data[-which(annotation_sel == ""), ]
    annotation_sel <- annotation_sel[-which(annotation_sel == "")]
  }
  
  a <- which(duplicated(annotation_sel))
  while(length(a) > 0){
    for(i in seq_along(unique(annotation_sel))){
      if(length(which(annotation_sel == unique(annotation_sel)[i])) > 1){
        m <- which.max(rowMeans(data[which(annotation_sel == unique(annotation_sel)[i]), ], na.rm = TRUE))
        data <- data[-which(annotation_sel == unique(annotation_sel)[i])[-m], ]
        annotation_sel <- annotation_sel[-which(annotation_sel == unique(annotation_sel)[i])[-m]]
      }
    }
    data <- data[which(is.na(annotation_sel) == FALSE), ]
    annotation_sel <- na.omit(annotation_sel)
    a <- which(duplicated(annotation_sel))
  }
  
  rownames(data) <- annotation_sel
  data
}

clamp01 <- function(x, min_p=1e-300){
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x[x <= 0] <- min_p
  x[x > 1]  <- 1
  x
}

build_desc2go_from_ego <- function(ego_list){
  tabs <- lapply(ego_list, function(x){
    if(is.null(x)) return(NULL)
    df <- as.data.frame(x)
    if(!all(c("ID","Description") %in% colnames(df))) return(NULL)
    df[, c("ID","Description"), drop=FALSE]
  })
  tabs <- tabs[!vapply(tabs, is.null, logical(1))]
  if(length(tabs) == 0) return(data.frame(ID=character(), Description=character()))
  bind_rows(tabs) %>% distinct(Description, .keep_all = TRUE)
}

goid_to_genesets_mouse <- function(go_ids){
  go_ids <- unique(na.omit(go_ids))
  if(length(go_ids) == 0) return(list())
  
  sel <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys     = go_ids,
    keytype  = "GOALL",
    columns  = c("SYMBOL")
  )
  sel <- sel %>%
    filter(!is.na(SYMBOL), SYMBOL != "", !is.na(GOALL)) %>%
    distinct(GOALL, SYMBOL)
  
  gs <- split(sel$SYMBOL, sel$GOALL)
  gs <- lapply(gs, unique)
  gs
}

genesets_to_index <- function(genesets, rownames_expr, min_size = 10){
  idx <- lapply(genesets, function(genes) which(rownames_expr %in% genes))
  idx <- idx[lengths(idx) >= min_size]
  idx
}

# -----------------------------
# Build the "opposite GO" sets from LSD vs Dementia
# ----------------------------

opp_UP_LSD_DOWN_DEM <- names(q_dem_for_LSDup)[which(!is.na(q_dem_for_LSDup) & q_dem_for_LSDup < thr_invitro)]
opp_DN_LSD_UP_DEM   <- names(q_dem_for_LSDdn)[which(!is.na(q_dem_for_LSDdn) & q_dem_for_LSDdn < thr_invitro)]

df_opp_dem <- bind_rows(
  tibble(OppositePattern = "UP in LSD / DOWN in Dementia", Description = opp_UP_LSD_DOWN_DEM,
         FDR_DementiaOpp = q_dem_for_LSDup[opp_UP_LSD_DOWN_DEM]),
  tibble(OppositePattern = "DOWN in LSD / UP in Dementia", Description = opp_DN_LSD_UP_DEM,
         FDR_DementiaOpp = q_dem_for_LSDdn[opp_DN_LSD_UP_DEM])
) %>%
  mutate(FDR_DementiaOpp = clamp01(FDR_DementiaOpp)) %>%
  arrange(OppositePattern, FDR_DementiaOpp)

write.csv(df_opp_dem, file=file.path(outdir_inv, "GO_opposite_LSD_vs_Dementia.csv"), row.names=FALSE)
saveRDS(df_opp_dem,  file=file.path(outdir_inv, "GO_opposite_LSD_vs_Dementia.rds"))


desc2go <- bind_rows(
  build_desc2go_from_ego(ego_up),
  build_desc2go_from_ego(ego_down)
) %>% distinct(Description, .keep_all = TRUE)

df_opp_dem <- df_opp_dem %>%
  left_join(desc2go, by="Description") %>%
  rename(GO_ID = ID)

write.csv(df_opp_dem, file=file.path(outdir_inv, "GO_opposite_LSD_vs_Dementia_with_GOIDs.csv"), row.names=FALSE)

go_ids <- unique(na.omit(df_opp_dem$GO_ID))
genesets_goid <- goid_to_genesets_mouse(go_ids)

desc_by_go <- df_opp_dem %>% filter(!is.na(GO_ID)) %>% distinct(GO_ID, Description) %>% deframe()

genesets_desc <- genesets_goid[names(genesets_goid) %in% names(desc_by_go)]
names(genesets_desc) <- desc_by_go[names(genesets_desc)]
genesets_desc <- genesets_desc[lengths(genesets_desc) >= 10]

saveRDS(genesets_desc, file=file.path(outdir_inv, "genesets_opposite_GO_desc_mouse.rds"))

# -----------------------------
# Load in vitro RNA-seq + ComBat-seq + log2 RPM (same as invitro_GO)
# -----------------------------
if(!file.exists(xlsx_metadata)) stop("Missing: ", xlsx_metadata)
if(!file.exists(tsv_counts))    stop("Missing: ", tsv_counts)

meta_iv <- openxlsx::read.xlsx(xlsx_metadata, rowNames = TRUE)

counts <- read.table(tsv_counts, header = TRUE, sep = "\t", check.names = FALSE)
rownames(counts) <- counts$gene_id
anno <- counts[, c(1:2)]
counts <- changenames(counts[, -c(1:2)], anno = anno)
counts <- round(counts)

rownames(meta_iv) <- colnames(counts)

meta_iv$Treatment  <- as.factor(meta_iv$Treatment)
meta_iv$Replicate  <- as.factor(meta_iv$Replicate)

meta_iv$Trt <- rep(c("Ct", "LSD", "Abeta", "Abeta", "Abeta",
                     "Abeta + LSD", "Abeta + LSD", "Abeta + LSD"), 3)
meta_iv$Trt <- factor(meta_iv$Trt, levels = c("Ct", "LSD", "Abeta", "Abeta + LSD"))

adj <- sva::ComBat_seq(as.matrix(counts), batch = meta_iv$Replicate, group = NULL)
expr <- log2(t(t(adj) / colSums(adj)) * 1e6 + 1)

grp <- droplevels(meta_iv$Trt)
keep <- grp %in% c("Abeta", "Abeta + LSD")
grp2 <- droplevels(grp[keep])
expr_sub <- expr[, keep, drop=FALSE]


# -----------------------------
# CAMERA (limma) on gene-level expression
# -----------------------------
design_s <- model.matrix(~ 0 + grp2)
colnames(design_s) <- make.names(levels(grp2))

contr <- limma::makeContrasts(Abeta...LSD - Abeta, levels = design_s)


idx <- genesets_to_index(genesets_desc, rownames(expr_sub), min_size = 10)

design_g <- model.matrix(~ 0 + grp2)
colnames(design_g) <- levels(grp2)

cam <- limma::camera(expr_sub, index = idx, design = design_g, contrast = contr[,1])
tab_camera <- cam %>%
  as.data.frame() %>%
  rownames_to_column("Description") %>%
  as_tibble() %>%
  arrange(PValue)

write.csv(tab_camera, file=file.path(outdir_inv, "camera_AbetaLSD_vs_Abeta.csv"), row.names=FALSE)

# -----------------------------
# Merge + save
# -----------------------------
merged_dem <- 
  tab_camera %>% transmute(Description, camera_Dir=Direction, camera_P=PValue, camera_FDR=FDR, camera_NGenes=NGenes) %>%
  left_join(df_opp_dem %>% select(OppositePattern, Description, GO_ID, FDR_DementiaOpp), by="Description") %>%
  arrange(camera_FDR)

write.csv(merged_dem, file=file.path(outdir_inv, "Merged_camera_oppositeGO_LSD_vs_Dementia.csv"), row.names=FALSE)
saveRDS(merged_dem,  file=file.path(outdir_inv, "Merged_camera_oppositeGO_LSD_vs_Dementia.rds"))



# -----------------------------
# Plotting: heatmap signed + top20 + scatter
# -----------------------------
plot_block <- function(term_q, prefix, lsd_pmeta_vec) {
  
  sel <- names(term_q)[which(!is.na(term_q) & term_q < 0.05)]
  if (length(sel) == 0) {
    message("No GO terms pass dementia directional FDR<0.05 for ", prefix)
    return(invisible(NULL))
  }
  
  # heatmap matrix: signed -log10(p) * sign(NES) per dataset demenza
  mat <- -log10(p_mat_dem[sel, , drop=FALSE]) * sign(nes_mat_dem[sel, , drop=FALSE])
  
  # ordina per significatività meta demenza
  ord <- order(-log10(term_q[sel]), decreasing = TRUE)
  mat <- mat[sel[ord], , drop=FALSE]
  
  pal <- make_div_breaks(mat, paletteLength = 50)
  
  anno <- data.frame(Dementia_meta = -log10(term_q[rownames(mat)]))
  rownames(anno) <- rownames(mat)
  
  pdf(file.path("Results/Figures", paste0(prefix, "_GO_shared.pdf")),
      width = 16, height = max(8, 0.25*nrow(mat)))
  pheatmap(mat,
           cluster_cols = FALSE, cluster_rows = FALSE,
           cellwidth = 12, cellheight = 12,
           breaks = pal$breaks, color = pal$color,
           annotation_row = anno)
  dev.off()
  
  # top20 by combined evidence: (-log10 LSD meta) + (-log10 dementia meta)
  common <- intersect(rownames(mat), names(lsd_pmeta_vec))
  comb <- (-log10(lsd_pmeta_vec[common])) + (-log10(term_q[common]))
  top20 <- names(sort(comb, decreasing = TRUE))[1:min(20, length(comb))]
  
  mat2 <- -log10(p_mat_dem[top20, , drop=FALSE]) * sign(nes_mat_dem[top20, , drop=FALSE])
  mat2 <- mat2[top20[order(comb[top20], decreasing=TRUE)], , drop=FALSE]
  pal2 <- make_div_breaks(mat2, paletteLength = 50)
  
  anno2 <- data.frame(Dementia_meta = -log10(term_q[rownames(mat2)]))
  rownames(anno2) <- rownames(mat2)
  
  pdf(file.path("Results/Figures", paste0(prefix, "_GO_shared_top20.pdf")), 12, 8)
  pheatmap(mat2,
           cluster_cols = FALSE, cluster_rows = FALSE,
           cellwidth = 12, cellheight = 12,
           breaks = pal2$breaks, color = pal2$color,
           annotation_row = anno2)
  dev.off()
  
  # scatter: x = -log10(FDR_LSD_meta), y = -log10(FDR_dementia_directional_meta)
  df_sc <- data.frame(
    x = -log10(lsd_pmeta_vec[rownames(mat)]),
    y = -log10(term_q[rownames(mat)]),
    term = rownames(mat)
  )
  
  pdf(file.path("Results/Figures", paste0(prefix, "_LSDvsDementia_scatter.pdf")), 7, 7)
  print(
    ggplot(df_sc, aes(x = x, y = y, label = term)) +
      geom_point() +
      geom_label_repel(label.size = 0.15, max.overlaps = 30) +
      theme_classic() +
      labs(x = "-log10(FDR LSD meta)", y = "-log10(FDR Dementia directional meta)")
  )
  dev.off()
  
  invisible(list(shared=rownames(mat), top20=rownames(mat2)))
}

# GO-up LSD (da ego_up meta): confronta con demenza "opposta"
plot_block(q_dem_for_LSDup, "GOup_LSD_vs_DEMENTIA", lsd_up$pmeta)

# GO-down LSD (da ego_down meta): confronta con demenza "opposta"
plot_block(q_dem_for_LSDdn, "GOdn_LSD_vs_DEMENTIA", lsd_dn$pmeta)

message("DONE.")
message(" - Dementia GSEA saved: Results/RData/gsea_GO_dementia_exactContrasts.RData")
message(" - Figures in: Results/Figures/")


inv_camera_csv <- file.path(outdir_inv, "camera_AbetaLSD_vs_Abeta.csv")
inv_merged_csv <- file.path(outdir_inv, "Merged_ssGSEA_camera_oppositeGO.csv")

  inv_cam <- read.csv(inv_camera_csv, stringsAsFactors = FALSE, check.names = FALSE) %>%
    as_tibble() %>%
    { if(!"Description" %in% colnames(.)) tibble::rownames_to_column(as.data.frame(.), "Description") else . } %>%
    transmute(
      Description,
      invitro_camera_P   = PValue,
      invitro_camera_FDR = FDR,
      invitro_camera_Dir = Direction,
      invitro_camera_N   = NGenes
    )

# -----------------------------
# helpers
# -----------------------------
min_p <- 1e-300

clamp01 <- function(x){
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x[x <= 0] <- min_p
  x[x > 1]  <- 1
  x
}

fdr_bin3 <- function(fdr){
  fdr <- clamp01(fdr)
  b <- dplyr::case_when(
    fdr < 0.05 ~ "FDR < 0.05",
    fdr < 0.10 ~ "FDR < 0.10",
    fdr < 0.15 ~ "FDR < 0.15",
    TRUE       ~ "FDR ≥ 0.15"
  )
  factor(b, levels = c("FDR ≥ 0.15","FDR < 0.15","FDR < 0.10","FDR < 0.05"))
}


# -----------------------------
# Build direction-specific data.frames
# Direction 1: LSD-UP (ego_up meta) vs Dementia expected DOWN
#   - dementia directional meta is q_dem_for_LSDup  (FDR-adjusted)
#   - LSD meta FDR is lsd_up$pmeta (already FDR-adjusted)
# -----------------------------
df_LSDup_DEMdown <- tibble(
  Description  = names(q_dem_for_LSDup),
  FDR_Dementia = as.numeric(q_dem_for_LSDup),
  FDR_LSD      = as.numeric(lsd_up$pmeta[names(q_dem_for_LSDup)])
) %>%
  left_join(inv_cam, by = "Description") %>%
  filter(!is.na(FDR_LSD), !is.na(FDR_Dementia))

# Direction 2: LSD-DOWN (ego_down meta) vs Dementia expected UP
df_LSDdn_DEMup <- tibble(
  Description  = names(q_dem_for_LSDdn),
  FDR_Dementia = as.numeric(q_dem_for_LSDdn),
  FDR_LSD      = as.numeric(lsd_dn$pmeta[names(q_dem_for_LSDdn)])
) %>%
  left_join(inv_cam, by = "Description") %>%
  filter(!is.na(FDR_LSD), !is.na(FDR_Dementia))

# save tables (useful for debugging / reporting)
write.csv(df_LSDup_DEMdown, file = "Results/Figures/Table_LSDup_vs_DementiaDown_withInVitro.csv", row.names = FALSE)
write.csv(df_LSDdn_DEMup,   file = "Results/Figures/Table_LSDdn_vs_DementiaUp_withInVitro.csv",   row.names = FALSE)

# -----------------------------
# Make plots
# -----------------------------
all_plots <- function(df_dir,
                      out_pdf,
                      # --- in vitro restriction ---
                      invitro_fdr_thr = 0.15,
                      # --- label controls ---
                      n_labels = 8,
                      top_n_each = 120,
                      hard_thr = 0.15,
                      label_max_chars = 70,
                      # --- sizing / layout ---
                      side_in_x = 5.2,
                      side_in_y = 6.2,
                      point_size = 2.2,
                      repel_seed = 1) {
  
  df_dir <- df_dir %>%
    mutate(
      FDR_LSD      = clamp01(FDR_LSD),
      FDR_Dementia = clamp01(FDR_Dementia),
      invitro_FDR  = clamp01(invitro_camera_FDR)
    ) %>%
    filter(!is.na(FDR_LSD), !is.na(FDR_Dementia), !is.na(invitro_FDR))
  
  # =========================================================
  # only in vitro significant
  # =========================================================
  df_plot <- df_dir %>%
    filter(invitro_FDR < invitro_fdr_thr) %>%
    mutate(
      x = -log10(FDR_Dementia),
      y = -log10(FDR_LSD),
      invitro_bin = fdr_bin3(invitro_FDR)
    ) %>%
    filter(!is.na(x), !is.na(y))
  
  if (nrow(df_plot) < 5) {
    stop("Too few points after in vitro filtering.")
  }
  
  df_filt <- df_plot %>%
    filter(FDR_LSD <= hard_thr,
           FDR_Dementia <= hard_thr,
           invitro_FDR <= invitro_fdr_thr)
  
  if (nrow(df_filt) < (n_labels + 5)) df_filt <- df_plot
  
  top_lsd <- df_filt %>% arrange(FDR_LSD)      %>% slice_head(n = top_n_each) %>% pull(Description)
  top_dem <- df_filt %>% arrange(FDR_Dementia) %>% slice_head(n = top_n_each) %>% pull(Description)
  top_inv <- df_filt %>% arrange(invitro_FDR)  %>% slice_head(n = top_n_each) %>% pull(Description)
  
  top_all3 <- Reduce(intersect, list(top_lsd, top_dem, top_inv))
  
  df_scored <- df_filt %>%
    mutate(score = -log10(FDR_LSD) +
             -log10(FDR_Dementia) +
             -log10(invitro_FDR))
  
  lab <- df_scored %>%
    filter(Description %in% top_all3) %>%
    arrange(desc(score)) %>%
    slice_head(n = n_labels)
  
  if (nrow(lab) < n_labels) {
    lab2 <- df_scored %>%
      filter(!Description %in% lab$Description) %>%
      arrange(desc(score)) %>%
      slice_head(n = (n_labels - nrow(lab)))
    lab <- bind_rows(lab, lab2)
  }
  
  lab <- lab %>%
    mutate(label = ifelse(nchar(Description) > label_max_chars,
                          paste0(substr(Description, 1, label_max_chars - 3), "..."),
                          Description))
  
  # =========================================================
  # Plot
  # =========================================================
  
  p <- ggplot(df_plot, aes(x = x, y = y)) +
    
    geom_point(
      aes(fill = invitro_bin),
      shape = 21,
      size = point_size,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    
    ggrepel::geom_label_repel(
      data = lab,
      aes(label = label),
      size = 4,                # <- labels più grandi
      min.segment.length = 0,
      box.padding = 0.5,
      point.padding = 0.35,
      label.size = 0.4,
      max.overlaps = n_labels,
      force = 8,
      force_pull = 0.8,
      seed = repel_seed
    ) +
    
    scale_fill_manual(
      values = c(
        "FDR ≥ 0.15" = "white",
        "FDR < 0.15" = "#FEE0D2",
        "FDR < 0.10" = "#FC9272",
        "FDR < 0.05" = "#DE2D26"
      ),
      name = "In vitro (CAMERA)\nFDR"
    ) +
    
    theme_classic(base_size = 16) +   # <- base più grande
    
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      
      # 🔹 rimuove completamente titolo e sottotitolo
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
      # 🔹 assi più grandi
      axis.title = element_text(size = 14, face = "bold"),
      axis.text  = element_text(size = 14),
      
      # 🔹 legend più leggibile
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 14)
    ) +
    
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    
    labs(
      x = expression(-log[10]~"(FDR Dementia directional meta)"),
      y = expression(-log[10]~"(FDR LSD meta)")
    )
  
  ggsave(file.path("Results/Figures", out_pdf),
         p,
         width = side_in_x,
         height = side_in_y)
  
  invisible(p)
}

all_plots(
  df_LSDup_DEMdown,
  out_pdf  = "Summary_LSDup_vs_DementiaDown_colorInVitro_CAMERA_FDRbins.pdf",
  invitro_fdr_thr = 0.15,
  n_labels = 8,
  top_n_each = 120,
  hard_thr = 0.15
)

all_plots(
  df_LSDdn_DEMup,
  out_pdf  = "Summary_LSDdn_vs_DementiaUp_colorInVitro_CAMERA_FDRbins.pdf",
  invitro_fdr_thr = 0.15,
  n_labels = 8,
  top_n_each = 120,
  hard_thr = 0.15
)

