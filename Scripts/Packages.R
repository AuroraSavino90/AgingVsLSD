library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Rn.eg.db)   
library(FactoMineR)
library(sva)
library(biomaRt)
library(openxlsx)
library(DESeq2)
library(ggrepel)
library(metap)
library(org.Hs.eg.db)
library(slinky)
library(GSVA)
library(GEOquery)

# Save versions of loaded packages
pkgs <- sessionInfo()$otherPkgs
pkg_versions <- data.frame(
  Package = names(pkgs),
  Version = sapply(pkgs, function(x) x$Version)
)

write.csv(pkg_versions, "package_versions.csv", row.names = FALSE)
