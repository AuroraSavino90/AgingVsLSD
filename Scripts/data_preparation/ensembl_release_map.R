################################################################################
# ensembl_release_map.R
#
# Ensembl gene ID -> gene symbol table of a fixed Ensembl release, identical to
#   getBM(attributes = c("ensembl_gene_id", "<symbol>_symbol"),
#         mart = useEnsembl("ensembl", "<species>_gene_ensembl", version = <release>))
#
# The original scripts queried the live Ensembl BioMart (useMart("ensembl")),
# so the symbols depend on the release that was current at the time
#
# The BioMart archive sites of old releases are often unreachable, so the table
# is rebuilt from the BioMart database dump of that release on the Ensembl FTP
# (permanent, https://ftp.ensembl.org/pub/release-<N>/mysql/ensembl_mart_<N>/):
#   <species>_gene_ensembl__gene__main.txt.gz      gene key -> Ensembl gene ID
#   <species>_gene_ensembl__ox_<db>__dm.txt.gz     gene key -> symbol (HGNC/MGI/RGD)
# Genes without a symbol get "" (as in getBM). The files are downloaded once
# into CACHE_DIR and the table is cached there as a tab-separated file.
################################################################################

ensembl_release_map <- function(release, species = c("hsapiens", "mmusculus", "rnorvegicus"),
                                cache_dir = "Data/ensembl") {
  species <- match.arg(species)
  db  <- c(hsapiens = "hgnc", mmusculus = "mgi", rnorvegicus = "rgd")[[species]]
  col <- paste0(db, "_symbol")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(cache_dir, sprintf("ensembl%d_%s_%s.tsv", release, species, col))
  if (file.exists(out))
    return(read.delim(out, colClasses = "character", na.strings = character(0)))

  get_table <- function(table) {
    f <- file.path(cache_dir, sprintf("ensembl_mart_%d_%s.txt.gz", release, table))
    if (!file.exists(f)) {
      old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
      download.file(sprintf("https://ftp.ensembl.org/pub/release-%d/mysql/ensembl_mart_%d/%s.txt.gz",
                            release, release, table), f, mode = "wb", quiet = TRUE)
    }
    read.delim(gzfile(f), header = FALSE, quote = "", comment.char = "",
               na.strings = "\\N", colClasses = "character")
  }
  genes <- get_table(paste0(species, "_gene_ensembl__gene__main"))   # col 2 = key, col 7 = stable ID
  xref  <- get_table(paste0(species, "_gene_ensembl__ox_", db, "__dm"))  # col 1 = key, col 3 = symbol
  xref  <- xref[!is.na(xref[, 3]), c(1, 3)]

  sym <- split(xref[, 2], factor(xref[, 1], levels = unique(xref[, 1])))
  s   <- sym[genes[, 2]]
  n   <- lengths(s)
  map <- data.frame(ensembl_gene_id = rep(genes[, 7], pmax(n, 1)),
                    symbol = unlist(lapply(s, function(v) if (length(v)) v else ""), use.names = FALSE),
                    stringsAsFactors = FALSE)
  names(map)[2] <- col
  map <- unique(map)
  rownames(map) <- NULL
  write.table(map, out, sep = "\t", quote = FALSE, row.names = FALSE)
  map
}
