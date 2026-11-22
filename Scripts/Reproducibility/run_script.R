# ============================================================
# Run ONE original pipeline script, unmodified, inside a sandbox.
#
#   Rscript Scripts/Reproducibility/run_script.R <sandbox_dir> <script.R> [seed]
#
# - working directory = <sandbox_dir> (contains Data/ junction, Results/ copy,
#   Scripts/ copy), so nothing in the repository is written
# - setwd("workdir") in the scripts is turned into a no-op via a shim attached
#   to the search path (survives the scripts' rm(list = ls()))
# - the script is sourced with print.eval = TRUE, i.e. top-level ggplot objects
#   are printed as when running the code interactively in RStudio
# - optional [seed]: set.seed(seed) right before sourcing (no code edits)
# - writes <sandbox>/logs/<script>.status (OK / ERROR + message + runtime) and
#   <sandbox>/logs/<script>.sessionInfo.txt
# ============================================================

# everything runs inside local(): the scripts' rm(list = ls()) only clears globalenv
status <- local({
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) >= 2)
  sandbox <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
  script  <- args[2]
  seed    <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else NA

  base::setwd(sandbox)
  dir.create("logs", showWarnings = FALSE)
  stem <- sub("\\.R$", "", basename(script))

  shims <- new.env()
  shims$setwd <- function(dir) {
    if (identical(dir, "workdir")) invisible(getwd()) else base::setwd(dir)
  }
  attach(shims, name = "repro_shims", warn.conflicts = FALSE)

  # ---- optional, opt-in shims (environment variables) ----------------------
  # They patch package NAMESPACES before the scripts call library(), so the
  # scripts themselves stay unmodified.
  patch <- function(pkg, fun, wrapper) {
    ns <- asNamespace(pkg)
    orig <- get(fun, envir = ns)
    unlockBinding(fun, ns)
    assign(fun, wrapper(orig), envir = ns)
    lockBinding(fun, ns)
  }
  # REPRO_ENSEMBL_HOST: pin biomaRt queries to an archived Ensembl release
  # (e.g. https://may2025.archive.ensembl.org) instead of the live site
  ens_host <- Sys.getenv("REPRO_ENSEMBL_HOST")
  if (nzchar(ens_host)) {
    patch("biomaRt", "useMart", function(orig) function(biomart, dataset, host, ...) {
      if (identical(biomart, "ensembl")) biomart <- "ENSEMBL_MART_ENSEMBL"
      orig(biomart = biomart, dataset = dataset, host = ens_host, ...)
    })
    patch("biomaRt", "useEnsembl", function(orig) function(biomart, dataset, host, mirror, version, ...)
      orig(biomart = biomart, dataset = dataset, host = ens_host, ...))
    message("[repro] biomaRt pinned to ", ens_host)
  }
  # REPRO_PRELOAD: comma-separated packages attached before the script runs
  # (e.g. data.table: 12_dementia_cohen.R calls data.table() before library(data.table),
  # which only works in a session where data.table is already attached)
  pre <- strsplit(Sys.getenv("REPRO_PRELOAD"), ",")[[1]]
  for (p in pre[nzchar(pre)]) suppressPackageStartupMessages(library(p, character.only = TRUE))
  if (length(pre)) message("[repro] preloaded: ", paste(pre, collapse = ", "))

  # REPRO_TRACE_FGSEA=1: record .Random.seed (and gene-set sizes) before every fgsea()
  # call, saved to logs/<script>_fgsea_trace.rds; results are not changed
  if (Sys.getenv("REPRO_TRACE_FGSEA") == "1") {
    trace_env <- new.env(); trace_env$calls <- list()
    patch("fgsea", "fgsea", function(orig) function(pathways, stats, ...) {
      rs <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
      trace_env$calls[[length(trace_env$calls) + 1]] <- list(seed = rs, n_path = length(pathways),
        sizes = lengths(pathways), n_stats = length(stats))
      orig(pathways, stats, ...)
    })
    reg.finalizer(trace_env, function(e) NULL)
    assign("repro_trace_env", trace_env, envir = shims)
    message("[repro] tracing RNG state before every fgsea() call")
  }
  # REPRO_SEED: set.seed(REPRO_SEED) right before EVERY fgsea / gseGO call,
  # so each stochastic GSEA is reproducible regardless of call order
  rs <- Sys.getenv("REPRO_SEED")
  if (nzchar(rs)) {
    rs <- as.integer(rs)
    patch("fgsea", "fgsea", function(orig) function(...) { set.seed(rs); orig(...) })
    if (requireNamespace("clusterProfiler", quietly = TRUE))
      patch("clusterProfiler", "gseGO", function(orig) function(...) { set.seed(rs); orig(...) })
    message("[repro] set.seed(", rs, ") before every fgsea()/gseGO() call")
  }

  t0 <- Sys.time()
  if (!is.na(seed)) set.seed(seed)
  res <- tryCatch({
    source(file.path("Scripts", script), print.eval = TRUE, echo = FALSE,
           encoding = "UTF-8", max.deparse.length = 200)
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))

  while (dev.cur() > 1) dev.off()
  if (exists("repro_trace_env", envir = shims))
    saveRDS(shims$repro_trace_env$calls, file.path(sandbox, "logs", paste0(stem, "_fgsea_trace.rds")))
  # keep the small objects the script left in the workspace (statistics quoted in
  # the text are often computed but never written to a file)
  try({
    objs <- ls(globalenv())
    small <- objs[vapply(objs, function(o) as.numeric(object.size(get(o, envir = globalenv()))) < 50e6, logical(1))]
    save(list = small, envir = globalenv(),
         file = file.path(sandbox, "logs", paste0(stem, "_objects.RData")))
  }, silent = TRUE)
  mins <-round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  writeLines(c(res, paste("minutes:", mins), paste("seed:", seed)),
             file.path(sandbox, "logs", paste0(stem, ".status")))
  writeLines(capture.output(sessionInfo()),
             file.path(sandbox, "logs", paste0(stem, ".sessionInfo.txt")))
  cat("\n[run_script]", script, "->", res, "(", mins, "min )\n")
  res
})
quit(status = if (identical(status, "OK")) 0 else 1)
