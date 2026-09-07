#!/usr/bin/env Rscript
# Aggregate per-seed K=1728 saturation trait-recovery GLS summaries
# (CLAMPfull_bp and CLAMPbase) into one long table.

suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) {
    if (i == length(args)) stop("Missing value for ", args[[i]])
    out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", gsub("_", "-", name))
  value
}

args <- parse_cli()
clampfull_dir <- required_arg(args, "clampfull_dir")
clampbase_dir <- required_arg(args, "clampbase_dir")
fdr_cutoff <- as.numeric(required_arg(args, "fdr"))
traits_out <- required_arg(args, "traits_out")

# Filenames are sat_rs<fraction>_k1728_seed<seed>.tsv.gz, optionally suffixed
# with the model name (not consistently applied upstream) -- the directory a
# file lives in, not its filename, is what determines its model here.
scan_dir <- function(dir, model_label) {
  paths <- list.files(
    dir, pattern = "^sat_rs[0-9]+_k1728_seed[0-9]+(_[A-Za-z0-9]+)?\\.tsv\\.gz$", full.names = TRUE
  )
  if (!length(paths)) stop("No K=1728 trait GLS summaries found under ", dir)
  rbindlist(lapply(paths, function(path) {
    m <- regmatches(
      basename(path),
      regexec("^sat_rs([0-9]+)_k1728_seed([0-9]+)(?:_[A-Za-z0-9]+)?\\.tsv\\.gz$", basename(path))
    )[[1]]
    d <- fread(path, select = c("phenotype", "fdr"))
    data.table(
      fraction = as.integer(m[2]),
      seed = as.integer(m[3]),
      model = model_label,
      k = 1728L,
      eligible_traits = uniqueN(d$phenotype),
      recovered_traits = uniqueN(d[fdr < fdr_cutoff, phenotype]),
      fdr = fdr_cutoff
    )
  }))
}

traits <- rbindlist(list(
  scan_dir(clampfull_dir, "CLAMPfull_bp"),
  scan_dir(clampbase_dir, "CLAMPbase")
))
if (anyDuplicated(traits[, .(fraction, seed, model)])) {
  stop("Duplicate fraction/seed/model saturation trait GLS summaries found")
}
setorder(traits, model, fraction, seed)

dir.create(dirname(traits_out), recursive = TRUE, showWarnings = FALSE)
fwrite(traits, traits_out)
