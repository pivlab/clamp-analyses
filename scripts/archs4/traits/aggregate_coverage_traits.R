#!/usr/bin/env Rscript
# Aggregate per-seed coverage trait-recovery GLS summaries (CLAMPfull_bp and
# CLAMPbase, at each model's native full rank) into one long table.

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
finals_clampfull_dir <- required_arg(args, "finals_clampfull_dir")
finals_clampbase_dir <- required_arg(args, "finals_clampbase_dir")
fdr_cutoff <- as.numeric(required_arg(args, "fdr"))
traits_out <- required_arg(args, "traits_out")
finals_out <- required_arg(args, "finals_out")

# Filenames are cov_rs<fraction>_seed<seed>.tsv.gz, optionally suffixed with
# the model name (not consistently applied upstream) -- the directory a file
# lives in, not its filename, is what determines its model here.
scan_dir <- function(dir, model_label) {
  paths <- list.files(
    dir, pattern = "^cov_rs[0-9]+_seed[0-9]+(_[A-Za-z0-9]+)?\\.tsv\\.gz$", full.names = TRUE
  )
  if (!length(paths)) stop("No coverage trait GLS summaries found under ", dir)
  rbindlist(lapply(paths, function(path) {
    m <- regmatches(
      basename(path),
      regexec("^cov_rs([0-9]+)_seed([0-9]+)(?:_[A-Za-z0-9]+)?\\.tsv\\.gz$", basename(path))
    )[[1]]
    d <- fread(path, select = c("phenotype", "fdr"))
    data.table(
      fraction = as.integer(m[2]),
      seed = as.integer(m[3]),
      model = model_label,
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
  stop("Duplicate fraction/seed/model coverage trait GLS summaries found")
}
setorder(traits, model, fraction, seed)

dir.create(dirname(traits_out), recursive = TRUE, showWarnings = FALSE)
fwrite(traits, traits_out)

# Full-data (100%) fits published per compendium: no fractions/seeds, one
# GLS summary per (dataset, model).
scan_finals_dir <- function(dir, model_label) {
  paths <- list.files(
    dir, pattern = "^final_(archs4|gtex|recount2)(_[A-Za-z0-9]+)?\\.tsv\\.gz$", full.names = TRUE
  )
  if (!length(paths)) stop("No final-model trait GLS summaries found under ", dir)
  rbindlist(lapply(paths, function(path) {
    m <- regmatches(
      basename(path),
      regexec("^final_(archs4|gtex|recount2)(?:_[A-Za-z0-9]+)?\\.tsv\\.gz$", basename(path))
    )[[1]]
    d <- fread(path, select = c("phenotype", "fdr"))
    data.table(
      dataset = m[2],
      model = model_label,
      eligible_traits = uniqueN(d$phenotype),
      recovered_traits = uniqueN(d[fdr < fdr_cutoff, phenotype]),
      fdr = fdr_cutoff
    )
  }))
}

finals <- rbindlist(list(
  scan_finals_dir(finals_clampfull_dir, "CLAMPfull_bp"),
  scan_finals_dir(finals_clampbase_dir, "CLAMPbase")
))
if (anyDuplicated(finals[, .(dataset, model)])) {
  stop("Duplicate dataset/model final-model trait GLS summaries found")
}
setorder(finals, model, dataset)

dir.create(dirname(finals_out), recursive = TRUE, showWarnings = FALSE)
fwrite(finals, finals_out)
