#!/usr/bin/env Rscript

# Fits ARCHS4 CLAMPbase or CLAMPfull models.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(dplyr)
  library(Matrix)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
meta_path <- required_arg(args, "metadata")
samples_path <- required_arg(args, "samples")
fbm_backingfile <- required_arg(args, "fbm")
svd_path <- required_arg(args, "svd")
k_path <- required_arg(args, "k")
stage <- required_arg(args, "stage")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

if (!stage %in% c("base", "full")) stop("--stage must be 'base' or 'full', got: ", stage)

meta <- readRDS(meta_path)
genes <- meta$gene_symbols_thin
n_genes <- meta$n_genes_thin
n_samples <- meta$n_samples

all_samples <- readRDS(samples_path)
sample_names <- all_samples[seq_len(n_samples)]

fbm <- attach_fbm(fbm_backingfile, nrow = n_genes, ncol = n_samples)
svd_res <- readRDS(svd_path)
clamp_k <- readRDS(k_path)
message("CLAMP", stage, ": ", n_genes, " genes x ", n_samples, " samples, K = ", clamp_k)

ensure_dir(out_dir)

if (stage == "base") {
  base_res <- CLAMPbase(
    Y = fbm,
    svdres = svd_res,
    trace = TRUE,
    clamp_k = clamp_k
  )
  write_clamp_model(
    base_res,
    model_dir = file.path(out_dir, "CLAMPbase"),
    genes = genes,
    sample_names = sample_names,
    rds_path = file.path(out_dir, "CLAMPbase.rds")
  )
  message("Wrote CLAMPbase to ", out_dir)
} else {
  base_rds <- required_arg(args, "base_model")
  prior_paths <- strsplit(required_arg(args, "prior_gmts"), ",", fixed = TRUE)[[1]]
  prior_names <- strsplit(required_arg(args, "prior_names"), ",", fixed = TRUE)[[1]]
  stopifnot(length(prior_paths) == length(prior_names))
  multiplier <- as.numeric(args$multiplier %||% 100)
  max_iter <- as.integer(args$max_iter %||% 5000L)

  gmts <- setNames(as.list(prior_paths), prior_names)
  prior <- build_prior(genes, gmts)
  saveRDS(prior$path_mat, file.path(out_dir, "pathMat.rds"))
  saveRDS(prior$matched, file.path(out_dir, "pathMat_matched.rds"))
  message("Prior matched: ", nrow(prior$matched), " genes x ", ncol(prior$matched), " gene sets")

  base_res <- readRDS(base_rds)
  full_res <- CLAMPfull(
    Y = fbm,
    svdres = svd_res,
    priorMat = prior$matched,
    clamp.base.result = base_res,
    use_cpp = TRUE,
    trace = TRUE,
    multiplier = multiplier,
    max.iter = max_iter,
    clamp_k = clamp_k
  )
  write_clamp_model(
    full_res,
    model_dir = file.path(out_dir, "CLAMPfull_hall"),
    genes = genes,
    sample_names = sample_names,
    rds_path = file.path(out_dir, "CLAMPfull_hall.rds")
  )
  message("Wrote CLAMPfull to ", out_dir)
}
