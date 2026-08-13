#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(bigstatsr)
  library(CLAMP)
  library(dplyr)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    key <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) {
    stop("Missing --", gsub("_", "-", name))
  }
  value
}

file_record <- function(path) {
  info <- file.info(path)
  list(path = normalizePath(path), bytes = unname(info$size))
}

output_record <- function(path) {
  info <- file.info(path)
  list(path = basename(path), bytes = unname(info$size))
}

attach_portable_fbm <- function(descriptor_path) {
  fbm <- readRDS(descriptor_path)
  recorded <- fbm$backingfile
  candidates <- unique(c(
    recorded,
    paste0(recorded, ".bk"),
    file.path(dirname(descriptor_path), basename(recorded)),
    file.path(dirname(descriptor_path), paste0(basename(recorded), ".bk"))
  ))
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop(
      "FBM backing file is missing. Recorded path: ", recorded,
      "; descriptor: ", descriptor_path
    )
  }
  if (!file.exists(recorded)) {
    fbm$backingfile <- normalizePath(existing[[1L]])
    message("Reattached FBM backing file at ", fbm$backingfile)
  }
  fbm
}

args <- parse_cli()
dataset <- required_arg(args, "dataset")
prior_gmt <- required_arg(args, "prior_gmt")
svd_path <- required_arg(args, "svd")
k_path <- required_arg(args, "k")
base_path <- required_arg(args, "base_model")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(required_arg(args, "seed"))
max_iter <- as.integer(args$max_iter %||% 5000L)
multiplier <- as.numeric(args$multiplier %||% 100)
fraction <- as.integer(args$fraction %||% 100L)
seed_index <- as.integer(args$seed_index %||% 1L)
requested_mem_mb <- as.integer(args$requested_mem_mb %||% NA_integer_)
requested_cpus <- as.integer(args$requested_cpus %||% NA_integer_)
requested_runtime_min <- as.integer(args$requested_runtime_min %||% NA_integer_)

set.seed(seed)
started <- Sys.time()

if (dir.exists(out_dir)) {
  stop("Refusing to overwrite an existing model directory: ", out_dir)
}
dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
tmp_dir <- paste0(out_dir, ".tmp.", Sys.getpid())
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
dir.create(tmp_dir, recursive = TRUE)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

if (!is.null(args$metadata)) {
  metadata <- readRDS(args$metadata)
  subsample <- readRDS(required_arg(args, "subsample_info"))
  genes <- metadata$gene_symbols_thin
  sample_names <- subsample$sample_names
  backing <- sub("\\.bk$", "", required_arg(args, "fbm_backing"))
  Y <- bigstatsr::FBM(
    nrow = length(genes),
    ncol = length(sample_names),
    backingfile = backing,
    create_bk = FALSE
  )
} else {
  genes <- readRDS(required_arg(args, "genes_rds"))
  sample_names <- readRDS(required_arg(args, "samples_rds"))
  Y <- attach_portable_fbm(required_arg(args, "fbm_rds"))
}

svd_res <- readRDS(svd_path)
clamp_k <- as.integer(readRDS(k_path))
base_res <- readRDS(base_path)

if (!identical(nrow(Y), length(genes)) || !identical(ncol(Y), length(sample_names))) {
  stop(
    "Expression dimensions ", paste(dim(Y), collapse = "x"),
    " do not match ", length(genes), " genes and ", length(sample_names), " samples"
  )
}

gene_sets <- CLAMP:::read_gmt(prior_gmt)
names(gene_sets) <- paste0("GOBP_", names(gene_sets))
path_mat <- CLAMP::gmtListToSparseMat(list(GOBP = gene_sets))
prior_mat <- CLAMP::getMatchedPathwayMat(path_mat, genes)
message(
  "Fitting ", dataset, " fraction=", fraction, " seed_index=", seed_index,
  ": ", nrow(Y), " genes x ", ncol(Y), " samples; K=", clamp_k,
  "; prior sets=", ncol(prior_mat)
)

full_res <- CLAMPfull(
  Y = Y,
  svdres = svd_res,
  priorMat = prior_mat,
  clamp.base.result = base_res,
  use_cpp = TRUE,
  trace = TRUE,
  multiplier = multiplier,
  max.iter = max_iter,
  clamp_k = clamp_k
)

full_res$Z <- data.frame(full_res$Z, check.names = FALSE)
rownames(full_res$Z) <- genes
full_res$B <- data.frame(full_res$B, check.names = FALSE)
colnames(full_res$B) <- sample_names
if (!is.null(full_res$summary)) {
  full_res$summary <- full_res$summary %>%
    rename(LV = LV_index) %>%
    mutate(LV = paste0("LV", LV))
}

if (!identical(dim(full_res$Z), c(length(genes), clamp_k))) {
  stop("Unexpected Z dimensions: ", paste(dim(full_res$Z), collapse = "x"))
}
if (!identical(dim(full_res$B), c(clamp_k, length(sample_names)))) {
  stop("Unexpected B dimensions: ", paste(dim(full_res$B), collapse = "x"))
}
if (anyNA(full_res$Z) || any(!is.finite(as.matrix(full_res$Z)))) {
  stop("Z contains missing or non-finite values")
}
if (anyNA(full_res$B) || any(!is.finite(as.matrix(full_res$B)))) {
  stop("B contains missing or non-finite values")
}

rds_path <- file.path(tmp_dir, "CLAMPfull_bp.rds")
b_path <- file.path(tmp_dir, "B.csv")
z_path <- file.path(tmp_dir, "Z.csv")
summary_path <- file.path(tmp_dir, "summary.csv")
saveRDS(full_res, rds_path)
write.csv(full_res$B, b_path)
write.csv(full_res$Z, z_path)
write.csv(full_res$summary, summary_path, row.names = FALSE)

finished <- Sys.time()
manifest <- list(
  schema_version = 1L,
  status = "complete",
  dataset = dataset,
  fraction = fraction,
  seed_index = seed_index,
  rng_seed = seed,
  genes = length(genes),
  samples = length(sample_names),
  latent_variables = clamp_k,
  prior = list(
    path = normalizePath(prior_gmt),
    md5 = unname(tools::md5sum(prior_gmt)),
    matched_sets = ncol(prior_mat)
  ),
  parameters = list(multiplier = multiplier, max_iter = max_iter),
  requested_resources = list(
    cpus = requested_cpus,
    mem_mb = requested_mem_mb,
    runtime_min = requested_runtime_min
  ),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = unname(as.numeric(difftime(finished, started, units = "secs"))),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  packages = list(
    R = as.character(getRversion()),
    CLAMP = as.character(packageVersion("CLAMP")),
    bigstatsr = as.character(packageVersion("bigstatsr"))
  ),
  inputs = list(
    svd = file_record(svd_path),
    k = file_record(k_path),
    base_model = file_record(base_path)
  ),
  outputs = list(
    model = output_record(rds_path),
    B = output_record(b_path),
    Z = output_record(z_path),
    summary = output_record(summary_path)
  )
)
jsonlite::write_json(manifest, file.path(tmp_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("complete", file.path(tmp_dir, "complete"))

if (!file.rename(tmp_dir, out_dir)) {
  stop("Could not atomically publish ", tmp_dir, " to ", out_dir)
}
published <- TRUE
message("Published validated model: ", out_dir)
