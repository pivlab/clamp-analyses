#!/usr/bin/env Rscript

# Fit CLAMPbase at a forced rank K for the saturation sweep.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(Matrix)
  library(jsonlite)
  library(CLAMP)
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

args <- parse_cli()
svd_path <- required_arg(args, "svd")
metadata_path <- required_arg(args, "metadata")
subsample_path <- required_arg(args, "subsample_info")
fbm_backing <- required_arg(args, "fbm_backing")
out_dir <- required_arg(args, "out_dir")
clamp_k <- as.integer(required_arg(args, "clamp_k"))
seed <- as.integer(required_arg(args, "seed"))
fraction <- as.integer(args$fraction %||% NA_integer_)
seed_index <- as.integer(args$seed_index %||% 1L)
requested_mem_mb <- as.integer(args$requested_mem_mb %||% NA_integer_)
requested_cpus <- as.integer(args$requested_cpus %||% NA_integer_)
requested_runtime_min <- as.integer(args$requested_runtime_min %||% NA_integer_)

if (is.na(clamp_k) || clamp_k < 1L) stop("Invalid CLAMP K: ", clamp_k)
set.seed(seed)
started <- Sys.time()

if (file.exists(file.path(out_dir, "CLAMPbase.rds"))) {
  stop("Refusing to overwrite an existing CLAMPbase: ", out_dir)
}
dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
tmp_dir <- paste0(out_dir, ".tmp.", Sys.getpid())
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
dir.create(file.path(tmp_dir, "CLAMPbase"), recursive = TRUE)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

metadata <- readRDS(metadata_path)
subsample <- readRDS(subsample_path)
genes <- metadata$gene_symbols_thin
sample_names <- subsample$sample_names
backing <- sub("\\.bk$", "", fbm_backing)
Y <- bigstatsr::FBM(
  nrow = length(genes),
  ncol = length(sample_names),
  backingfile = backing,
  create_bk = FALSE
)

svd_res <- readRDS(svd_path)

if (!identical(nrow(Y), length(genes)) || !identical(ncol(Y), length(sample_names))) {
  stop(
    "Expression dimensions ", paste(dim(Y), collapse = "x"),
    " do not match ", length(genes), " genes and ", length(sample_names), " samples"
  )
}

message(
  "CLAMPbase saturation fraction=", fraction, " K=", clamp_k,
  " seed_index=", seed_index, ": ",
  nrow(Y), " genes x ", ncol(Y), " samples"
)

base_res <- CLAMPbase(
  Y = Y,
  svdres = svd_res,
  trace = TRUE,
  clamp_k = clamp_k
)

if (!identical(dim(base_res$Z), c(length(genes), clamp_k))) {
  stop("Unexpected Z dimensions: ", paste(dim(base_res$Z), collapse = "x"))
}
if (!identical(dim(base_res$B), c(clamp_k, length(sample_names)))) {
  stop("Unexpected B dimensions: ", paste(dim(base_res$B), collapse = "x"))
}
if (anyNA(base_res$Z) || any(!is.finite(as.matrix(base_res$Z)))) {
  stop("Z contains missing or non-finite values")
}
if (anyNA(base_res$B) || any(!is.finite(as.matrix(base_res$B)))) {
  stop("B contains missing or non-finite values")
}

rownames(base_res$Z) <- genes
colnames(base_res$B) <- sample_names

saveRDS(base_res, file.path(tmp_dir, "CLAMPbase.rds"))
write.csv(base_res$Z, file.path(tmp_dir, "CLAMPbase", "Z.csv"))
write.csv(base_res$B, file.path(tmp_dir, "CLAMPbase", "B.csv"))

finished <- Sys.time()
manifest <- list(
  schema_version = 1L,
  status = "complete",
  model = "CLAMPbase",
  dataset = "archs4",
  fraction = fraction,
  seed_index = seed_index,
  rng_seed = seed,
  genes = length(genes),
  samples = length(sample_names),
  latent_variables = clamp_k,
  prior = NULL,
  requested_resources = list(
    mem_mb = requested_mem_mb,
    cpus = requested_cpus,
    runtime_min = requested_runtime_min
  ),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  package_versions = list(
    R = as.character(getRversion()),
    CLAMP = as.character(packageVersion("CLAMP"))
  )
)
jsonlite::write_json(
  manifest, file.path(tmp_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null"
)

if (dir.exists(out_dir)) {
  for (entry in list.files(tmp_dir, all.files = TRUE, no.. = TRUE)) {
    file.rename(file.path(tmp_dir, entry), file.path(out_dir, entry))
  }
  unlink(tmp_dir, recursive = TRUE)
} else if (!file.rename(tmp_dir, out_dir)) {
  stop("Could not atomically publish CLAMPbase directory")
}
published <- TRUE
message("Published CLAMPbase: ", out_dir)
