#!/usr/bin/env Rscript
# Fits CLAMPbase on the filtered recount2 matrix.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
fbm_rds <- required_arg(args, "fbm_rds")
genes_rds <- required_arg(args, "genes_rds")
samples_rds <- required_arg(args, "samples_rds")
svd_path <- required_arg(args, "svd")
k_path <- required_arg(args, "k")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

fbm <- readRDS(fbm_rds)
genes <- readRDS(genes_rds)
samples <- readRDS(samples_rds)
svd_res <- readRDS(svd_path)
clamp_k <- as.integer(readRDS(k_path))
if (is.na(clamp_k) || clamp_k < 1L) stop("Invalid CLAMP K: ", clamp_k)

message("CLAMPbase: ", nrow(fbm), " genes x ", ncol(fbm), " samples, K = ", clamp_k)
base_res <- CLAMPbase(Y = fbm, svdres = svd_res, trace = TRUE, clamp_k = clamp_k)

base_res$Z <- data.frame(base_res$Z, check.names = FALSE)
rownames(base_res$Z) <- genes
base_res$B <- data.frame(base_res$B, check.names = FALSE)
colnames(base_res$B) <- samples

if (anyDuplicated(genes)) stop("recount2 gene symbols are not unique")

ensure_dir(out_dir)
model_dir <- ensure_dir(file.path(out_dir, "CLAMPbase"))
saveRDS(base_res, file.path(out_dir, "CLAMPbase.rds"))
write.csv(base_res$B, file.path(model_dir, "B.csv"))
write.csv(base_res$Z, file.path(model_dir, "Z.csv"))
message("Wrote CLAMPbase to ", out_dir)
