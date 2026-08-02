#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CLAMP)
  library(data.table)
  library(rsvd)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
input <- required_arg(args, "input")
norm_out <- required_arg(args, "norm")
cpm_filt_out <- required_arg(args, "cpm_filt")
stats_out <- required_arg(args, "row_stats")
k_out <- required_arg(args, "k")
diagnostics_out <- required_arg(args, "diagnostics")
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

# norm
counts <- read_csv_matrix(input)
if (!all(is.finite(counts)) || any(counts < 0)) stop("Input contains non-finite or negative values")
cpm <- CLAMP::cpmCLAMP(counts)
prep <- CLAMP::preprocessCLAMP(cpm, mean_cutoff = mean_cutoff, var_cutoff = var_cutoff)
norm <- CLAMP::zscoreCLAMP(prep$Y_filtered, prep$rowStats)

# SVD
n_genes <- nrow(norm)
n_samples <- ncol(norm)
svd_k <- floor((min(n_genes, n_samples) - 1) / 4)
svdres <- rsvd::rsvd(norm, k = svd_k)

# model number of components
k <- CLAMP::num.pc(data = norm, method = "elbow") * 2
k <- max(as.integer(k), 2L)
k <- min(k, svd_k)  # cannot use more singular vectors than computed

diagnostics <- data.frame(
  component = seq_along(svdres$d), singular_value = svdres$d,
  k = k, svd_k = svd_k, seed = seed
)

write_csv_matrix(norm, norm_out)
write_csv_matrix(prep$Y_filtered, cpm_filt_out)
write_csv_matrix(prep$rowStats, stats_out)
ensure_parent(k_out); fwrite(data.frame(k = k), k_out)
ensure_parent(diagnostics_out); fwrite(diagnostics, diagnostics_out)
message("Preprocessed ", nrow(norm), " genes x ", ncol(norm), " samples; k=", k)
