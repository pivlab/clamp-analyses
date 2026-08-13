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

# Log-transform before filtering, matching the GTEx pipeline.
#
# The two pipelines enter CLAMP through different doors: GTEx builds an FBM and
# calls preprocessCLAMPFBM, which runs cleanFBM and applies log2(x + 1) whenever
# the matrix maximum is >= 100.  Pseudobulk calls the in-memory preprocessCLAMP,
# which does no transform at all and assumes the caller passes data on a suitable
# scale.  Without this line the pipelines differed in two ways: GTEx modelled
# log2 data while pseudobulk modelled raw CPM, and -- more subtly -- the shared
# mean_cutoff/var_cutoff were applied to incomparable scales, so identical config
# values selected genes by different criteria (linear CPM row variances here run
# to 1e7, against a var_cutoff of 0.1).
#
# cleanFBM's >= 100 condition is always true for CPM, so it is applied
# unconditionally rather than reproducing a branch that never takes the else.
log_cpm <- log2(cpm + 1)

prep <- CLAMP::preprocessCLAMP(log_cpm, mean_cutoff = mean_cutoff, var_cutoff = var_cutoff)
norm <- CLAMP::zscoreCLAMP(prep$Y_filtered, prep$rowStats)

# SVD
n_genes <- nrow(norm)
n_samples <- ncol(norm)
svd_k <- max(floor((min(n_genes, n_samples) - 1) / 4), 2L)
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
# NOTE: despite the file name, this is now the filtered log2(CPM + 1) matrix, not
# raw CPM.  It is the pre-z-score matrix CoGAPS consumes, and CoGAPS no longer
# log-transforms it itself.  The name is kept because renaming it ripples through
# pseudobulk.smk, donor_bulk.smk, computational_timing.smk and the timing driver.
write_csv_matrix(prep$Y_filtered, cpm_filt_out)
write_csv_matrix(prep$rowStats, stats_out)
ensure_parent(k_out); fwrite(data.frame(k = k), k_out)
ensure_parent(diagnostics_out); fwrite(diagnostics, diagnostics_out)
message("Preprocessed ", nrow(norm), " genes x ", ncol(norm), " samples; k=", k)
