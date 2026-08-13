#!/usr/bin/env Rscript
#
# Randomized SVD of the filtered ARCHS4 FBM, plus the CLAMP rank estimate.
#
# Extracted from nbs/01_model_building/04_archs4/02_archs4_svd.ipynb, which does
# not run as committed: it referenced an undefined `SVD_K`, read the core count
# from `config$ARCHS4$PLIER_PARAMS` (the block is called CLAMP_PARAMS), and never
# loaded bigstatsr. `--svd-k` is now an explicit argument.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
meta_path <- required_arg(args, "metadata")
fbm_backingfile <- required_arg(args, "fbm")
svd_out <- required_arg(args, "svd")
k_out <- required_arg(args, "k")
# Optional: the raw decomposition before NaN singular values are dropped. Tens
# of GB, and not retained for the shipped model, so it is written only on request.
svd_full_out <- args$svd_full
n_cores <- as.integer(args$n_cores %||% 1L)
k_multiplier <- as.integer(args$k_multiplier %||% 2L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

meta <- readRDS(meta_path)
n_genes <- meta$n_genes_thin
n_samples <- meta$n_samples

fbm <- attach_fbm(fbm_backingfile, nrow = n_genes, ncol = n_samples)

# Default matches scripts/pseudobulk/preprocess.R: as many components as the
# smaller dimension supports, capped so the randomized solver stays tractable.
svd_k <- as.integer(args$svd_k %||% max(floor((min(n_genes, n_samples) - 1L) / 4L), 2L))
message("Randomized SVD: ", n_genes, " x ", n_samples, ", k = ", svd_k, ", cores = ", n_cores)

svd_res <- with_single_threaded_blas(
  n_cores,
  big_randomSVD(fbm, k = svd_k, ncores = n_cores)
)
if (!is.null(svd_full_out) && is.character(svd_full_out)) {
  saveRDS(svd_res, ensure_parent(svd_full_out))
}

# big_randomSVD can return NaN singular values for the trailing components on a
# matrix this shape; drop them before anything downstream ranks on `d`.
valid_idx <- which(!is.nan(svd_res$d))
if (length(valid_idx) < length(svd_res$d)) {
  message("Dropping ", length(svd_res$d) - length(valid_idx), " NaN singular values")
}
svd_res$d <- svd_res$d[valid_idx]
svd_res$u <- svd_res$u[, valid_idx, drop = FALSE]
svd_res$v <- svd_res$v[, valid_idx, drop = FALSE]
saveRDS(svd_res, ensure_parent(svd_out))

clamp_k <- as.integer(num.pc(list(d = svd_res$d)) * k_multiplier)
message("Inferred CLAMP K = ", clamp_k)
saveRDS(clamp_k, ensure_parent(k_out))
