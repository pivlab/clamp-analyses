# Shared ARCHS4 helpers used by the production scripts.
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
  if (is.null(value) || !nzchar(as.character(value))) stop("Missing --", gsub("_", "-", name))
  value
}

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

ensure_parent <- function(path) {
  ensure_dir(dirname(path))
  path
}

# bigstatsr parallelises with forked workers, which deadlock against a threaded
# BLAS. Wrap any multi-core bigstatsr call in this rather than repeating the
# option dance at every call site (the notebooks repeated it four times).
with_single_threaded_blas <- function(n_cores, expr) {
  if (n_cores > 1) {
    old_check <- getOption("bigstatsr.check.parallel.blas")
    old_nproc <- getOption("default.nproc.blas")
    options(bigstatsr.check.parallel.blas = FALSE, default.nproc.blas = NULL)
    on.exit({
      options(bigstatsr.check.parallel.blas = old_check, default.nproc.blas = old_nproc)
    }, add = TRUE)
  }
  force(expr)
}

# Attach an existing FBM backing file (.bk) without recreating it.
attach_fbm <- function(backingfile, nrow, ncol) {
  bigstatsr::FBM(nrow = nrow, ncol = ncol, backingfile = backingfile, create_bk = FALSE)
}

# TPM from a counts matrix laid out samples x genes, using a named vector of
# gene lengths in bp.
#
# Why TPM here and not CLAMP's cpmCLAMP(), which the pseudobulk pipeline uses:
#   * This is what built every ARCHS4 model currently on disk. Both the main and
#     bp_coverage copies of 01_archs4_preprocess.ipynb normalize with TPM and
#     neither ever calls cpmCLAMP(), so changing it would mean refitting the
#     whole compendium rather than adopting the existing fits.
#   * ARCHS4 ships raw gene counts, whereas GTEx enters the repo already in TPM
#     (the download is gene_tpm.gct.gz) and gets no CPM step -- scripts/gtex/
#     clamp.R goes straight to preprocessCLAMPFBM() + zscoreCLAMPFBM(). TPM
#     therefore lands ARCHS4 on the same footing as GTEx.
#   * CLAMP ships no TPM equivalent (only cpmCLAMP/cpmCLAMPFBM), so there is
#     nothing in the package to call here.
# Note the choice matters less than it looks: both pipelines z-score per gene
# afterwards, which cancels the constant per-gene length factor. What TPM and
# CPM still disagree on is the per-sample denominator (length-corrected library
# size vs raw library size), and that does survive z-scoring.
#
# `gene_lengths` must already be aligned to the columns of `counts`: scaling is
# positional, so a mismatch would silently divide each gene by another gene's
# length.
tpm_norm <- function(counts, gene_lengths) {
  if (!is.matrix(counts)) stop("`counts` must be a matrix.")
  if (is.null(names(gene_lengths))) stop("`gene_lengths` must be a named numeric vector.")
  if (ncol(counts) != length(gene_lengths)) {
    stop("`counts` has ", ncol(counts), " columns but `gene_lengths` has ", length(gene_lengths), " entries.")
  }
  if (!is.null(colnames(counts)) && !identical(colnames(counts), names(gene_lengths))) {
    stop("`counts` columns are not in the same order as `gene_lengths`.")
  }

  lengths_kb <- gene_lengths / 1e3
  rpk <- sweep(counts, 2, lengths_kb, FUN = "/")
  per_sample_sum <- rowSums(rpk)
  sweep(rpk, 1, per_sample_sum, FUN = "/") * 1e6
}

# Build the CLAMPfull prior from several GMT collections, prefixing each set
# name with its collection so names stay unique after the merge.
# gmts: named character vector, collection -> gmt path.
build_prior <- function(genes, gmts) {
  gmt_list <- lapply(names(gmts), function(collection) {
    sets <- CLAMP:::read_gmt(gmts[[collection]])
    names(sets) <- paste0(collection, "_", names(sets))
    sets
  })
  names(gmt_list) <- names(gmts)

  path_mat <- CLAMP::gmtListToSparseMat(gmt_list)
  matched <- CLAMP::getMatchedPathwayMat(path_mat, genes)
  list(path_mat = path_mat, matched = matched)
}

# CLAMP result -> the on-disk model layout every downstream analysis expects:
# <model_dir>/{B,Z,summary}.csv plus the serialized object beside it.
write_clamp_model <- function(res, model_dir, genes, sample_names, rds_path = NULL) {
  res$Z <- data.frame(res$Z)
  rownames(res$Z) <- genes
  res$B <- data.frame(res$B)
  colnames(res$B) <- sample_names

  if (!is.null(res$summary)) {
    res$summary <- dplyr::mutate(
      dplyr::rename(res$summary, LV = LV_index),
      LV = paste0("LV", LV)
    )
  }

  ensure_dir(model_dir)
  write.csv(res$B, file.path(model_dir, "B.csv"))
  write.csv(res$Z, file.path(model_dir, "Z.csv"))
  if (!is.null(res$summary)) {
    write.csv(res$summary, file.path(model_dir, "summary.csv"))
  }
  if (!is.null(rds_path)) saveRDS(res, ensure_parent(rds_path))

  invisible(res)
}
