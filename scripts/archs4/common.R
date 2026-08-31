#!/usr/bin/env Rscript

# Shared ARCHS4 helpers.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Parses named command-line arguments.
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

# Returns a required command-line argument.
required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) stop("Missing --", gsub("_", "-", name))
  value
}

# Creates a directory when it does not exist.
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

# Creates the parent directory for a file.
ensure_parent <- function(path) {
  ensure_dir(dirname(path))
  path
}

# Runs an expression with one BLAS thread per worker.
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

# Attaches an existing FBM backing file.
attach_fbm <- function(backingfile, nrow, ncol) {
  bigstatsr::FBM(nrow = nrow, ncol = ncol, backingfile = backingfile, create_bk = FALSE)
}

# Converts sample-by-gene counts to TPM.
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

# Builds a CLAMPfull prior from GMT collections.
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

# Writes a CLAMP result in the standard model layout.
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
