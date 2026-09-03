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

# Projection workflow helpers -------------------------------------------------

read_projection_model <- function(path) {
  mdl <- readRDS(path)
  keep <- list(Z = as.matrix(mdl$Z), L2 = mdl$L2, summary = mdl$summary)
  rm(mdl)
  gc()
  if (is.null(keep$Z)) stop("Model has no Z matrix: ", path)
  if (is.null(keep$L2)) stop("Model has no L2 value: ", path)
  keep
}

permute_model_z <- function(mdl, seed) {
  set.seed(seed)
  gene_names <- rownames(mdl$Z)
  mdl$Z <- mdl$Z[sample(nrow(mdl$Z)), , drop = FALSE]
  rownames(mdl$Z) <- gene_names
  mdl
}

begin_publish <- function(out_dir) {
  if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir)
  ensure_dir(dirname(out_dir))
  tmp_dir <- paste0(out_dir, ".tmp.", Sys.getpid())
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  ensure_dir(tmp_dir)
  tmp_dir
}

finish_publish <- function(tmp_dir, out_dir) {
  writeLines("complete", file.path(tmp_dir, "complete"))
  if (!file.rename(tmp_dir, out_dir)) stop("Could not atomically publish: ", out_dir)
  message("Published: ", out_dir)
  invisible(out_dir)
}

file_md5 <- function(path) unname(tools::md5sum(normalizePath(path, mustWork = FALSE)))

write_manifest <- function(dir, method, parameters = list(), inputs = list(),
                           extra = list(), started = NULL) {
  finished <- Sys.time()
  manifest <- c(
    list(
      schema_version = 1L,
      status = "complete",
      method = method,
      backend = if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "slurm" else "local",
      inputs = lapply(inputs, function(p) list(path = normalizePath(p, mustWork = FALSE),
                                               md5 = file_md5(p))),
      parameters = parameters
    ),
    extra,
    list(
      started_utc = if (is.null(started)) NA_character_ else format(started, tz = "UTC", usetz = TRUE),
      finished_utc = format(finished, tz = "UTC", usetz = TRUE),
      elapsed_seconds = if (is.null(started)) NA_real_ else as.numeric(difftime(finished, started, units = "secs")),
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
      package_versions = list(R = as.character(getRversion()),
                              CLAMP = as.character(utils::packageVersion("CLAMP")))
    )
  )
  jsonlite::write_json(manifest, file.path(dir, "manifest.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(manifest)
}

read_matrix_csv <- function(path) {
  dt <- data.table::fread(path, check.names = FALSE)
  m <- as.matrix(dt[, -1L])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(dt[[1L]])
  m
}

write_matrix_csv <- function(m, path) {
  out <- data.table::data.table(rn = rownames(m))
  data.table::setnames(out, "rn", "")
  data.table::fwrite(cbind(out, data.table::as.data.table(m)), path)
  invisible(path)
}

apply_annotation_rules <- function(x, spec) {
  df <- data.frame(.source = x, stringsAsFactors = FALSE)
  for (nm in names(spec$rules)) {
    rule <- spec$rules[[nm]]
    df[[nm]] <- if (!is.null(rule$join)) {
      do.call(paste, c(lapply(unlist(rule$join), function(k) df[[k]]),
                       list(sep = rule$sep %||% "_")))
    } else if (!is.null(rule$extract)) {
      sub(rule$extract, "\\1", x)
    } else if (!is.null(rule$pattern)) {
      ifelse(grepl(rule$pattern, x), rule$match, rule$nomatch)
    } else {
      stop("Unrecognised annotation rule for column: ", nm)
    }
  }
  df$sample <- if (!is.null(spec$id)) {
    vapply(seq_len(nrow(df)), function(i) {
      s <- spec$id
      for (nm in names(df)) s <- gsub(paste0("{", nm, "}"), df[[nm]][i], s, fixed = TRUE)
      s
    }, character(1))
  } else {
    x
  }
  df
}

ensembl_to_symbol <- function(counts) {
  ids <- sub("\\..*", "", rownames(counts))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids),
                               column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  vec <- unname(sym[ids])
  keep <- !is.na(vec) & nzchar(vec)
  counts <- counts[keep, , drop = FALSE]
  rownames(counts) <- vec[keep]
  rowsum(counts, rownames(counts))
}

normalize_counts <- function(counts, mean_cutoff, var_cutoff) {
  cpm <- CLAMP::cpmCLAMP(counts)
  prep <- CLAMP::preprocessCLAMP(Y = cpm, mean_cutoff = mean_cutoff, var_cutoff = var_cutoff)
  CLAMP::zscoreCLAMP(Y_filtered = prep$Y_filtered, rowStats = prep$rowStats)
}
