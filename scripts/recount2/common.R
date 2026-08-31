# Shared recount2 helpers.

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
  if (is.null(value) || !nzchar(as.character(value))) {
    stop("Missing --", gsub("_", "-", name))
  }
  value
}

# Creates a directory when it does not exist.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Attaches an existing FBM backing file.
attach_fbm <- function(backingfile, nrow, ncol) {
  bigstatsr::FBM(nrow = nrow, ncol = ncol, backingfile = backingfile, create_bk = FALSE)
}

# Runs an expression with one BLAS thread per worker.
with_single_threaded_blas <- function(n_cores, expr) {
  if (n_cores <= 1L) return(force(expr))
  old_check <- getOption("bigstatsr.check.parallel.blas")
  old_nproc <- getOption("default.nproc.blas")
  options(bigstatsr.check.parallel.blas = FALSE, default.nproc.blas = NULL)
  on.exit({
    options(bigstatsr.check.parallel.blas = old_check)
    options(default.nproc.blas = old_nproc)
  }, add = TRUE)
  force(expr)
}
