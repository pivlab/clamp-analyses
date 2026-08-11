# Shared pseudobulk helpers used by multiple production scripts and reports.
# Batch scripts must never fall back to R's implicit `Rplots.pdf` device.
# Explicit devices (for example ggsave/cairo_pdf) are unaffected.
options(device = function(...) grDevices::pdf(file = NULL))

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

ensure_parent <- function(path) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

read_csv_matrix <- function(path) {
  dt <- data.table::fread(path)
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt[[1]]
  mat
}

write_csv_matrix <- function(x, path) {
  ensure_parent(path)
  write.csv(as.data.frame(x), path)
}

read_gmt_file <- function(path, prefix = "BP_") {
  lines <- readLines(path, warn = FALSE)
  pieces <- strsplit(lines, "\t", fixed = TRUE)
  result <- lapply(pieces, function(x) unique(x[-c(1L, 2L)]))
  names(result) <- paste0(prefix, vapply(pieces, `[[`, character(1), 1L))
  result
}

read_k <- function(k_path) as.integer(data.table::fread(k_path)$k[[1L]])

read_norm_and_k <- function(norm_path, k_path) {
  list(norm = read_csv_matrix(norm_path), k = read_k(k_path))
}

# Drop cell types with mean prevalence below `min_prevalence` in small
# studies (nrow <= min_nrow_for_filter), where a handful of samples makes a
# rare cell type's fraction estimate too noisy to trust. Works on both a
# matrix (e.g. from read_csv_matrix()) and a data.frame (e.g. from
# read.csv(..., row.names = 1)) since colMeans()/`[` behave the same on both.
filter_rare_cell_types <- function(x, min_nrow_for_filter = 100, min_prevalence = 0.005) {
  if (nrow(x) <= min_nrow_for_filter) {
    keep <- colMeans(x, na.rm = TRUE) >= min_prevalence
    x <- x[, keep, drop = FALSE]
  }
  x
}

# Greedy one-to-one assignment used by the pseudobulk benchmark.
# Rows and columns can each be selected at most once.
oneToOneMask <- function(cc) {
  if (!is.matrix(cc)) cc <- as.matrix(cc)
  mask <- matrix(NA_real_, nrow(cc), ncol(cc), dimnames = dimnames(cc))
  cc_work <- cc
  cc_work[!is.finite(cc_work)] <- NA_real_

  while (any(is.finite(cc_work))) {
    best <- which(cc_work == max(cc_work, na.rm = TRUE), arr.ind = TRUE)
    if (!nrow(best)) break
    row_index <- best[1L, 1L]
    column_index <- best[1L, 2L]
    mask[row_index, column_index] <- cc_work[row_index, column_index]
    cc_work[row_index, ] <- NA_real_
    cc_work[, column_index] <- NA_real_
  }

  mask
}

# Read a factorization's B matrix and transpose it to samples x LVs, the
# orientation score_per_ct/assign_with_margins expect. B.csv is stored on
# disk as LVs x samples.
# path: path to a B.csv written by a CLAMP/PLIER/NMF/... factorization run.
# Returns: numeric matrix, samples x LVs.
read_B <- function(path) t(read_csv_matrix(path))

# Score a factorization's LVs against ground-truth cell-type fractions, using
# oneToOneMask so each LV can be claimed by at most one cell type -- a method
# can't get credit for the same LV matching multiple cell types. When there
# are more cell types than LVs, some cell types will have no LV left to
# claim; that's a real limitation of the factorization for that dataset, not
# an artifact to mask.
# fac: samples x LVs numeric matrix (e.g. from read_B()).
# truth: samples x cell-types numeric matrix of ground-truth fractions.
# ct_names: character vector of cell-type names to report a score for.
# Returns: named numeric vector (one entry per ct_names) with the Pearson
# correlation of the assigned LV for each cell type, or NA where fac and
# truth share no samples, fac has no non-constant columns, or the cell type
# had no LV left to claim.
score_per_ct <- function(fac, truth, ct_names) {
  shared <- intersect(rownames(fac), rownames(truth))
  if (length(shared) == 0) return(setNames(rep(NA_real_, length(ct_names)), ct_names))
  fac_s   <- fac[shared, , drop = FALSE]
  truth_s <- truth[shared, , drop = FALSE]

  fac_s <- fac_s[, apply(fac_s, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(fac_s) == 0) return(setNames(rep(NA_real_, length(ct_names)), ct_names))

  cc  <- cor(fac_s, truth_s, method = "pearson", use = "pairwise.complete.obs")
  cc  <- oneToOneMask(cc)
  v   <- apply(cc, 2, function(x) if (all(is.na(x))) NA_real_ else x[which(!is.na(x))[1]])
  v[!is.finite(v)] <- NA_real_
  setNames(v, ct_names)
}

# Same one-to-one assignment as score_per_ct, but also returns the full
# unmasked correlation matrix and per-assignment specificity margins: how
# much higher the assigned pair's correlation is than the next-best
# competing LV/cell type. Used to check whether related cell types map to
# distinct LVs (e.g. by 02_disentangle.ipynb).
# fac, truth: as in score_per_ct().
# ct_names: cell types to compute an assignment and margin for.
# Returns: a list with
#   cc: long data.frame (LV, cell_type, cor) of the full unmasked correlation
#     matrix, or NULL if fac and truth share no samples or fac has no
#     non-constant columns.
#   assignment: data.frame with one row per assigned cell type -- cell_type,
#     LV, cor, r_next_best_ct/r_next_best_lv (next-best competing
#     correlation for the same cell type / LV), margin_ct/margin_lv (how
#     much the assigned pair beats that next-best competitor). NULL under
#     the same conditions as cc, or if no cell type in ct_names got an
#     assignment.
assign_with_margins <- function(fac, truth, ct_names) {
  shared <- intersect(rownames(fac), rownames(truth))
  if (length(shared) == 0) return(list(cc = NULL, assignment = NULL))
  fac_s   <- fac[shared, , drop = FALSE]
  truth_s <- truth[shared, , drop = FALSE]
  fac_s   <- fac_s[, apply(fac_s, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(fac_s) == 0) return(list(cc = NULL, assignment = NULL))

  cc   <- cor(fac_s, truth_s, method = "pearson", use = "pairwise.complete.obs")
  mask <- oneToOneMask(cc)

  assignment <- do.call(rbind, lapply(ct_names, function(ct) {
    if (!ct %in% colnames(mask)) return(NULL)
    col_vals <- mask[, ct]
    assigned <- which(!is.na(col_vals))
    if (length(assigned) == 0) return(NULL)
    lv <- rownames(mask)[assigned[1]]
    r  <- col_vals[assigned[1]]

    # Next-best LV for this cell type (excluding the assigned LV itself).
    ct_col         <- cc[, ct]
    ct_col_other   <- ct_col[names(ct_col) != lv]
    r_next_best_ct <- if (length(ct_col_other) > 0) max(ct_col_other, na.rm = TRUE) else NA_real_

    # Next-best cell type for this LV (excluding the assigned cell type itself).
    lv_row         <- cc[lv, ]
    lv_row_other   <- lv_row[names(lv_row) != ct]
    r_next_best_lv <- if (length(lv_row_other) > 0) max(lv_row_other, na.rm = TRUE) else NA_real_

    data.frame(
      cell_type      = ct,
      LV             = lv,
      cor            = r,
      r_next_best_ct = r_next_best_ct,
      r_next_best_lv = r_next_best_lv,
      margin_ct      = r - r_next_best_ct,
      margin_lv      = r - r_next_best_lv,
      stringsAsFactors = FALSE
    )
  }))

  cc_long <- as.data.frame(as.table(cc), stringsAsFactors = FALSE)
  names(cc_long) <- c("LV", "cell_type", "cor")

  list(cc = cc_long, assignment = assignment)
}
