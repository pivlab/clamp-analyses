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

prepare_linear_cpm <- function(x, input_scale, label = "expression input") {
  input_scale <- tolower(as.character(input_scale))
  if (!input_scale %in% c("counts", "cpm")) {
    stop(label, ": --input-scale must be either 'counts' or 'cpm'")
  }
  if (!all(is.finite(x)) || any(x < 0)) {
    stop(label, ": input contains non-finite or negative values")
  }

  library_sizes <- colSums(x)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop(label, ": every sample must have a positive finite library size")
  }
  if (input_scale == "counts") return(CLAMP::cpmCLAMP(x))

  if (any(abs(library_sizes - 1e6) > 1e4)) {
    stop(label, ": input declared as CPM but sample totals are not approximately 1e6")
  }
  x
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

filter_analysis_cell_types <- function(x, dataset, excluded_targets, strict = TRUE) {
  excluded <- as.character(excluded_targets[[dataset]] %||% character())
  missing <- setdiff(excluded, colnames(x))
  if (strict && length(missing)) {
    stop(dataset, ": configured cell-type exclusions absent from truth: ",
         paste(missing, collapse = ", "))
  }
  x[, setdiff(colnames(x), excluded), drop = FALSE]
}

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

read_B <- function(path) t(read_csv_matrix(path))

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

    ct_col         <- cc[, ct]
    ct_col_other   <- ct_col[names(ct_col) != lv]
    r_next_best_ct <- if (length(ct_col_other) > 0) max(ct_col_other, na.rm = TRUE) else NA_real_

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
