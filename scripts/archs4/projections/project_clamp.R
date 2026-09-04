#!/usr/bin/env Rscript
# Project one or more datasets into a CLAMP model with CLAMP::projectCLAMP.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(CLAMP)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
model_path <- required_arg(args, "model")
targets_path <- required_arg(args, "targets")
out_dir <- required_arg(args, "out_dir")
scale <- as.numeric(args$scale %||% 1)
null_seeds <- if (is.null(args$null_seeds)) integer(0) else
  as.integer(strsplit(as.character(args$null_seeds), ",")[[1]])
null_base_seed <- as.integer(args$null_base_seed %||% 0L)
null_pattern <- args$null_out_pattern %||% ""
if (length(null_seeds) && !nzchar(null_pattern)) stop("--null-seeds needs --null-out-pattern")

started <- Sys.time()
targets <- fread(targets_path)
for (col in c("id", "input_path", "input_kind", "out_path")) {
  if (!col %in% names(targets)) stop("targets TSV is missing column: ", col)
}

tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

message("Loading CLAMP model: ", model_path)
model <- read_projection_model(model_path)
message("Model: ", nrow(model$Z), " genes x ", ncol(model$Z), " LVs | L2 = ", signif(model$L2, 6))

load_newdata <- function(row) {
  if (identical(row$input_kind, "fbm")) {
    fbm <- readRDS(row$input_path)
    if (inherits(fbm, "list") && !is.null(fbm$fbm)) fbm <- fbm$fbm
    fbm$is_read_only <- TRUE
    genes <- as.character(readRDS(row$genes_rds))
    samples <- as.character(readRDS(row$samples_rds))
    if (nrow(fbm) != length(genes) || ncol(fbm) != length(samples)) {
      stop("FBM dimensions do not match genes/samples for target ", row$id)
    }
    return(list(x = fbm, genes = genes, samples = samples, is_fbm = TRUE))
  }
  m <- read_matrix_csv(row$input_path)
  list(x = m, genes = rownames(m), samples = colnames(m), is_fbm = FALSE)
}

align_to_model <- function(mdl, nd) {
  common <- intersect(rownames(mdl$Z), nd$genes)
  if (!length(common)) stop("No common genes between model Z and newdata")
  mdl$Z <- mdl$Z[common, , drop = FALSE]
  idx <- match(common, nd$genes)
  x <- if (nd$is_fbm) bigstatsr::big_copy(nd$x, ind.row = idx) else nd$x[idx, , drop = FALSE]
  list(mdl = mdl, newdata = x, samples = nd$samples, common = common)
}

project_one <- function(mdl, newdata, samples, out_path, label) {
  stopifnot(nrow(mdl$Z) == nrow(newdata))
  B <- CLAMP::projectCLAMP(CLAMPres = mdl, newdata = newdata, scale = scale)
  rownames(B) <- colnames(mdl$Z)
  colnames(B) <- samples
  ensure_parent(out_path)
  write_matrix_csv(B, out_path)
  message("  ", label, " -> ", out_path, "  (", nrow(B), " LVs x ", ncol(B), " samples)")
  data.table(n_lvs = nrow(B), n_samples = ncol(B))
}

summaries <- list()
for (i in seq_len(nrow(targets))) {
  row <- as.list(targets[i])
  message("[", i, "/", nrow(targets), "] ", row$id)
  nd <- load_newdata(row)
  n_input_genes <- length(nd$genes)
  al <- align_to_model(model, nd)
  model_aligned <- al$mdl
  newdata <- al$newdata
  samples <- al$samples
  common <- length(al$common)
  message("  aligned on ", common, " common genes (",
          sprintf("%.1f%%", 100 * common / nrow(model$Z)), " of model genes)")
  st <- project_one(model_aligned, newdata, samples, row$out_path, "real")
  summaries[[length(summaries) + 1L]] <- data.table(
    id = row$id, kind = "real", seed = NA_integer_, out_path = row$out_path,
    n_input_genes = n_input_genes, n_common_genes = common,
    n_lvs = st$n_lvs, n_samples = st$n_samples)

  for (s in null_seeds) {
    null_model <- permute_model_z(model_aligned, null_base_seed + s)
    out_path <- gsub("{id}", row$id, gsub("{seed}", s, null_pattern, fixed = TRUE), fixed = TRUE)
    st <- project_one(null_model, newdata, samples, out_path, paste0("null seed ", s))
    summaries[[length(summaries) + 1L]] <- data.table(
      id = row$id, kind = "null", seed = s, out_path = out_path,
      n_input_genes = n_input_genes, n_common_genes = common,
      n_lvs = st$n_lvs, n_samples = st$n_samples)
    rm(null_model); gc(verbose = FALSE)
  }
  rm(newdata, nd, model_aligned, al); gc(verbose = FALSE)
}

summary_dt <- rbindlist(summaries)
fwrite(summary_dt, file.path(tmp_dir, "targets_summary.csv"))

write_manifest(
  tmp_dir, method = "CLAMP::projectCLAMP",
  inputs = list(model = model_path, targets = targets_path),
  parameters = list(scale = scale, align = "manual_intersect",
                    null_method = "permute_z_rows",
                    null_base_seed = null_base_seed,
                    null_seeds = null_seeds),
  extra = list(model_n_genes = nrow(model$Z), model_n_lvs = ncol(model$Z),
               model_l2 = model$L2, n_targets = nrow(targets),
               n_outputs = nrow(summary_dt)),
  started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
