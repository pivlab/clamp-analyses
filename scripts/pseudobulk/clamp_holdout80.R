#!/usr/bin/env Rscript
# Fit one deterministic training-only CLAMPfull model and project held-out samples.
suppressPackageStartupMessages({
  library(CLAMP)
  library(data.table)
  library(rsvd)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()

dataset <- required_arg(args, "dataset")
counts <- read_csv_matrix(required_arg(args, "bulk"))
truth <- read_csv_matrix(required_arg(args, "truth"))
info <- fread(required_arg(args, "patient_info"), data.table = FALSE)
seed <- as.integer(args$seed %||% 42L)
max_iter <- as.integer(args$max_iter %||% 500L)
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
out <- required_arg(args, "out_dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- sort(Reduce(intersect, list(colnames(counts), rownames(truth), as.character(info$sample))))
if (length(samples) < 5L) stop(dataset, ": too few aligned pseudobulk samples")

# Group by donor so a donor's samples don't split across train/test
group <- setNames(samples, samples)
split_unit <- "sample"
if ("donor_id" %in% names(info)) {
  donors <- setNames(as.character(info$donor_id), as.character(info$sample))
  usable <- !is.na(donors[samples]) & nzchar(donors[samples])
  group[usable] <- donors[samples][usable]
  split_unit <- "donor_id"
}

# 80/20 split at the group level
units <- sort(unique(unname(group)))
set.seed(seed)
n_train <- min(length(units) - 1L, max(1L, floor(0.8 * length(units))))
train_units <- sample(units, n_train, replace = FALSE)
train_samples <- sort(samples[group[samples] %in% train_units])
test_samples <- sort(setdiff(samples, train_samples))
if (length(test_samples) < 2L) stop(dataset, ": held-out partition has fewer than two samples")

# bulk_expr.csv is already CPM; subset it without a second normalization.
# Preprocess train only; reuse its row_stats to z-score test (no leakage).
train_cpm <- counts[, train_samples, drop = FALSE]
prep <- CLAMP::preprocessCLAMP(train_cpm, mean_cutoff = mean_cutoff, var_cutoff = var_cutoff)
train_norm <- CLAMP::zscoreCLAMP(prep$Y_filtered, prep$rowStats)
test_cpm <- counts[, test_samples, drop = FALSE]
genes <- intersect(rownames(train_norm), rownames(test_cpm))
train_norm <- train_norm[genes, , drop = FALSE]
row_stats <- prep$rowStats[genes, , drop = FALSE]
test_norm <- CLAMP::zscoreCLAMP(test_cpm[genes, , drop = FALSE], row_stats)

gmt <- read_gmt_file(required_arg(args, "gmt"))
pathways <- CLAMP::gmtListToSparseMat(list(BP = gmt))
set.seed(seed)

# SVD (needed by CLAMPbase and CLAMPfull). Sized from the training partition
# only, so it never leaks the held-out samples.
svd_k <- max(2L, floor((min(nrow(train_norm), ncol(train_norm)) - 1) / 4))
svdres <- rsvd::rsvd(train_norm, k = svd_k)

# Rank selection, mirroring preprocess.R's elbow-based selection but scoped
# to the training partition only (no leakage from held-out samples).
fold_k <- CLAMP::num.pc(data = train_norm, method = "elbow") * 2
fold_k <- max(as.integer(fold_k), 2L)
fold_k <- min(fold_k, svd_k)  # cannot use more singular vectors than computed

# CLAMPbase
base <- CLAMP::CLAMPbase(Y = train_norm, svdres = svdres, clamp_k = fold_k, trace = FALSE)

# Match BP prior to dataset genes
matched <- CLAMP::getMatchedPathwayMat(pathways, rownames(train_norm))

# CLAMPfull
model <- CLAMP::CLAMPfull(
  Y                 = train_norm,
  svdres            = svdres,
  priorMat          = matched,
  clamp.base.result = base,
  use_cpp           = TRUE,
  trace             = FALSE,
  max.iter          = max_iter,
  clamp_k           = fold_k
)
model$Z <- as.data.frame(model$Z); rownames(model$Z) <- rownames(train_norm)
model$B <- as.data.frame(model$B); colnames(model$B) <- colnames(train_norm)
model$Z <- as.matrix(model$Z[genes, , drop = FALSE])

# Project held-out samples onto the LVs learned from training only
test_B <- CLAMP::projectCLAMP(model, test_norm)

membership <- data.frame(
  dataset = dataset, sample = samples, group_id = unname(group[samples]), split_unit = split_unit,
  split = ifelse(samples %in% train_samples, "train", "test"), seed = seed
)
write.csv(model$B, file.path(out, "train_B.csv"))
write.csv(model$Z, file.path(out, "train_Z.csv"))
write.csv(test_B, file.path(out, "test_B.csv"))
write.csv(row_stats, file.path(out, "row_stats.csv"))
write.csv(truth[train_samples, , drop = FALSE], file.path(out, "train_truth.csv"))
write.csv(truth[test_samples, , drop = FALSE], file.path(out, "test_truth.csv"))
fwrite(membership, file.path(out, "split_membership.csv"))
fwrite(data.frame(
  dataset = dataset, seed = seed, split_unit = split_unit, n_train = length(train_samples),
  n_test = length(test_samples), k = fold_k, svd_k = svd_k, n_genes = length(genes)
), file.path(out, "summary.csv"))
saveRDS(model, file.path(out, "CLAMPfull_train80.rds"))
message(dataset, ": holdout model complete (", length(train_samples), "/", length(test_samples), ")")
