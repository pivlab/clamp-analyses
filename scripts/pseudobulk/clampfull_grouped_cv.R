#!/usr/bin/env Rscript
# Fit CLAMPfull on the grouped training folds and project the untouched fold.
suppressPackageStartupMessages({
  library(CLAMP)
  library(data.table)
  library(rsvd)
})

script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()

dataset <- required_arg(args, "dataset")
fold <- as.integer(required_arg(args, "fold"))
counts <- read_csv_matrix(required_arg(args, "bulk"))
truth <- read_csv_matrix(required_arg(args, "truth"))
membership <- fread(required_arg(args, "membership"), data.table = FALSE)
membership$sample <- as.character(membership$sample)
membership$group_id <- as.character(membership$group_id)
base_seed <- as.integer(args$seed %||% 123L)
model_seed <- base_seed + fold
max_iter <- as.integer(args$max_iter %||% 500L)
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
out <- required_arg(args, "out_dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (!all(c("dataset", "sample", "group_id", "fold", "split_seed") %in% names(membership))) {
  stop(dataset, ": malformed fold membership file")
}
# Restrict to this dataset's fold assignment, then split into this fold's
# held-out test samples vs. every other sample for training.
membership <- membership[membership$dataset == dataset, , drop = FALSE]
split_seed <- unique(membership$split_seed)
if (length(split_seed) != 1L) stop(dataset, ": split seed is not unique")
samples <- sort(Reduce(intersect, list(colnames(counts), rownames(truth), membership$sample)))
test_samples <- sort(intersect(samples, membership$sample[membership$fold == fold]))
train_samples <- sort(setdiff(samples, test_samples))
if (length(test_samples) < 2L) stop(dataset, " fold ", fold, ": fewer than two held-out samples")
if (length(train_samples) < 2L) stop(dataset, " fold ", fold, ": fewer than two training samples")

# Guard against a group (e.g. patient) appearing on both sides of the split
train_groups <- unique(membership$group_id[membership$sample %in% train_samples])
test_groups <- unique(membership$group_id[membership$sample %in% test_samples])
if (length(intersect(train_groups, test_groups)) > 0L) stop(dataset, " fold ", fold, ": group leakage")

# bulk_expr.csv is already CPM; subset it without a second normalization.
# Every data-derived preprocessing decision is learned from the training fold.
train_cpm <- counts[, train_samples, drop = FALSE]
prep <- CLAMP::preprocessCLAMP(train_cpm, mean_cutoff = mean_cutoff, var_cutoff = var_cutoff)
train_norm <- CLAMP::zscoreCLAMP(prep$Y_filtered, prep$rowStats)
test_cpm <- counts[, test_samples, drop = FALSE]
genes <- intersect(rownames(train_norm), rownames(test_cpm))
if (length(genes) < 2L) stop(dataset, " fold ", fold, ": insufficient retained genes")
train_norm <- train_norm[genes, , drop = FALSE]
row_stats <- prep$rowStats[genes, , drop = FALSE]
test_norm <- CLAMP::zscoreCLAMP(test_cpm[genes, , drop = FALSE], row_stats)

gmt <- read_gmt_file(required_arg(args, "gmt"))
pathways <- CLAMP::gmtListToSparseMat(list(BP = gmt))
set.seed(model_seed)

# SVD (needed by CLAMPbase and CLAMPfull). Sized from this fold's own
# training data, not the full dataset, so it never leaks held-out samples.
svd_k <- max(2L, floor((min(nrow(train_norm), ncol(train_norm)) - 1) / 4))
svdres <- rsvd::rsvd(train_norm, k = svd_k)

# Rank selection, mirroring preprocess.R's elbow-based selection but scoped
# to this fold's training data only (no leakage from held-out samples).
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

# Project the held-out fold onto the LVs learned from training only
test_B <- CLAMP::projectCLAMP(model, test_norm)

write.csv(model$B, file.path(out, "train_B.csv"))
write.csv(model$Z, file.path(out, "train_Z.csv"))
write.csv(test_B, file.path(out, "test_B.csv"))
write.csv(row_stats, file.path(out, "row_stats.csv"))
write.csv(truth[train_samples, , drop = FALSE], file.path(out, "train_truth.csv"))
write.csv(truth[test_samples, , drop = FALSE], file.path(out, "test_truth.csv"))
fwrite(data.frame(
  dataset = dataset,
  method = "CLAMPfull",
  fold = fold,
  split_seed = split_seed,
  model_seed = model_seed,
  n_train = length(train_samples),
  n_test = length(test_samples),
  n_train_groups = length(train_groups),
  n_test_groups = length(test_groups),
  k = fold_k,
  svd_k = svd_k,
  n_genes = length(genes)
), file.path(out, "summary.csv"))
saveRDS(model, file.path(out, "CLAMPfull.rds"))
message(dataset, " fold ", fold, ": CLAMPfull complete (", length(train_samples), "/", length(test_samples), ")")
