#!/usr/bin/env Rscript
# Create one deterministic, group-preserving fold assignment per pseudobulk cohort.
suppressPackageStartupMessages(library(data.table))

script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()

dataset <- required_arg(args, "dataset")
bulk_header <- readLines(required_arg(args, "bulk"), n = 1L, warn = FALSE)
count_samples <- scan(
  text = bulk_header,
  what = character(),
  sep = ",",
  quote = '"',
  quiet = TRUE,
  blank.lines.skip = FALSE
)
truth <- read_csv_matrix(required_arg(args, "truth"))
info <- fread(required_arg(args, "patient_info"), data.table = FALSE)
out <- required_arg(args, "output")
summary_out <- required_arg(args, "summary")
n_folds <- as.integer(args$n_folds %||% 5L)
seed <- as.integer(args$seed %||% 42L)
group_col <- as.character(args$group_col %||% "sample")

if (!"sample" %in% names(info)) stop(dataset, ": patient_info is missing the sample column")
if (!group_col %in% names(info)) stop(dataset, ": patient_info is missing group column ", group_col)
if (n_folds < 2L) stop("n_folds must be at least 2")
if (anyDuplicated(info$sample)) stop(dataset, ": patient_info contains duplicate sample rows")

samples <- sort(Reduce(intersect, list(count_samples, rownames(truth), as.character(info$sample))))
if (length(samples) < n_folds) stop(dataset, ": fewer aligned samples than folds")

sample_to_group <- setNames(as.character(info[[group_col]]), as.character(info$sample))
groups <- unname(sample_to_group[samples])
if (anyNA(groups) || any(!nzchar(groups))) stop(dataset, ": missing group identifiers among aligned samples")
if (length(unique(groups)) < n_folds) stop(dataset, ": fewer unique groups than folds")

# Seeded greedy bin-packing keeps whole groups together while balancing the
# number of samples per fold. Randomized tie ordering makes the seed effective.
group_sizes <- table(groups)
set.seed(seed)
shuffled_groups <- sample(names(group_sizes), length(group_sizes), replace = FALSE)
ordered_groups <- shuffled_groups[order(-as.integer(group_sizes[shuffled_groups]))]
fold_load <- integer(n_folds)
group_fold <- setNames(integer(length(ordered_groups)), ordered_groups)
for (group_id in ordered_groups) {
  lightest <- which(fold_load == min(fold_load))
  chosen <- if (length(lightest) == 1L) lightest[[1L]] else sample(lightest, 1L)
  group_fold[[group_id]] <- chosen
  fold_load[[chosen]] <- fold_load[[chosen]] + as.integer(group_sizes[[group_id]])
}

membership <- data.frame(
  dataset = dataset,
  sample = samples,
  group_id = groups,
  group_col = group_col,
  fold = unname(group_fold[groups]),
  split_seed = seed,
  stringsAsFactors = FALSE
)
membership <- membership[order(membership$sample), , drop = FALSE]

group_fold_counts <- as.data.table(membership)[, uniqueN(fold), by = group_id]
if (any(group_fold_counts$V1 != 1L)) stop(dataset, ": a group crosses folds")
if (!identical(sort(unique(membership$fold)), seq_len(n_folds))) stop(dataset, ": at least one fold is empty")

fold_summary <- as.data.table(membership)[, .(
  n_samples = .N,
  n_groups = uniqueN(group_id)
), by = .(dataset, fold, group_col, split_seed)]
fold_summary[, `:=`(n_folds = n_folds, n_total_samples = length(samples), n_total_groups = length(unique(groups)))]

ensure_parent(out)
ensure_parent(summary_out)
fwrite(membership, out)
fwrite(fold_summary, summary_out)
message(dataset, ": wrote ", n_folds, " grouped folds for ", length(samples), " samples")
