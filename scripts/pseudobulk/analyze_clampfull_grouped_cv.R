#!/usr/bin/env Rscript
# Compute CLAMPfull grouped-CV calibrations, OOF metrics, bootstrap summaries, and LV stability.
# Plotting is intentionally kept in the holdout notebook.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(yaml)
})
source(here('scripts', 'pseudobulk', 'common.R'))

cfg <- yaml::read_yaml(here('workflow', 'config', 'pseudobulk.yaml'))
analysis_cfg <- yaml::read_yaml(here('workflow', 'config', 'cell_type_analysis.yaml'))
DATASETS <- names(cfg$datasets)
N_FOLDS <- as.integer(cfg$grouped_cv$n_folds)
MIN_TRAIN_LV_COR <- as.numeric(cfg$grouped_cv$min_train_lv_cor)
BOOTSTRAP_SEED <- as.integer(cfg$grouped_cv$bootstrap_seed)
BOOTSTRAP_REPLICATES <- as.integer(cfg$grouped_cv$bootstrap_replicates)
PROD <- here(cfg$paths$production)
OUT <- here(cfg$paths$production, 'grouped_cv_analysis')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
METHOD <- 'CLAMPfull'

read_matrix <- function(path) {
  x <- suppressWarnings(fread(path))
  m <- as.matrix(x[, -1, with = FALSE])
  storage.mode(m) <- 'numeric'
  rownames(m) <- x[[1L]]
  m
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L || var(x[keep]) <= 0 || var(y[keep]) <= 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = 'pearson'))
}

numeric_lv_order <- function(x) {
  value <- suppressWarnings(as.integer(sub('^.*?([0-9]+)$', '\\1', x)))
  value[is.na(value)] <- .Machine$integer.max
  order(value, x)
}

match_reference_lvs <- function(fold_z, reference_z) {
  fold_lvs <- colnames(fold_z)
  reference_lvs <- setNames(rep(NA_character_, length(fold_lvs)), fold_lvs)
  match_abs_cor <- setNames(rep(NA_real_, length(fold_lvs)), fold_lvs)
  genes <- intersect(rownames(fold_z), rownames(reference_z))
  if (length(genes) < 3L) {
    return(list(reference_lv = reference_lvs, abs_cor = match_abs_cor))
  }

  similarity <- abs(cor(
    fold_z[genes, , drop = FALSE],
    reference_z[genes, , drop = FALSE],
    use = 'pairwise.complete.obs', method = 'pearson'
  ))
  matched <- oneToOneMask(similarity)
  for (row_index in seq_len(nrow(matched))) {
    column_index <- which(is.finite(matched[row_index, ]))
    if (!length(column_index)) next
    fold_lv <- rownames(matched)[[row_index]]
    reference_lvs[[fold_lv]] <- colnames(matched)[[column_index[[1L]]]]
    match_abs_cor[[fold_lv]] <- matched[row_index, column_index[[1L]]]
  }
  list(reference_lv = reference_lvs, abs_cor = match_abs_cor)
}


calibration_rows <- list()
prediction_rows <- list()
membership_rows <- list()
row_index <- 0L

for (dataset in DATASETS) {
  cv_root <- file.path(PROD, dataset, 'grouped_cv')
  membership <- fread(file.path(cv_root, 'fold_membership.csv'))
  membership[, `:=`(sample = as.character(sample), group_id = as.character(group_id))]
  required_membership <- c('dataset', 'sample', 'group_id', 'fold', 'split_seed')
  if (!all(required_membership %in% names(membership))) stop(dataset, ': malformed membership table')
  if (membership[, uniqueN(fold), by = group_id][V1 != 1L, .N] > 0L) stop(dataset, ': group leakage')
  if (membership[, uniqueN(sample)] != nrow(membership)) stop(dataset, ': duplicate membership rows')
  membership_rows[[dataset]] <- membership

  # Align selected LVs across training folds using fold 1 as the reference.
  # This measures CV selection stability without depending on a full-data fit.
  reference_z <- read_matrix(file.path(cv_root, 'fold1', METHOD, 'train_Z.csv'))

  for (fold in seq_len(N_FOLDS)) {
    fold_root <- file.path(cv_root, paste0('fold', fold), METHOD)
    train_b <- t(read_matrix(file.path(fold_root, 'train_B.csv')))
    test_b <- t(read_matrix(file.path(fold_root, 'test_B.csv')))
    fold_z <- read_matrix(file.path(fold_root, 'train_Z.csv'))
    train_truth <- read_matrix(file.path(fold_root, 'train_truth.csv'))
    test_truth <- read_matrix(file.path(fold_root, 'test_truth.csv'))
    train_truth <- filter_analysis_cell_types(
      train_truth, dataset, analysis_cfg$cell_type_analysis$excluded_targets,
      strict = FALSE
    )
    test_truth <- filter_analysis_cell_types(
      test_truth, dataset, analysis_cfg$cell_type_analysis$excluded_targets,
      strict = FALSE
    )

    train_ids <- sort(Reduce(intersect, list(rownames(train_b), rownames(train_truth))))
    test_ids <- sort(Reduce(intersect, list(rownames(test_b), rownames(test_truth))))
    fold_id <- fold
    expected_test <- sort(membership[fold == fold_id, sample])
    if (!identical(test_ids, expected_test)) stop(dataset, ' fold ', fold, ': held-out samples disagree with membership')
    if (length(intersect(train_ids, test_ids))) stop(dataset, ' fold ', fold, ': train/test overlap')
    if (!identical(sort(colnames(train_truth)), sort(colnames(test_truth)))) {
      stop(dataset, ' fold ', fold, ': train/test truth columns differ')
    }

    train_b <- train_b[train_ids, , drop = FALSE]
    test_b <- test_b[test_ids, , drop = FALSE]
    train_truth <- train_truth[train_ids, , drop = FALSE]
    test_truth <- test_truth[test_ids, , drop = FALSE]
    correlations <- cor(train_b, train_truth, use = 'pairwise.complete.obs', method = 'pearson')
    assignment_mask <- oneToOneMask(correlations)
    assignments <- setNames(vapply(seq_len(ncol(assignment_mask)), function(column_index) {
      row_index <- which(is.finite(assignment_mask[, column_index]))
      if (length(row_index)) rownames(assignment_mask)[[row_index[[1L]]]] else NA_character_
    }, character(1)), colnames(assignment_mask))
    reference_matches <- match_reference_lvs(fold_z, reference_z)

    for (cell_type in colnames(train_truth)) {
      row_index <- row_index + 1L
      y_train <- train_truth[, cell_type]
      y_test <- test_truth[, cell_type]
      status <- 'ok'
      selected_lv <- NA_character_
      train_lv_cor <- NA_real_
      test_lv_cor <- NA_real_
      alpha <- mean(y_train, na.rm = TRUE)
      beta <- 0
      train_prediction <- rep(alpha, length(y_train))
      test_prediction <- rep(alpha, length(y_test))
      test_lv_score <- rep(NA_real_, length(y_test))

      if (sum(is.finite(y_train)) < 2L || var(y_train, na.rm = TRUE) <= 0) {
        status <- 'constant_training_truth'
      } else {
        selected_lv <- assignments[[cell_type]]
        if (is.na(selected_lv)) {
          status <- 'no_one_to_one_lv'
        } else {
          x_train <- train_b[, selected_lv]
          x_test <- test_b[, selected_lv]
          train_lv_cor <- safe_cor(x_train, y_train)
          test_lv_cor <- safe_cor(x_test, y_test)
          fit_rows <- is.finite(x_train) & is.finite(y_train)
          if (sum(fit_rows) < 2L || var(x_train[fit_rows]) <= 0) {
            status <- 'constant_selected_lv'
          } else {
            calibration <- lm(y_train[fit_rows] ~ x_train[fit_rows])
            coefficients <- coef(calibration)
            if (length(coefficients) < 2L || any(!is.finite(coefficients))) {
              status <- 'calibration_failure'
            } else {
              alpha <- unname(coefficients[[1L]])
              beta <- unname(coefficients[[2L]])
              train_prediction <- alpha + beta * x_train
              test_prediction <- alpha + beta * x_test
              test_lv_score <- x_test
            }
          }
        }
      }

      reference_lv <- if (is.na(selected_lv)) NA_character_ else reference_matches$reference_lv[[selected_lv]]
      reference_match_abs_cor <- if (is.na(selected_lv)) NA_real_ else reference_matches$abs_cor[[selected_lv]]
      eligible_train_signal <- status == 'ok' && is.finite(train_lv_cor) &&
        train_lv_cor >= MIN_TRAIN_LV_COR
      calibration_rows[[row_index]] <- data.table(
        dataset = dataset, method = METHOD, fold = fold, cell_type = cell_type,
        selected_lv = selected_lv, train_lv_cor = train_lv_cor,
        test_lv_cor = test_lv_cor,
        train_abs_cor = abs(train_lv_cor),
        train_calibration_cor = safe_cor(train_prediction, y_train),
        alpha = alpha, beta = beta, n_train = length(train_ids), n_test = length(test_ids),
        reference_lv = reference_lv,
        reference_match_abs_cor = reference_match_abs_cor,
        eligible_train_signal = eligible_train_signal,
        min_train_lv_cor = MIN_TRAIN_LV_COR, status = status
      )
      prediction_rows[[row_index]] <- data.table(
        dataset = dataset, method = METHOD, fold = fold, sample = test_ids,
        group_id = membership[match(test_ids, sample), group_id],
        cell_type = cell_type, selected_lv = selected_lv, lv_score = test_lv_score,
        observed = y_test, predicted = test_prediction, train_lv_cor = train_lv_cor,
        eligible_train_signal = eligible_train_signal,
        min_train_lv_cor = MIN_TRAIN_LV_COR, status = status
      )
    }
  }
}

fold_membership <- rbindlist(membership_rows, use.names = TRUE, fill = TRUE)
fold_calibrations <- rbindlist(calibration_rows, use.names = TRUE, fill = TRUE)
oof_predictions <- rbindlist(prediction_rows, use.names = TRUE, fill = TRUE)
cat('Fold calibrations:', nrow(fold_calibrations), '\n')
cat('Out-of-fold predictions:', nrow(oof_predictions), '\n')


prediction_key <- c('dataset', 'method', 'sample', 'cell_type')
if (oof_predictions[, anyDuplicated(.SD), .SDcols = prediction_key] > 0L) stop('duplicate out-of-fold predictions')
coverage <- oof_predictions[, .(n_predictions = .N), by = .(dataset, method, sample)]
expected_coverage <- fold_calibrations[, .(n_cell_types = uniqueN(cell_type)), by = dataset]
coverage <- merge(coverage, expected_coverage, by = 'dataset', all.x = TRUE)
if (coverage[n_predictions != n_cell_types, .N] > 0L) stop('incomplete cell-type prediction coverage')
predicted_sample_keys <- unique(oof_predictions[, paste(dataset, sample, sep = '\r')])
member_sample_keys <- unique(fold_membership[, paste(dataset, sample, sep = '\r')])
if (!setequal(predicted_sample_keys, member_sample_keys)) stop('not every member sample has an OOF prediction')

oof_metrics <- oof_predictions[, {
  keep <- is.finite(observed) & is.finite(predicted)
  observed_valid <- observed[keep]
  predicted_valid <- predicted[keep]
  sst <- sum((observed_valid - mean(observed_valid))^2)
  list(
    n_oof = sum(keep),
    pearson_r = safe_cor(predicted_valid, observed_valid),
    predictive_r2 = if (length(observed_valid) && sst > 0) 1 - sum((observed_valid - predicted_valid)^2) / sst else NA_real_,
    mae = if (length(observed_valid)) mean(abs(observed_valid - predicted_valid)) else NA_real_,
    n_outside_unit_interval = sum(predicted_valid < 0 | predicted_valid > 1)
  )
}, by = .(dataset, method, cell_type)]

training_summary <- fold_calibrations[, .(
  median_train_lv_cor = if (any(is.finite(train_lv_cor))) median(train_lv_cor, na.rm = TRUE) else NA_real_,
  median_train_abs_cor = if (any(is.finite(train_abs_cor))) median(train_abs_cor, na.rm = TRUE) else NA_real_,
  mean_train_abs_cor = if (any(is.finite(train_abs_cor))) mean(train_abs_cor, na.rm = TRUE) else NA_real_,
  n_successful_folds = sum(status == 'ok'), n_folds = .N
), by = .(dataset, method, cell_type)]
oof_metrics <- merge(oof_metrics, training_summary, by = c('dataset', 'method', 'cell_type'), all.x = TRUE)

test_lv_summary <- fold_calibrations[, {
  valid <- is.finite(test_lv_cor)
  correlations <- test_lv_cor[valid]
  weights <- pmax(n_test[valid] - 3, 1)
  bounded <- pmin(pmax(correlations, -1 + 1e-12), 1 - 1e-12)
  list(
    median_test_lv_cor = if (length(correlations)) median(correlations) else NA_real_,
    mean_test_lv_cor = if (length(correlations)) mean(correlations) else NA_real_,
    fisher_test_lv_cor = if (length(correlations)) tanh(weighted.mean(atanh(bounded), weights)) else NA_real_,
    n_test_cor_folds = length(correlations)
  )
}, by = .(dataset, method, cell_type)]
oof_metrics <- merge(oof_metrics, test_lv_summary, by = c('dataset', 'method', 'cell_type'), all.x = TRUE)

thresholded_metrics <- oof_predictions[, {
  finite_prediction <- is.finite(observed) & is.finite(predicted)
  n_total <- sum(is.finite(observed) & is.finite(predicted))
  folds_meeting_threshold <- uniqueN(fold[eligible_train_signal])
  total_folds <- uniqueN(fold)
  valid_all_folds <- total_folds == N_FOLDS && folds_meeting_threshold == total_folds
  valid_prediction <- valid_all_folds & finite_prediction
  observed_valid <- observed[valid_prediction]
  predicted_valid <- predicted[valid_prediction]
  valid_sst <- sum((observed_valid - mean(observed_valid))^2)
  n_valid <- sum(valid_prediction)
  list(
    min_train_lv_cor = MIN_TRAIN_LV_COR,
    minimum_observed_train_lv_cor = if (any(is.finite(train_lv_cor))) {
      min(train_lv_cor, na.rm = TRUE)
    } else NA_real_,
    n_eligible_folds = folds_meeting_threshold,
    n_total_folds = total_folds,
    valid_all_folds = valid_all_folds,
    n_eligible_test = n_valid,
    n_total_test = n_total,
    prediction_coverage = if (n_total > 0L) n_valid / n_total else NA_real_,
    coverage_class = fifelse(valid_all_folds, 'valid', 'invalid'),
    pearson_r = safe_cor(predicted_valid, observed_valid),
    predictive_r2 = if (n_valid > 0L && valid_sst > 0) {
      1 - sum((observed_valid - predicted_valid)^2) / valid_sst
    } else NA_real_,
    mae = if (n_valid > 0L) mean(abs(observed_valid - predicted_valid)) else NA_real_,
    n_outside_unit_interval = sum(
      predicted_valid < 0 | predicted_valid > 1, na.rm = TRUE
    )
  )
}, by = .(dataset, method, cell_type)]

valid_thresholded <- thresholded_metrics[valid_all_folds == TRUE]
thresholded_summary <- thresholded_metrics[, .(
  min_train_lv_cor = MIN_TRAIN_LV_COR,
  n_dataset_cell_types = .N,
  n_valid = sum(valid_all_folds),
  n_invalid = sum(!valid_all_folds),
  n_eligible_test = sum(n_eligible_test),
  n_total_test = sum(n_total_test)
), by = method]
thresholded_summary[, `:=`(
  mean_pearson_r = mean(valid_thresholded$pearson_r, na.rm = TRUE),
  mean_predictive_r2 = mean(valid_thresholded$predictive_r2, na.rm = TRUE),
  mean_mae = mean(valid_thresholded$mae, na.rm = TRUE)
)]

lv_selection_stability <- fold_calibrations[, {
  valid_reference <- reference_lv[status == 'ok' & !is.na(reference_lv)]
  valid_similarity <- reference_match_abs_cor[status == 'ok' & is.finite(reference_match_abs_cor)]
  if (!length(valid_reference)) {
    list(most_frequent_aligned_lv = NA_character_, fraction_folds_most_frequent_lv = NA_real_,
         n_distinct_aligned_lvs_across_folds = 0L, median_alignment_abs_cor = NA_real_,
         n_folds_with_aligned_lv = 0L, n_folds = .N)
  } else {
    counts <- sort(table(valid_reference), decreasing = TRUE)
    tied <- names(counts)[counts == max(counts)]
    modal <- tied[numeric_lv_order(tied)][1L]
    list(most_frequent_aligned_lv = modal,
         fraction_folds_most_frequent_lv = as.numeric(max(counts)) / .N,
         n_distinct_aligned_lvs_across_folds = uniqueN(valid_reference),
         median_alignment_abs_cor = median(valid_similarity),
         n_folds_with_aligned_lv = length(valid_reference), n_folds = .N)
  }
}, by = .(dataset, method, cell_type)]

hierarchical_bootstrap <- function(data, value_col, statistic, replicates, seed) {
  usable <- data[is.finite(get(value_col))]
  dataset_ids <- unique(usable$dataset)
  if (!nrow(usable) || !length(dataset_ids)) return(c(NA_real_, NA_real_))
  set.seed(seed)
  estimates <- replicate(replicates, {
    sampled_datasets <- sample(dataset_ids, length(dataset_ids), replace = TRUE)
    values <- unlist(lapply(sampled_datasets, function(dataset_id) {
      dataset_values <- usable[dataset == dataset_id, get(value_col)]
      sample(dataset_values, length(dataset_values), replace = TRUE)
    }), use.names = FALSE)
    statistic(values, na.rm = TRUE)
  })
  quantile(estimates, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

bootstrap_mean_interval <- function(values, replicates, seed) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(NA_real_, NA_real_))
  set.seed(seed)
  estimates <- replicate(replicates, mean(sample(values, length(values), replace = TRUE)))
  quantile(estimates, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

# Dataset-level summaries for compact paper figures. Intervals describe
# variation across the valid cell types within each dataset.
thresholded_dataset_summary <- thresholded_metrics[, {
  valid_rows <- .SD[valid_all_folds == TRUE]
  seed_offset <- match(.BY$dataset, DATASETS) * 10L
  pearson_ci <- bootstrap_mean_interval(
    valid_rows$pearson_r, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + seed_offset
  )
  r2_ci <- bootstrap_mean_interval(
    valid_rows$predictive_r2, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + seed_offset + 1L
  )
  mae_ci <- bootstrap_mean_interval(
    valid_rows$mae, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + seed_offset + 2L
  )
  list(
    n_valid = nrow(valid_rows), n_total = .N,
    mean_pearson_r = mean(valid_rows$pearson_r, na.rm = TRUE),
    pearson_ci_low = pearson_ci[[1L]], pearson_ci_high = pearson_ci[[2L]],
    mean_predictive_r2 = mean(valid_rows$predictive_r2, na.rm = TRUE),
    predictive_r2_ci_low = r2_ci[[1L]], predictive_r2_ci_high = r2_ci[[2L]],
    mean_mae = mean(valid_rows$mae, na.rm = TRUE),
    mae_ci_low = mae_ci[[1L]], mae_ci_high = mae_ci[[2L]],
    bootstrap_replicates = BOOTSTRAP_REPLICATES,
    bootstrap_seed = BOOTSTRAP_SEED + seed_offset
  )
}, by = .(dataset, method)]

# Overall intervals preserve the dataset hierarchy before resampling cell types.
thresholded_pearson_ci <- hierarchical_bootstrap(
  valid_thresholded, 'pearson_r', mean, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 100L
)
thresholded_r2_ci <- hierarchical_bootstrap(
  valid_thresholded, 'predictive_r2', mean, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 101L
)
thresholded_mae_ci <- hierarchical_bootstrap(
  valid_thresholded, 'mae', mean, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 102L
)
thresholded_summary[, `:=`(
  mean_pearson_ci_low = thresholded_pearson_ci[[1L]],
  mean_pearson_ci_high = thresholded_pearson_ci[[2L]],
  mean_predictive_r2_ci_low = thresholded_r2_ci[[1L]],
  mean_predictive_r2_ci_high = thresholded_r2_ci[[2L]],
  mean_mae_ci_low = thresholded_mae_ci[[1L]],
  mean_mae_ci_high = thresholded_mae_ci[[2L]],
  bootstrap_replicates = BOOTSTRAP_REPLICATES,
  bootstrap_seed = BOOTSTRAP_SEED + 100L
)]

valid_metrics <- oof_metrics[is.finite(pearson_r)]
median_ci <- hierarchical_bootstrap(valid_metrics, 'pearson_r', median, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED)
mean_ci <- hierarchical_bootstrap(valid_metrics, 'pearson_r', mean, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 1L)
valid_test_lv <- oof_metrics[is.finite(median_test_lv_cor)]
test_lv_median_ci <- hierarchical_bootstrap(valid_test_lv, 'median_test_lv_cor', median, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 2L)
test_lv_mean_ci <- hierarchical_bootstrap(valid_test_lv, 'median_test_lv_cor', mean, BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED + 3L)
method_summary <- valid_metrics[, .(
  n_dataset_cell_types = .N, n_datasets = uniqueN(dataset),
  mean_pearson_r = mean(pearson_r), median_pearson_r = median(pearson_r),
  n_missing_pearson_r = nrow(oof_metrics) - .N
), by = method]
method_summary[, `:=`(
  mean_ci_low = mean_ci[[1L]], mean_ci_high = mean_ci[[2L]],
  median_ci_low = median_ci[[1L]], median_ci_high = median_ci[[2L]],
  bootstrap_replicates = BOOTSTRAP_REPLICATES, bootstrap_seed = BOOTSTRAP_SEED
)]
method_summary[, `:=`(
  n_test_lv_results = nrow(valid_test_lv),
  mean_test_lv_cor = mean(valid_test_lv$median_test_lv_cor),
  median_test_lv_cor = median(valid_test_lv$median_test_lv_cor),
  mean_test_lv_ci_low = test_lv_mean_ci[[1L]],
  mean_test_lv_ci_high = test_lv_mean_ci[[2L]],
  median_test_lv_ci_low = test_lv_median_ci[[1L]],
  median_test_lv_ci_high = test_lv_median_ci[[2L]]
)]

fwrite(fold_membership, file.path(OUT, 'fold_membership.csv'))
fwrite(fold_calibrations, file.path(OUT, 'fold_calibrations.csv'))
fwrite(oof_predictions, file.path(OUT, 'oof_predictions.csv'))
fwrite(oof_metrics, file.path(OUT, 'oof_metrics.csv'))
fwrite(thresholded_metrics, file.path(OUT, 'thresholded_metrics.csv'))
fwrite(thresholded_summary, file.path(OUT, 'thresholded_summary.csv'))
fwrite(thresholded_dataset_summary, file.path(OUT, 'thresholded_dataset_summary.csv'))
fwrite(lv_selection_stability, file.path(OUT, 'lv_selection_stability.csv'))
fwrite(method_summary, file.path(OUT, 'method_summary.csv'))

method_summary
thresholded_summary
thresholded_dataset_summary[order(dataset)]
thresholded_metrics[order(dataset, cell_type)]
oof_metrics[order(dataset, cell_type)]
lv_selection_stability[order(dataset, cell_type)]
