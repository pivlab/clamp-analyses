#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) {
    if (i == length(args)) stop("Missing value for ", args[[i]])
    out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", gsub("_", "-", name))
  value
}

args <- parse_cli()
ora_root <- required_arg(args, "ora_root")
saturation_out <- required_arg(args, "saturation_out")
panel_out <- required_arg(args, "panel_out")

# Only published ORA directories count.  run_ora.R builds a complete result
# set, `complete` marker included, inside `<out>.tmp.<pid>` and publishes it
# with an atomic rename, so an interrupted run leaves behind a directory that
# looks finished.  Matching the canonical rs<f>/k<K>/seed<s>/<model>/<database>
# layout skips those instead of double-counting them.
summary_files <- list.files(ora_root, pattern = "^summary\\.csv$", recursive = TRUE, full.names = TRUE)
canonical <- "/rs[0-9]+/k[0-9]+/seed[0-9]+/[^./]+/[^./]+/summary\\.csv$"
summary_files <- grep(canonical, summary_files, value = TRUE)
if (!length(summary_files)) stop("No published ORA summaries found below ", ora_root)
complete_files <- file.path(dirname(summary_files), "complete")
if (any(!file.exists(complete_files))) {
  stop("Incomplete ORA directories: ", paste(dirname(summary_files)[!file.exists(complete_files)], collapse = ", "))
}

saturation <- rbindlist(lapply(summary_files, function(path) {
  row <- fread(path)
  # run_ora.R records dataset/fraction/seed but has no notion of K, so it is
  # recovered from the path.  `latent_variables` in the summary is the
  # independent check that the two agree.
  # .../rs{f}/k{K}/seed{s}/{model}/{database}/summary.csv
  #        -5    -4     -3      -2         -1            0
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  k_part <- parts[length(parts) - 4L]
  if (!grepl("^k[0-9]+$", k_part)) stop("Cannot read K from path: ", path)
  row[, k := as.integer(sub("^k", "", k_part))]
  row
}), fill = TRUE)

required <- c(
  "dataset", "fraction", "seed", "model", "database", "database_label",
  "n_samples", "eligible_pathways", "recovered_pathways", "recovered_percent",
  "fdr", "latent_variables", "k"
)
if (!all(required %in% names(saturation))) {
  stop("ORA summaries are missing: ", paste(setdiff(required, names(saturation)), collapse = ", "))
}
if (anyDuplicated(saturation[, .(fraction, k, seed, model, database)])) {
  stop("Duplicate fraction/K/seed/model/database ORA summaries found")
}
# The rank the model actually carries must match the rank its path claims;
# a mismatch would mean a cell was scored against the wrong model.
mismatch <- saturation[latent_variables != k]
if (nrow(mismatch)) {
  stop(
    "K from path disagrees with latent_variables in: ",
    paste(sprintf("rs%s/k%s/seed%s", mismatch$fraction, mismatch$k, mismatch$seed), collapse = ", ")
  )
}
setorder(saturation, database, fraction, k, model, seed)

panel <- saturation[, .(
  mean_recovered = mean(recovered_pathways),
  sd_recovered = if (.N > 1L) sd(recovered_pathways) else NA_real_,
  mean_percent = mean(recovered_percent),
  eligible_pathways = unique(eligible_pathways),
  n_observations = .N,
  median_samples = as.integer(median(n_samples))
), by = .(fraction, k, model, database, database_label, fdr)]
setorder(panel, database, fraction, k, model)

dir.create(dirname(saturation_out), recursive = TRUE, showWarnings = FALSE)
fwrite(saturation, saturation_out)
fwrite(panel, panel_out)
