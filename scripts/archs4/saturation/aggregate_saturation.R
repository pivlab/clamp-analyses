#!/usr/bin/env Rscript
# Aggregate per-cell saturation ORA summaries into long/panel-ready tables.

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

summary_files <- list.files(ora_root, pattern = "^summary\\.csv$", recursive = TRUE, full.names = TRUE)
canonical <- "/rs[0-9]+/k[0-9]+/seed[0-9]+/[^./]+/[^./]+/summary\\.csv$"
summary_files <- grep(canonical, summary_files, value = TRUE)
# Saturation scores against canonical + CellMarker only; drop any Reactome
# directories left on disk from earlier runs.
summary_files <- summary_files[!grepl("/reactome/", summary_files, fixed = TRUE)]
if (!length(summary_files)) stop("No published ORA summaries found below ", ora_root)
complete_files <- file.path(dirname(summary_files), "complete")
if (any(!file.exists(complete_files))) {
  stop("Incomplete ORA directories: ", paste(dirname(summary_files)[!file.exists(complete_files)], collapse = ", "))
}

saturation <- rbindlist(lapply(summary_files, function(path) {
  row <- fread(path)
  # K is recovered from the path: .../rs{f}/k{K}/seed{s}/{model}/{database}/summary.csv
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
