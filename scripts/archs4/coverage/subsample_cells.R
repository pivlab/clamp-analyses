#!/usr/bin/env Rscript

# Builds one study-level sample of the filtered ARCHS4 compendium

suppressPackageStartupMessages({
  library(bigstatsr)
  library(hdf5r)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "..", "common.R"))

args <- parse_cli()
h5_path <- required_arg(args, "h5")
metadata_path <- required_arg(args, "metadata_filtered")
fbm_filtered_bk <- required_arg(args, "fbm_filtered")
out_dir <- required_arg(args, "out_dir")
fraction <- as.numeric(required_arg(args, "fraction"))
seed <- as.integer(required_arg(args, "seed"))
run_index <- as.integer(args$run_index %||% 1L)
sc_max <- as.numeric(args$single_cell_probability_max %||% 0.5)
block_size <- as.integer(args$block_size %||% 1000L)

if (!is.finite(fraction) || fraction <= 0 || fraction > 100) {
  stop("--fraction must be in (0, 100]; got ", fraction)
}
ensure_dir(out_dir)

metadata <- readRDS(metadata_path)
n_genes <- length(metadata$gene_symbols_thin)

h5 <- H5File$new(h5_path, mode = "r")
on.exit(h5$close_all(), add = TRUE)
all_samples <- h5[["/meta/samples/geo_accession"]]$read()
all_series <- h5[["/meta/samples/series_id"]]$read()
sc_probability <- h5[["/meta/samples/singlecellprobability"]]$read()
h5$close_all()

keep <- which(sc_probability < sc_max)
sample_names_all <- all_samples[keep]
series_all <- all_series[keep]
n_total <- length(keep)

fbm_filtered <- attach_fbm(sub("\\.bk$", "", fbm_filtered_bk), n_genes, n_total)
fbm_filtered$is_read_only <- TRUE

studies <- unique(series_all)
n_studies_total <- length(studies)
n_samples_target <- floor(fraction / 100 * n_total)

message("Filtered universe: ", n_total, " samples across ", n_studies_total, " studies")
message("Target for ", fraction, "%: ", n_samples_target, " samples")

if (fraction == 100) {
  selected_studies <- studies
  sample_idx <- seq_len(n_total)
} else {
  set.seed(seed)
  shuffled <- sample(studies)
  study_sizes <- table(series_all)
  sizes <- as.integer(study_sizes[shuffled])
  take <- which(cumsum(sizes) >= n_samples_target)[1L]
  if (is.na(take)) take <- length(shuffled)
  selected_studies <- shuffled[seq_len(take)]
  sample_idx <- sort(which(series_all %in% selected_studies))
}

sample_names <- sample_names_all[sample_idx]
n_samples <- length(sample_idx)
message("Drew ", n_samples, " samples from ", length(selected_studies), " studies")

if (n_samples < 1L) stop("Subsample is empty at fraction ", fraction)

subsample_info <- list(
  run = run_index,
  seed = seed,
  coverage = fraction / 100,
  n_studies_total = n_studies_total,
  n_studies = length(selected_studies),
  study_ids = selected_studies,
  n_samples = n_samples,
  n_samples_target = n_samples_target,
  sample_idx = sample_idx,
  sample_names = sample_names,
  sampling_method = "study_sampling",
  sampling_universe = "filtered_singlecellprobability",
  single_cell_probability_max = sc_max,
  filtered_universe_samples = n_total
)
saveRDS(subsample_info, file.path(out_dir, "subsample_info.rds"))

saveRDS(
  list(
    gene_symbols_thin = metadata$gene_symbols_thin,
    gene_lengths = metadata$gene_lengths,
    n_genes_thin = n_genes,
    n_samples = n_samples
  ),
  file.path(out_dir, "metadata_cell.rds")
)
saveRDS(sample_names, file.path(out_dir, "samples.rds"))

if (fraction == 100) {
  message("Fraction 100: reusing ", fbm_filtered_bk, "; no subsampled copy written")
} else {
  backing <- file.path(out_dir, "fbm_subsampled")
  if (file.exists(paste0(backing, ".bk"))) file.remove(paste0(backing, ".bk"))
  subsampled <- FBM(nrow = n_genes, ncol = n_samples,
                    backingfile = backing, create_bk = TRUE)
  n_blocks <- ceiling(n_samples / block_size)
  for (i in seq_len(n_blocks)) {
    lo <- (i - 1L) * block_size + 1L
    hi <- min(i * block_size, n_samples)
    subsampled[, lo:hi] <- fbm_filtered[, sample_idx[lo:hi], drop = FALSE]
  }
  message("Wrote ", backing, ".bk (", n_genes, " genes x ", n_samples, " samples)")
}
