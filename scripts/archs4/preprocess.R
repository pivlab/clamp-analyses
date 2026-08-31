#!/usr/bin/env Rscript

# Builds the filtered, z-scored ARCHS4 expression matrix.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(hdf5r)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
h5_path <- required_arg(args, "h5")
gene_lengths_path <- required_arg(args, "gene_lengths")
out_dir <- required_arg(args, "out_dir")
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
sc_max <- as.numeric(args$single_cell_probability_max %||% 0.5)
block_size <- as.integer(args$block_size %||% 100L)
n_cores <- as.integer(args$n_cores %||% 1L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

ensure_dir(out_dir)

h5 <- H5File$new(h5_path, mode = "r")
on.exit(h5$close_all(), add = TRUE)

dset <- h5[["/data/expression"]]
gene_symbols <- h5[["/meta/genes/symbol"]]$read()
sample_names <- h5[["/meta/samples/geo_accession"]]$read()
sc_probability <- h5[["/meta/samples/singlecellprobability"]]$read()

n_genes_raw <- length(gene_symbols)
n_samples_raw <- length(sample_names)
message("ARCHS4 raw: ", n_genes_raw, " genes x ", n_samples_raw, " samples")

saveRDS(
  list(gene_symbols = gene_symbols, sample_names = sample_names,
       single_cell_probability = sc_probability),
  file.path(out_dir, "metadata.rds")
)
saveRDS(sample_names, file.path(out_dir, "all_samples.rds"))

gene_lengths <- readRDS(gene_lengths_path)

gene_symbols_idx <- which(gene_symbols %in% names(gene_lengths))
gene_symbols_thin <- gene_symbols[gene_symbols_idx]
gene_symbols_unique <- sort(unique(gene_symbols_thin))
gene_lengths <- gene_lengths[gene_symbols_unique]
n_genes_thin <- length(gene_symbols_unique)
stopifnot(!anyNA(gene_lengths), n_genes_thin > 0L)
message("Genes with a known length, after collapsing duplicates: ", n_genes_thin)

fbm_file <- file.path(out_dir, "fbm")
if (file.exists(paste0(fbm_file, ".bk"))) file.remove(paste0(fbm_file, ".bk"))

fbm_obj <- FBM(nrow = n_genes_thin, ncol = n_samples_raw,
               backingfile = fbm_file, create_bk = TRUE)

n_blocks <- ceiling(n_samples_raw / block_size)
pb <- txtProgressBar(min = 0, max = n_blocks, style = 3)
n_bad_blocks <- 0L

for (i in seq_len(n_blocks)) {
  setTxtProgressBar(pb, i)
  start_col <- (i - 1L) * block_size + 1L
  end_col <- min(i * block_size, n_samples_raw)

  raw_block <- tryCatch(dset[start_col:end_col, ], error = function(e) NULL)

  if (is.null(raw_block)) {
    n_bad_blocks <- n_bad_blocks + 1L
    message("\nUnreadable block ", i, " (samples ", start_col, "-", end_col, "); substituting zeros")
    raw_block <- matrix(1e-10, nrow = end_col - start_col + 1L, ncol = n_genes_raw)
  }

  raw_block <- raw_block[, gene_symbols_idx, drop = FALSE]
  raw_block_summed <- rowsum(t(raw_block), group = gene_symbols_thin)
  if (i == 1L) stopifnot(identical(rownames(raw_block_summed), gene_symbols_unique))

  tpm <- tpm_norm(t(raw_block_summed), gene_lengths)
  fbm_obj[, start_col:end_col] <- as.matrix(t(tpm))
}
close(pb)
if (n_bad_blocks > 0L) message("Unreadable blocks substituted: ", n_bad_blocks)

with_single_threaded_blas(n_cores, cleanFBM(fbm_obj, ncores = n_cores))
row_stats <- with_single_threaded_blas(n_cores, computeRowStatsFBM(fbm_obj, ncores = n_cores))

samples_idx <- which(sc_probability < sc_max)
n_samples <- length(samples_idx)
message("Samples below single-cell probability ", sc_max, ": ", n_samples)

filter_result <- filterFBM(
  fbm_obj,
  row_stats,
  keep_samples_idx = samples_idx,
  mean_cutoff = mean_cutoff,
  var_cutoff = var_cutoff,
  backingfile = paste0(fbm_file, "_filtered")
)
fbm_filtered <- filter_result$fbm_filtered

gene_symbols_kept <- gene_symbols_unique[filter_result$kept_rows]
message("Filtered dataset: ", nrow(fbm_filtered), " genes x ", ncol(fbm_filtered), " samples")

saveRDS(
  list(
    gene_symbols_thin = gene_symbols_kept,
    gene_lengths = gene_lengths,
    n_genes_thin = length(gene_symbols_kept),
    n_samples = n_samples
  ),
  file.path(out_dir, "metadata_filtered.rds")
)

row_stats$row_means <- row_stats$row_means[filter_result$kept_rows]
row_stats$row_variances <- row_stats$row_variances[filter_result$kept_rows]
zscoreFBM(fbm_filtered, rowStats = row_stats, chunk_size = block_size)

message("Wrote ", file.path(out_dir, "fbm_filtered.bk"))
