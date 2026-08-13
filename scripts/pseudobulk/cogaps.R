#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(CoGAPS)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()
loaded <- read_norm_and_k(required_arg(args, "cpm_filt"), required_arg(args, "k"))
n_iterations <- as.integer(args$n_iterations %||% 5000L)
n_threads <- as.integer(args$n_threads %||% 4L)
parallel_workers <- if (!is.null(args$parallel_workers)) {
  as.integer(args$parallel_workers)
} else {
  NULL
}
seed <- as.integer(args$seed %||% 123L)

# cpm_filt is the filtered log2(CPM + 1) matrix from preprocess.R: non-negative and
# already on the log scale CoGAPS's Gibbs sampler expects (without a log transform it
# warns "Large values detected, is data log transformed?").  The log1p that used to be
# applied here was removed when preprocess.R started log-transforming before filtering,
# to match the GTEx pipeline -- keeping it would have transformed the data twice.
mat <- loaded$norm

# Distributed genome-wide CoGAPS: each subset should have 1000-5000 genes
nSets <- max(ceiling(nrow(mat) / 2500), 2L)
params <- CogapsParams(
  nPatterns = loaded$k,
  nIterations = n_iterations,
  seed = seed,
  distributed = "genome-wide"
)
params <- setDistributedParams(params, nSets = nSets)
message("nPatterns=", loaded$k, " nSets=", nSets, " nIterations=", n_iterations,
        " parallelWorkers=", parallel_workers %||% nSets)

bpparam <- if (is.null(parallel_workers)) {
  NULL
} else {
  BiocParallel::MulticoreParam(workers = parallel_workers)
}
cogapsresult <- CoGAPS(
  mat,
  params,
  nThreads = n_threads,
  BPPARAM = bpparam,
  outputFrequency = 10000
)

# LVs x samples / genes x LVs, matching CLAMP's B/Z convention
B <- t(cogapsresult@sampleFactors)
colnames(B) <- colnames(mat); rownames(B) <- paste0("LV", seq_len(nrow(B)))
Z <- cogapsresult@featureLoadings
rownames(Z) <- rownames(mat); colnames(Z) <- paste0("LV", seq_len(ncol(Z)))

out <- required_arg(args, "out_dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(B, file.path(out, "B.csv")); write.csv(Z, file.path(out, "Z.csv"))
saveRDS(cogapsresult, file.path(out, "cogaps_model.rds"))
message("CoGAPS complete")
