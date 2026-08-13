#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CoGAPS)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
gtex_data <- readRDS(required_arg(args, "df_gtex_fbm_filt"))
K <- as.integer(readRDS(required_arg(args, "k")))
out_dir <- required_arg(args, "out_dir")
n_iterations <- as.integer(args$n_iterations %||% 5000L)
n_threads <- as.integer(args$n_threads %||% 4L)
seed <- as.integer(args$seed %||% 123L)

cat("Dimensions:", dim(gtex_data), "\n")
cat("nPatterns (K):", K, "\n")

# CoGAPS needs a non-negative matrix.  The pseudobulk pipeline feeds it
# log1p(CPM) (scripts/pseudobulk/cogaps.R); GTEx cannot match that, and this
# divergence is deliberate rather than an oversight:
#   * GTEx source data is TPM, not CPM, and CLAMP::cleanFBM already applies
#     log2(x + 1) upstream, so the pseudobulk expression is not reproducible here;
#   * clamp_gtex z-scores the matrix in place (scripts/gtex/clamp.R) and never
#     writes a pre-z-score filtered copy, so no log-scale matrix exists on disk to
#     hand to this rule.  Producing one means re-running clamp_gtex and
#     regenerating CLAMPbase/CLAMPfull, cascading into every GTEx analysis.
# Since GTEx CoGAPS does not converge inside its 7-day budget and is excluded
# from full_models_gtex, harmonizing the input would cost a full GTEx refit to
# change a model nobody consumes.  A global min-shift is used instead.
gtex_data_shifted <- gtex_data - min(gtex_data)

# Calculate nSets: each subset should have 1000-5000 genes
n_genes <- nrow(gtex_data_shifted)
nSets <- ceiling(n_genes / 2500)
cat("Number of genes:", n_genes, "\n")
cat("nSets (subsets):", nSets, "\n")

params <- CogapsParams(
  nPatterns = K,
  nIterations = n_iterations,
  seed = seed,
  distributed = "genome-wide"
)
params <- setDistributedParams(params, nSets = nSets)

cogapsresult <- CoGAPS(gtex_data_shifted, params, nThreads = n_threads, outputFrequency = 10000)

B <- t(cogapsresult@sampleFactors) # K x samples
sample_names <- colnames(gtex_data)
rownames(B) <- paste0("LV", seq_len(nrow(B)))
colnames(B) <- sample_names

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(B, file = file.path(out_dir, "gtex_B.csv"), quote = FALSE)
saveRDS(cogapsresult, file.path(out_dir, "cogaps_model.rds"))

message("CoGAPS saved -> ", out_dir)
