#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CLAMP)
  library(data.table)
  library(dplyr)
  library(rsvd)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
norm <- read_csv_matrix(required_arg(args, "norm"))
k <- read_k(required_arg(args, "k"))
gmt <- read_gmt_file(required_arg(args, "gmt"))
pathways <- CLAMP::gmtListToSparseMat(list(BP = gmt))
base_dir <- required_arg(args, "base_dir")
full_dir <- required_arg(args, "full_dir")
max_iter <- as.integer(args$max_iter %||% 500L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

genes <- rownames(norm)
samples <- colnames(norm)
message("norm: ", nrow(norm), " genes x ", ncol(norm), " samples")
message("k = ", k)

# SVD (needed by CLAMPbase and CLAMPfull)
n_genes <- nrow(norm)
n_samples <- ncol(norm)
svd_k <- floor((min(n_genes, n_samples) - 1) / 4)
svdres <- rsvd::rsvd(norm, k = svd_k)

# CLAMPbase
message("Running CLAMPbase ...")
base <- CLAMP::CLAMPbase(Y = norm, svdres = svdres, clamp_k = k, trace = FALSE)
base$Z <- data.frame(base$Z); rownames(base$Z) <- genes
base$B <- data.frame(base$B); colnames(base$B) <- samples

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(base$B, file.path(base_dir, "B.csv"))
write.csv(base$Z, file.path(base_dir, "Z.csv"))
saveRDS(base, file.path(base_dir, "CLAMPbase.rds"))
message("CLAMPbase saved -> ", base_dir)

# Match BP prior to dataset genes
matched <- CLAMP::getMatchedPathwayMat(pathways, genes)
message("Prior matched: ", nrow(matched), " genes x ", ncol(matched), " pathways")

# CLAMPfull
message("Running CLAMPfull ...")
full <- CLAMP::CLAMPfull(
  Y                 = norm,
  svdres            = svdres,
  priorMat          = matched,
  clamp.base.result = base,
  use_cpp           = TRUE,
  trace             = FALSE,
  max.iter          = max_iter,
  clamp_k           = k
)
full$Z <- data.frame(full$Z); rownames(full$Z) <- genes
full$B <- data.frame(full$B); colnames(full$B) <- samples

dir.create(full_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(full$B, file.path(full_dir, "B.csv"))
write.csv(full$Z, file.path(full_dir, "Z.csv"))
fwrite(data.frame(L2 = as.numeric(full$L2)), file.path(full_dir, "L2.csv"))

if (!is.null(full$summary)) {
  sumdf <- as.data.frame(full$summary) %>%
    dplyr::rename(LV = LV_index) %>%
    dplyr::mutate(LV = paste0("LV", LV))
} else {
  sumdf <- data.frame()
}
write.csv(sumdf, file.path(full_dir, "summary.csv"), row.names = FALSE)
saveRDS(full, file.path(full_dir, "CLAMPfull.rds"))
message("CLAMPfull saved -> ", full_dir)
