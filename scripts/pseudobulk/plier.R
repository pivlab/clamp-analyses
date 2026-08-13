#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CLAMP)
  library(PLIER)
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
out <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

samples <- colnames(norm)
genes <- rownames(norm)
message("norm: ", nrow(norm), " genes x ", ncol(norm), " samples")
message("k = ", k)

# SVD (for PLIER warm start)
n_genes <- nrow(norm)
n_samples <- ncol(norm)
svd_k <- max(floor((min(n_genes, n_samples) - 1) / 4), k, 2L)
svdres <- rsvd::rsvd(norm, k = svd_k)

# Match BP prior to dataset genes
matched <- CLAMP::getMatchedPathwayMat(pathways, genes)
chatObj <- CLAMP::getChat(matched)
message("Prior matched: ", nrow(matched), " genes x ", ncol(matched), " pathways")

# PLIER
message("Running PLIER ...")
plier_res <- PLIER::PLIER(
  norm,
  as.matrix(matched),
  svdres     = svdres,
  Chat       = as.matrix(chatObj),
  doCrossval = TRUE,
  k          = k
)

colnames(plier_res$Z) <- paste0("LV", seq_len(ncol(plier_res$Z)))
colnames(plier_res$B) <- samples

if (!is.null(plier_res$summary)) {
  plier_res$summary <- plier_res$summary %>%
    dplyr::rename(LV = `LV index`) %>%
    dplyr::mutate(LV = paste0("LV", LV))
} else {
  plier_res$summary <- data.frame()
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(plier_res$B, file.path(out, "B.csv"))
write.csv(plier_res$Z, file.path(out, "Z.csv"))
write.csv(plier_res$summary, file.path(out, "summary.csv"), row.names = FALSE)
saveRDS(plier_res, file.path(out, "PLIER.rds"))
message("PLIER saved -> ", out)
