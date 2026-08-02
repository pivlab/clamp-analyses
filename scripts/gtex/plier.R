#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(PLIER)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
fbm_filt <- readRDS(required_arg(args, "fbm_filt"))
svdRes <- readRDS(required_arg(args, "svd_res"))
CLAMP_K_gtex <- readRDS(required_arg(args, "k"))
gtex_genes <- readRDS(required_arg(args, "genes"))
samples <- readRDS(required_arg(args, "samples"))
gmt_path <- required_arg(args, "gmt")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Prepare pathway priors (pinned GO:BP file, see common.R::build_go_bp_prior)
prior <- build_go_bp_prior(gtex_genes, gmt_path)
gtex_matched <- prior$matched
gtex_chatObj <- prior$chat

# PLIER
message("Running PLIER ...")
gtex_plier <- PLIER::PLIER(
  fbm_filt[],
  as.matrix(gtex_matched),
  svdres = svdRes,
  Chat = as.matrix(gtex_chatObj),
  doCrossval = TRUE,
  k = CLAMP_K_gtex
)

colnames(gtex_plier$Z) <- paste0("LV", seq_len(ncol(gtex_plier$Z)))
gtex_plier$summary <- gtex_plier$summary %>%
  dplyr::rename(LV = `LV index`) %>%
  dplyr::mutate(LV = paste0("LV", LV))
colnames(gtex_plier$B) <- samples

saveRDS(gtex_plier, file = file.path(out_dir, "PLIER.rds"))
model_dir <- file.path(out_dir, "PLIER")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(gtex_plier$B, file.path(model_dir, "B.csv"))
write.csv(gtex_plier$Z, file.path(model_dir, "Z.csv"))
write.csv(gtex_plier$summary, file.path(model_dir, "summary.csv"))

message("PLIER saved -> ", out_dir)
