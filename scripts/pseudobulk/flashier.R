#!/usr/bin/env Rscript
# flashier pseudobulk models

suppressPackageStartupMessages({
  library(data.table)
  library(flashier)
  library(ebnm)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()
loaded <- read_norm_and_k(required_arg(args, "norm"), required_arg(args, "k"))
backfit <- as.integer(if (is.null(args$backfit_maxiter)) 10L else args$backfit_maxiter)
set.seed(as.integer(if (is.null(args$seed)) 123L else args$seed))
fl <- flashier::flash_init(t(loaded$norm), var_type = 0L)
fl <- flashier::flash_greedy(
  fl, Kmax = loaded$k, ebnm_fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_normal),
  verbose = 0L
)
fl <- flashier::flash_backfit(fl, verbose = 0L, maxiter = backfit)

B <- t(fl$L_pm); colnames(B) <- colnames(loaded$norm); rownames(B) <- paste0("LV", seq_len(nrow(B)))
Z <- fl$F_pm; rownames(Z) <- rownames(loaded$norm); colnames(Z) <- paste0("LV", seq_len(ncol(Z)))
out <- required_arg(args, "out_dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(B, file.path(out, "B.csv")); write.csv(Z, file.path(out, "Z.csv"))
saveRDS(fl, file.path(out, "flashier_model.rds"))
message("flashier complete")
