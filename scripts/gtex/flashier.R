#!/usr/bin/env Rscript
# flashier GTEx models
suppressPackageStartupMessages({
  library(flashier)
  library(ebnm)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
gtex_data <- readRDS(required_arg(args, "df_gtex_fbm_filt"))
K <- as.integer(readRDS(required_arg(args, "k")))
out_dir <- required_arg(args, "out_dir")
backfit_maxiter <- as.integer(args$backfit_maxiter %||% 20L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

X <- t(as.matrix(gtex_data))

stopifnot(is.numeric(K), length(K) == 1)
stopifnot(K <= min(nrow(X), ncol(X)))

fl <- flash_init(X, var_type = 0L)
fl <- flash_greedy(
  fl,
  Kmax = as.integer(K),
  ebnm_fn = c(ebnm_point_exponential, ebnm_point_normal),
  verbose = 1L
)
fl <- flash_backfit(fl, verbose = 1L, maxiter = backfit_maxiter)

fit <- flash_fit(fl)
pm1 <- flash_fit_get_pm(fit, n = 1)
F_scores <- pm1

B <- t(F_scores)

sample_names <- colnames(gtex_data)
stopifnot(nrow(F_scores) == length(sample_names))

colnames(F_scores) <- paste0("LV", seq_len(ncol(F_scores)))
rownames(F_scores) <- sample_names

rownames(B) <- colnames(F_scores)
colnames(B) <- sample_names

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(B, file = file.path(out_dir, "gtex_B.csv"), quote = FALSE)
saveRDS(fl, file.path(out_dir, "flashier_model.rds"))

message("flashier saved -> ", out_dir)
