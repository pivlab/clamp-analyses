#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(reticulate)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
gtex_rds <- required_arg(args, "df_gtex_fbm_filt")
k_rds <- required_arg(args, "k")
out_dir <- required_arg(args, "out_dir")
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

# Name of the env you expect users to have
env_name <- "clamp-analyses"

if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  message("Using RETICULATE_PYTHON = ", Sys.getenv("RETICULATE_PYTHON"))
} else {
  conda <- Sys.which("conda")
  if (nzchar(conda)) {
    cmd <- sprintf(
      '%s run -n %s python -c "import sys; print(sys.executable)"',
      shQuote(conda), shQuote(env_name)
    )
    py <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
    if (length(py) == 1 && nzchar(py) && file.exists(py)) {
      Sys.setenv(RETICULATE_PYTHON = py)
      message("Auto-set RETICULATE_PYTHON = ", py)
    } else {
      message("Could not resolve env python via conda. Falling back to Sys.which('python').")
      Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
    }
  } else {
    message("conda not found on PATH. Falling back to Sys.which('python').")
    Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
  }
}
py_config()

mfl <- import("mofaflex", delay_load = FALSE)
ad <- import("anndata", delay_load = FALSE)
np <- import("numpy", delay_load = FALSE)
pd <- import("pandas", delay_load = FALSE)

stopifnot(file.exists(gtex_rds), file.exists(k_rds))

gtex_data <- readRDS(gtex_rds) # genes x samples
K <- readRDS(k_rds) # scalar

stopifnot(is.numeric(K), length(K) == 1)
stopifnot(!is.null(rownames(gtex_data)), !is.null(colnames(gtex_data)))

# AnnData samples x genes
X <- t(as.matrix(gtex_data))
stopifnot(K <= min(nrow(X), ncol(X)))

# lower memory than float64
X_np <- np$array(X, dtype = "float32")
adata <- ad$AnnData(X_np)

# preserve names
adata$obs_names <- pd$Index(colnames(gtex_data)) # samples
adata$var_names <- pd$Index(rownames(gtex_data)) # genes

# MOFA-FLEX expects nested dict: group -> view -> AnnData
data_list <- dict(group_1 = dict(view_1 = adata))

model <- mfl$MOFAFLEX(
  data_list,
  mfl$DataOptions(
    scale_per_group = FALSE, # already z-scored
    plot_data_overview = FALSE,
    remove_constant_features = TRUE
  ),
  mfl$ModelOptions(
    n_factors = as.integer(K),
    likelihoods = "Normal",
    weight_prior = "Laplace"
  ),
  mfl$TrainingOptions(
    seed = as.integer(seed),
    max_epochs = 2000L,
    batch_size = 1000L
  )
)

# extract factors Z (samples x K) then convert to B (K x samples)
Z_dict <- model$get_factors(return_type = "numpy", ordered = FALSE)
Z <- py_to_r(Z_dict[["group_1"]]) # samples x K

B <- t(Z) # K x samples
rownames(B) <- paste0("LV", seq_len(nrow(B)))
colnames(B) <- colnames(gtex_data)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(B, file = file.path(out_dir, "gtex_B.csv"), quote = FALSE)
py_save_object(model, file.path(out_dir, "mofaflex_model.pkl"))

message("MOFA-FLEX base saved -> ", out_dir)
