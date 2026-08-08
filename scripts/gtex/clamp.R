#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(bigstatsr)
  library(data.table)
  library(dplyr)
  library(rsvd)
  library(CLAMP)
  library(PCAtools)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
raw_gct <- required_arg(args, "raw_gct")
gmt_path <- required_arg(args, "gmt")
out_dir <- required_arg(args, "out_dir")
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
block_size <- as.integer(args$chunk_size %||% 100L)
n_cores <- as.integer(args$n_cores %||% 1L)
max_iter <- as.integer(args$max_iter %||% 500L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Aggregate raw GTEx TPM matrix to gene level
exprs_data <- read.table(raw_gct, header = TRUE, sep = "\t", skip = 2, check.names = FALSE)
gtex <- as.data.table(exprs_data)
aggregated_gtex <- gtex[, lapply(.SD, sum), by = Description, .SDcols = is.numeric]

genes <- aggregated_gtex$Description
samples <- colnames(aggregated_gtex[, -1])
data_mat <- as.matrix(aggregated_gtex[, -1])

# Create the FBM
fbm_files <- c(
  file.path(out_dir, "FBMgtex.bk"),
  file.path(out_dir, "FBMgtex_preproc.bk"),
  file.path(out_dir, "FBMgtex_preproc_filtered.bk"),
  file.path(out_dir, "FBMgtex_preproc_filtered.rds")
)
for (f in fbm_files) {
  if (file.exists(f)) {
    res <- unlink(f)
    if (identical(res, 0L)) {
      message("Removed existing FBM file: ", f)
    } else {
      warning("Failed to remove: ", f)
    }
  }
}

fbm_file <- file.path(out_dir, "FBMgtex")
gtexFBM <- FBM(nrow = nrow(data_mat), ncol = ncol(data_mat), backingfile = fbm_file, create_bk = TRUE)

n_blocks <- ceiling(nrow(aggregated_gtex) / block_size)
for (i in 1:n_blocks) {
  start_row <- (i - 1) * block_size + 1
  end_row <- min(i * block_size, nrow(data_mat))
  gtexFBM[start_row:end_row, ] <- as.matrix(data_mat[start_row:end_row, ])
}

# Preprocess and z-score FBM
prep_gtex <- preprocessCLAMPFBM(
  fbm = gtexFBM,
  mean_cutoff = mean_cutoff,
  var_cutoff = var_cutoff,
  ncores = n_cores
)

gtex_fbm_filt <- prep_gtex$fbm_filtered
gtex_rowStats <- prep_gtex$rowStats

zscoreCLAMPFBM(gtex_fbm_filt, gtex_rowStats, ncores = n_cores)
gtex_genes <- genes[prep_gtex$kept_rows]

saveRDS(samples, file = file.path(out_dir, "gtex_samples.rds"))
saveRDS(gtex_genes, file = file.path(out_dir, "gtex_genes.rds"))
saveRDS(gtex_fbm_filt, file = file.path(out_dir, "gtex_fbm_filt.rds"))

df_gtex_fbm_filt <- as.data.frame(as.matrix(gtex_fbm_filt[]))
colnames(df_gtex_fbm_filt) <- samples
rownames(df_gtex_fbm_filt) <- gtex_genes

saveRDS(df_gtex_fbm_filt, file = file.path(out_dir, "df_gtex_fbm_filt.rds"))
write.csv(df_gtex_fbm_filt, file = file.path(out_dir, "df_gtex_fbm_filt.csv"))

# SVD computation
g_fb <- nrow(gtex_fbm_filt)
samples_fb <- ncol(gtex_fbm_filt)
SVD_K_gtex <- min(g_fb, samples_fb) - 1
SVD_K_gtex <- floor(SVD_K_gtex / 4)
message("Using SVD K = ", SVD_K_gtex)

gtex_svdRes <- rsvd(gtex_fbm_filt[], k = SVD_K_gtex)
saveRDS(gtex_svdRes, file = file.path(out_dir, "gtex_svdRes.rds"))

# Estimate K for CLAMP
n_genes_gtex <- nrow(gtex_fbm_filt)
n_samples_gtex <- ncol(gtex_fbm_filt)
eigenvalues <- sort(gtex_svdRes$d^2 / (n_samples_gtex - 1), decreasing = TRUE)
noise_gd <- median(eigenvalues)
CLAMP_K_gtex <- PCAtools::chooseGavishDonoho(
  .dim = c(n_genes_gtex, n_samples_gtex),
  var.explained = eigenvalues,
  noise = noise_gd
) * 2
message("Inferred CLAMP K = ", CLAMP_K_gtex)

saveRDS(CLAMP_K_gtex, file = file.path(out_dir, "CLAMP_K_gtex.rds"))
write.csv(
  as.data.frame(CLAMP_K_gtex),
  file = file.path(out_dir, "CLAMP_K_gtex.csv"),
  row.names = TRUE
)

# CLAMPbase
message("Running CLAMPbase ...")
gtex_baseRes <- CLAMPbase(
  Y = gtex_fbm_filt,
  svdres = gtex_svdRes,
  trace = TRUE,
  clamp_k = CLAMP_K_gtex
)
gtex_baseRes$Z <- data.frame(gtex_baseRes$Z)
rownames(gtex_baseRes$Z) <- gtex_genes
gtex_baseRes$B <- data.frame(gtex_baseRes$B)
colnames(gtex_baseRes$B) <- samples

saveRDS(gtex_baseRes, file = file.path(out_dir, "CLAMPbase.rds"))
base_dir <- file.path(out_dir, "CLAMPbase")
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(gtex_baseRes$B, file.path(base_dir, "B.csv"))
write.csv(gtex_baseRes$Z, file.path(base_dir, "Z.csv"))

# Prepare pathway priors (pinned GO:BP file, see common.R::build_go_bp_prior)
prior <- build_go_bp_prior(gtex_genes, gmt_path)
gtex_matched <- prior$matched

# CLAMPfull
message("Running CLAMPfull ...")
gtex_fullRes <- CLAMPfull(
  Y = gtex_fbm_filt,
  svdres = gtex_svdRes,
  priorMat = gtex_matched,
  clamp.base.result = gtex_baseRes,
  use_cpp = TRUE,
  trace = TRUE,
  max.iter = max_iter,
  clamp_k = CLAMP_K_gtex
)
gtex_fullRes$Z <- data.frame(gtex_fullRes$Z)
rownames(gtex_fullRes$Z) <- gtex_genes
gtex_fullRes$B <- data.frame(gtex_fullRes$B)
colnames(gtex_fullRes$B) <- samples
gtex_fullRes$summary <- gtex_fullRes$summary %>%
  dplyr::rename(LV = LV_index) %>%
  dplyr::mutate(LV = paste0("LV", LV))
colnames(gtex_fullRes$B) <- samples

saveRDS(gtex_fullRes, file = file.path(out_dir, "CLAMPfull.rds"))
full_dir <- file.path(out_dir, "CLAMPfull")
dir.create(full_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(gtex_fullRes$B, file.path(full_dir, "B.csv"))
write.csv(gtex_fullRes$Z, file.path(full_dir, "Z.csv"))
write.csv(gtex_fullRes$summary, file.path(full_dir, "summary.csv"))

message("GTEx CLAMP done -> ", out_dir)
