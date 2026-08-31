#!/usr/bin/env Rscript
# Builds the filtered, z-scored recount2 matrix used by coverage fits.

suppressPackageStartupMessages({
  library(bigstatsr)
  library(dplyr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
plier_rds <- required_arg(args, "plier_rds")
rpkm_rds <- required_arg(args, "rpkm_rds")
out_dir <- required_arg(args, "out_dir")
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)
block_size <- as.integer(args$block_size %||% 100L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

ensure_dir(out_dir)

prep <- readRDS(plier_rds)
gene_symbols <- rownames(prep$rpkm.cm)
samples <- colnames(prep$rpkm.cm)
rm(prep)
message("recount2 target axes: ", length(gene_symbols), " genes x ", length(samples), " samples")

rpkm <- readRDS(rpkm_rds)
rpkm$ensembl_gene_id <- vapply(strsplit(rpkm$ENSG, ".", fixed = TRUE), `[[`, character(1), 1L)

suppressPackageStartupMessages(library(org.Hs.eg.db))
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(rpkm$ensembl_gene_id),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)
gene_map <- gene_map %>%
  dplyr::rename(ensembl_gene_id = ENSEMBL, hgnc_symbol = SYMBOL) %>%
  dplyr::filter(stats::complete.cases(.))

rpkm <- dplyr::inner_join(gene_map, rpkm, by = "ensembl_gene_id",
                          relationship = "many-to-many")
rownames(rpkm) <- make.names(rpkm$hgnc_symbol, unique = TRUE)
rpkm <- rpkm %>% dplyr::select(-c(ensembl_gene_id:ENSG))

data_mat <- rpkm[gene_symbols, samples]
rm(rpkm)
if (!identical(dim(data_mat), c(length(gene_symbols), length(samples)))) {
  stop("Expression matrix does not match the requested gene/sample axes")
}

fbm_prefix <- file.path(out_dir, "FBMrecount2")
for (suffix in c(".bk", "_preproc.bk", "_preproc_filtered.bk", "_preproc_filtered.rds")) {
  path <- paste0(fbm_prefix, suffix)
  if (file.exists(path)) unlink(path)
}

n_genes <- length(gene_symbols)
n_samples <- length(samples)
fbm <- FBM(nrow = n_genes, ncol = n_samples, backingfile = fbm_prefix, create_bk = TRUE)

n_blocks <- ceiling(n_genes / block_size)
for (i in seq_len(n_blocks)) {
  lo <- (i - 1L) * block_size + 1L
  hi <- min(i * block_size, n_genes)
  fbm[lo:hi, ] <- as.matrix(data_mat[lo:hi, ])
}
rm(data_mat)

prep_recount2 <- preprocessCLAMPFBM(
  fbm = fbm,
  mean_cutoff = mean_cutoff,
  var_cutoff = var_cutoff
)
fbm_filt <- prep_recount2$fbm_filtered
zscoreCLAMPFBM(fbm_filt, prep_recount2$rowStats)

recount2_genes <- gene_symbols[prep_recount2$kept_rows]
message("Filtered recount2: ", nrow(fbm_filt), " genes x ", ncol(fbm_filt), " samples")

saveRDS(samples, file.path(out_dir, "recount2_samples.rds"))
saveRDS(recount2_genes, file.path(out_dir, "recount2_genes.rds"))
saveRDS(fbm_filt, file.path(out_dir, "recount2_fbm_filt.rds"))

message("Wrote recount2 inputs to ", out_dir)
