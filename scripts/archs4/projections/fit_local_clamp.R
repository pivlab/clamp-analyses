#!/usr/bin/env Rscript
# Fit a small CLAMP model on one dataset's own RNA-seq, for comparison against
# a projection of the same data into the ARCHS4 model.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
  library(rsvd)
  library(CLAMP)
  library(PCAtools)
  library(readxl)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
norm_path <- required_arg(args, "norm")
out_dir <- required_arg(args, "out_dir")
prior_names <- strsplit(required_arg(args, "prior_names"), ",")[[1]]
prior_gmts <- strsplit(required_arg(args, "prior_gmts"), ",")[[1]]
if (length(prior_names) != length(prior_gmts)) stop("--prior-names and --prior-gmts differ in length")
cellmarker_path <- required_arg(args, "cellmarker_file")
cellmarker_sheet <- args$cellmarker_sheet %||% "human"
cellmarker_term <- args$cellmarker_term_column %||% "cell_name"
cellmarker_gene <- args$cellmarker_gene_column %||% "Symbol"
k_rule <- args$k_rule %||% "n_samples_minus_1"
max_iter <- as.integer(args$max_iter %||% 200L)
seed <- as.integer(args$seed %||% 123L)

started <- Sys.time()
set.seed(seed)

norm <- read_matrix_csv(norm_path)
genes <- rownames(norm)
samples <- colnames(norm)
n_genes <- nrow(norm); n_samples <- ncol(norm)
message("norm: ", n_genes, " genes x ", n_samples, " samples")

clamp_k <- if (identical(k_rule, "n_samples_minus_1")) {
  min(n_genes, n_samples) - 1L
} else if (identical(k_rule, "gavish_donoho")) {
  svd_probe <- rsvd::rsvd(norm, k = max(floor((min(n_genes, n_samples) - 1) / 4), 2L))
  eig <- sort(svd_probe$d^2 / (n_samples - 1), decreasing = TRUE)
  as.integer(PCAtools::chooseGavishDonoho(.dim = c(n_genes, n_samples),
                                          var.explained = eig,
                                          noise = median(eig)) * 2)
} else {
  stop("Unknown --k-rule: ", k_rule)
}
svd_k <- clamp_k
message("K = ", clamp_k, " (rule: ", k_rule, ")")

tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

svdres <- rsvd::rsvd(norm, k = svd_k)

message("Running CLAMPbase ...")
base_res <- CLAMP::CLAMPbase(Y = norm, svdres = svdres, trace = TRUE, clamp_k = clamp_k)
base_res$Z <- data.frame(base_res$Z); rownames(base_res$Z) <- genes
base_res$B <- data.frame(base_res$B); colnames(base_res$B) <- samples

# Build the prior: GMT collections plus a CellMarker cell-type set.
gmt_list <- lapply(seq_along(prior_gmts), function(i) {
  sets <- CLAMP:::read_gmt(prior_gmts[[i]])
  names(sets) <- paste0(prior_names[[i]], "_", names(sets))
  sets
})
names(gmt_list) <- prior_names
marker <- data.table::as.data.table(readxl::read_excel(cellmarker_path, sheet = cellmarker_sheet))
if (!all(c(cellmarker_term, cellmarker_gene) %in% names(marker))) {
  stop("CellMarker columns are missing: ", cellmarker_term, ", ", cellmarker_gene)
}
marker <- marker[!is.na(get(cellmarker_term)) & !is.na(get(cellmarker_gene)),
                 .(term = as.character(get(cellmarker_term)), gene = as.character(get(cellmarker_gene)))]
marker_sets <- split(marker$gene, marker$term)
marker_sets <- lapply(marker_sets, unique)
names(marker_sets) <- paste0("CELLMARKER_", make.unique(names(marker_sets)))
gmt_list$CELLMARKER <- marker_sets
path_mat <- CLAMP::gmtListToSparseMat(gmt_list)
prior <- list(path_mat = path_mat,
              matched = CLAMP::getMatchedPathwayMat(path_mat, genes))
message("Prior matched: ", nrow(prior$matched), " genes x ", ncol(prior$matched), " pathways")

message("Running CLAMPfull ...")
full_res <- CLAMP::CLAMPfull(
  Y = norm, svdres = svdres, priorMat = prior$matched,
  clamp.base.result = base_res, use_cpp = TRUE, trace = TRUE,
  max.iter = max_iter, clamp_k = clamp_k
)

model_dir <- file.path(tmp_dir, "CLAMPfull")
write_clamp_model(full_res, model_dir, genes, samples,
                  rds_path = file.path(tmp_dir, "CLAMPfull.rds"))
fwrite(data.frame(L2 = as.numeric(full_res$L2)), file.path(model_dir, "L2.csv"))

n_annot <- if (is.null(full_res$summary)) NA_integer_ else
  sum(full_res$summary$FDR < 0.05 & full_res$summary$AUC > 0.7, na.rm = TRUE)
message("CLAMPfull: ", ncol(full_res$Z), " LVs | pathway annotations (AUC>0.7, FDR<0.05): ", n_annot)

fwrite(data.table(n_genes = n_genes, n_samples = n_samples, clamp_k = clamp_k,
                  svd_k = svd_k, n_lvs = ncol(full_res$Z),
                  n_prior_pathways = ncol(prior$matched),
                  n_annotated_lv_pathway = n_annot),
       file.path(tmp_dir, "fit_summary.csv"))

write_manifest(tmp_dir, method = "CLAMP::CLAMPfull (local)",
               inputs = c(list(norm = norm_path, cellmarker = cellmarker_path),
                          setNames(as.list(prior_gmts), prior_names)),
               parameters = list(k_rule = k_rule, clamp_k = clamp_k, svd_k = svd_k,
                                 max_iter = max_iter, seed = seed,
                                 prior_collections = c(prior_names, "CELLMARKER")),
               extra = list(n_genes = n_genes, n_samples = n_samples,
                            n_lvs = ncol(full_res$Z), l2 = as.numeric(full_res$L2)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
