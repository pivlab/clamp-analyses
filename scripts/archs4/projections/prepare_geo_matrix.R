#!/usr/bin/env Rscript
# Turn a GEO processed count matrix into counts.rds/norm.csv/samples.csv, the
# canonical triple every downstream projection rule consumes.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(CLAMP)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
counts_file <- required_arg(args, "counts_file")
out_dir <- required_arg(args, "out_dir")
annotation <- fromJSON(required_arg(args, "annotation"), simplifyVector = FALSE)
id_type <- args$id_type %||% "ensembl"
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)

started <- Sys.time()
tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

dt <- fread(counts_file)
counts <- as.matrix(dt[, -1L, with = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- as.character(dt[[1L]])
n_raw <- nrow(counts)

drop_prefix <- annotation$drop_row_prefix %||% ""
if (nzchar(drop_prefix)) {
  counts <- counts[!startsWith(rownames(counts), drop_prefix), , drop = FALSE]
}

if (identical(id_type, "ensembl")) counts <- ensembl_to_symbol(counts)
if (anyDuplicated(rownames(counts))) counts <- rowsum(counts, rownames(counts))

samples <- apply_annotation_rules(colnames(counts), annotation)
colnames(counts) <- samples$sample

norm <- normalize_counts(counts, mean_cutoff, var_cutoff)
message("Prepared: ", nrow(norm), " genes x ", ncol(norm), " samples")

saveRDS(counts, file.path(tmp_dir, "counts.rds"))
write_matrix_csv(as.matrix(norm), file.path(tmp_dir, "norm.csv"))
fwrite(samples, file.path(tmp_dir, "samples.csv"))
fwrite(data.table(n_genes_raw = n_raw, n_genes_mapped = nrow(counts),
                  n_genes_normalized = nrow(norm), n_samples = ncol(norm)),
       file.path(tmp_dir, "prepare_summary.csv"))

write_manifest(tmp_dir, method = "prepare_geo_matrix",
               inputs = list(counts = counts_file),
               parameters = list(id_type = id_type, mean_cutoff = mean_cutoff,
                                 var_cutoff = var_cutoff,
                                 drop_row_prefix = drop_prefix),
               extra = list(n_genes_raw = n_raw, n_genes_mapped = nrow(counts),
                            n_genes_normalized = nrow(norm), n_samples = ncol(norm)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
