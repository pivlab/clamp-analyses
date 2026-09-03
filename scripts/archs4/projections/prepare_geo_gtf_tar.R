#!/usr/bin/env Rscript
# Sum per-transcript FluxCapacitor GTF read counts to gene level and emit the
# same counts.rds/norm.csv/samples.csv triple as prepare_geo_matrix.R.

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
raw_dir <- required_arg(args, "raw_dir")
out_dir <- required_arg(args, "out_dir")
annotation <- fromJSON(required_arg(args, "annotation"), simplifyVector = FALSE)
tar_file <- args$tar_file %||% ""
id_type <- args$id_type %||% "ensembl"
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)

started <- Sys.time()

gtf_files <- sort(list.files(raw_dir, pattern = "\\.gtf\\.gz$", full.names = TRUE))
if (!length(gtf_files)) {
  if (!nzchar(tar_file)) stop("No *.gtf.gz in ", raw_dir, " and no --tar-file given")
  message("Extracting ", tar_file)
  untar(tar_file, exdir = raw_dir)
  gtf_files <- sort(list.files(raw_dir, pattern = "\\.gtf\\.gz$", full.names = TRUE))
}
if (!length(gtf_files)) stop("Still no *.gtf.gz after extracting ", tar_file)
message("Parsing ", length(gtf_files), " GTF files")

tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

parse_flux_gtf <- function(path) {
  raw <- fread(cmd = paste("zcat", shQuote(path)), sep = "\t", header = FALSE, quote = "",
               col.names = c("seqname", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attributes"))
  raw <- raw[feature == "transcript"]
  raw[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
  raw[, reads_val := as.numeric(sub(".*; reads ([0-9.]+);.*", "\\1", attributes))]
  raw[is.na(reads_val), reads_val := 0]
  raw[, gene_id_clean := sub("\\..*", "", gene_id)]
  raw[, .(gene_reads = sum(reads_val, na.rm = TRUE)), by = gene_id_clean]
}

samples <- apply_annotation_rules(gtf_files, annotation)

expr_list <- lapply(seq_along(gtf_files), function(i) {
  dt <- parse_flux_gtf(gtf_files[i])
  setnames(dt, "gene_reads", samples$sample[i])
  dt
})
merged <- Reduce(function(a, b) merge(a, b, by = "gene_id_clean", all = TRUE), expr_list)
for (col in names(merged)[-1]) set(merged, which(is.na(merged[[col]])), col, 0)

counts <- as.matrix(merged[, -1L, with = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- merged$gene_id_clean
n_raw <- nrow(counts)

if (identical(id_type, "ensembl")) counts <- ensembl_to_symbol(counts)
norm <- normalize_counts(counts, mean_cutoff, var_cutoff)
message("Prepared: ", nrow(norm), " genes x ", ncol(norm), " samples")

saveRDS(counts, file.path(tmp_dir, "counts.rds"))
write_matrix_csv(as.matrix(norm), file.path(tmp_dir, "norm.csv"))
fwrite(samples, file.path(tmp_dir, "samples.csv"))
fwrite(data.table(n_genes_raw = n_raw, n_genes_mapped = nrow(counts),
                  n_genes_normalized = nrow(norm), n_samples = ncol(norm),
                  n_gtf_files = length(gtf_files)),
       file.path(tmp_dir, "prepare_summary.csv"))

write_manifest(tmp_dir, method = "prepare_geo_gtf_tar",
               inputs = if (nzchar(tar_file)) list(tar = tar_file) else list(),
               parameters = list(id_type = id_type, mean_cutoff = mean_cutoff,
                                 var_cutoff = var_cutoff),
               extra = list(n_gtf_files = length(gtf_files), n_genes_raw = n_raw,
                            n_genes_mapped = nrow(counts),
                            n_genes_normalized = nrow(norm), n_samples = ncol(norm)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
