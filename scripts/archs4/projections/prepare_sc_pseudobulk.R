#!/usr/bin/env Rscript
# Pseudobulk a single-cell dataset (donor x cell type) and emit the same
# canonical triple as the bulk prep scripts.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(Matrix)
  library(CLAMP)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
raw_dir <- required_arg(args, "raw_dir")
kind <- required_arg(args, "kind")
spec <- fromJSON(required_arg(args, "spec"), simplifyVector = FALSE)
out_dir <- required_arg(args, "out_dir")
mean_cutoff <- as.numeric(args$mean_cutoff %||% 0.5)
var_cutoff <- as.numeric(args$var_cutoff %||% 0.1)

started <- Sys.time()
pb <- spec$pseudobulk %||% list()
min_cells <- as.integer(pb$min_cells %||% 20L)

map_cell_type <- function(x, rules, fallback = "Other") {
  out <- rep(fallback, length(x))
  lc <- tolower(trimws(x))
  for (r in rev(rules)) out[grepl(r$pattern, lc, ignore.case = TRUE)] <- r$label
  out[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  out
}

pool_cells <- function(mat, group) {
  keep <- !is.na(group)
  mat <- mat[, keep, drop = FALSE]; group <- factor(group[keep])
  tab <- table(group)
  ok <- names(tab)[tab >= min_cells]
  sel <- group %in% ok
  mat <- mat[, sel, drop = FALSE]; group <- factor(as.character(group[sel]), levels = ok)
  ind <- Matrix::sparseMatrix(i = seq_along(group), j = as.integer(group), x = 1,
                              dims = c(length(group), length(ok)))
  counts <- as.matrix(mat %*% ind)
  colnames(counts) <- ok
  list(counts = counts, n_cells = as.integer(tab[ok]))
}

if (identical(kind, "ebi_sc_atlas")) {
  mc <- spec$metadata_columns
  meta_raw <- fread(file.path(raw_dir, spec$metadata_file), sep = "\t")
  if (!is.null(spec$filter)) {
    meta_raw <- meta_raw[tolower(get(spec$filter$column)) == tolower(spec$filter$equals)]
  }
  meta <- data.table(cell = as.character(meta_raw[[mc$id]]),
                     donor = as.character(meta_raw[[mc$donor]]),
                     raw_label = as.character(meta_raw[[mc$label]]))
  base <- sub("\\.cell_metadata\\.tsv$", "", spec$metadata_file)
  genes <- fread(file.path(raw_dir, paste0(base, ".aggregated_filtered_counts.mtx_rows.gz")),
                 header = FALSE, sep = "\t")[[1]]
  cells <- fread(file.path(raw_dir, paste0(base, ".aggregated_filtered_counts.mtx_cols.gz")),
                 header = FALSE, sep = "\t")[[1]]
  message("Reading sparse matrix ...")
  mat <- Matrix::readMM(gzcon(file(file.path(raw_dir, paste0(base, ".aggregated_filtered_counts.mtx.gz")), open = "rb")))
  rownames(mat) <- genes; colnames(mat) <- cells
  keep <- intersect(colnames(mat), meta$cell)
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[match(keep, cell)]
} else {
  stop("Unknown --kind: ", kind)
}

annotation_method <- spec$annotation$method %||% "authors_metadata"
annotation_source <- spec$annotation$source %||% "source metadata"
meta[, cell_type := map_cell_type(raw_label, spec$cell_type_map)]
meta[, `:=`(annotation_method = annotation_method, annotation_source = annotation_source)]
meta[, sample := paste0(donor, "__", gsub("[^A-Za-z0-9]+", "_", cell_type))]
message("Typed cells: ", sum(!is.na(meta$cell_type)), " of ", nrow(meta))

tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

pooled <- pool_cells(mat, ifelse(is.na(meta$cell_type), NA_character_, meta$sample))
counts <- pooled$counts
message("Pseudobulk: ", nrow(counts), " genes x ", ncol(counts), " samples (>= ", min_cells, " cells)")

samples <- unique(meta[sample %chin% colnames(counts), .(sample, donor, cell_type)])
samples <- samples[match(colnames(counts), sample)]
samples[, n_cells := pooled$n_cells]
cf <- spec$condition_from
if (!is.null(cf)) {
  for (nm in names(cf$rules)) {
    r <- cf$rules[[nm]]
    samples[[nm]] <- ifelse(grepl(r$pattern, samples[[cf$column]], ignore.case = TRUE),
                            r$match, r$nomatch)
  }
}
for (ct in spec$cell_type_contrasts %||% list()) {
  samples[[ct$column]] <- ifelse(samples$cell_type == ct$label, ct$label, "Other")
}

if (identical(spec$id_type, "ensembl")) counts <- ensembl_to_symbol(counts)
if (anyDuplicated(rownames(counts))) counts <- rowsum(counts, rownames(counts))
norm <- normalize_counts(counts, mean_cutoff, var_cutoff)
message("Prepared: ", nrow(norm), " genes x ", ncol(norm), " samples")

saveRDS(counts, file.path(tmp_dir, "counts.rds"))
write_matrix_csv(as.matrix(norm), file.path(tmp_dir, "norm.csv"))
fwrite(samples, file.path(tmp_dir, "samples.csv"))
fwrite(meta, file.path(tmp_dir, "cell_metadata.csv"))
fwrite(data.table(n_cells_total = nrow(meta), n_cells_typed = sum(!is.na(meta$cell_type)),
                  n_pseudobulk_samples = ncol(counts), n_genes_mapped = nrow(counts),
                  n_genes_normalized = nrow(norm), annotation_method = annotation_method),
       file.path(tmp_dir, "prepare_summary.csv"))

write_manifest(tmp_dir, method = paste0("prepare_sc_pseudobulk (", kind, ")"),
               parameters = list(kind = kind, min_cells = min_cells,
                                 mean_cutoff = mean_cutoff, var_cutoff = var_cutoff,
                                 annotation_method = annotation_method,
                                 annotation_source = annotation_source),
               extra = list(n_cells_total = nrow(meta),
                            n_cells_typed = sum(!is.na(meta$cell_type)),
                            n_pseudobulk_samples = ncol(counts),
                            n_genes_normalized = nrow(norm)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
