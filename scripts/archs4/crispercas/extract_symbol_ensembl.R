#!/usr/bin/env Rscript
# Map gene symbols to Ensembl identifiers for gene-level TWAS lookups.

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: extract_symbol_ensembl.R OUTPUT_CSV")
}

out <- args[[1L]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  columns = "ENSEMBL",
  keytype = "SYMBOL"
)

write.csv(mapping, out, row.names = FALSE)
