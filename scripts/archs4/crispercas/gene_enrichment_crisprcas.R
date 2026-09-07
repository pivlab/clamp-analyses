#!/usr/bin/env Rscript
# Score CLAMP LVs against CRISPR-Cas9 lipid-screen gene sets with fgsea.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(fgsea)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    key <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) {
    stop("Missing --", gsub("_", "-", name))
  }
  value
}

args <- parse_cli()
z_path <- required_arg(args, "z")
lipid_deg_path <- required_arg(args, "lipid_deg")
out_path <- required_arg(args, "out")
n_reps <- as.integer(args$n_reps %||% 10L)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

message("Loading Z matrix: ", z_path)
z <- read.csv(z_path, row.names = 1, check.names = FALSE)
z_gene_names <- rownames(z)

lipid_deg <- read.csv(lipid_deg_path, stringsAsFactors = FALSE)

orig_deg_gene_sets <- list()
for (r in unique(lipid_deg$rank)) {
  if (r == 0) next
  orig_deg_gene_sets[[paste0("gene_set_", r)]] <- lipid_deg$gene_name[lipid_deg$rank == r]
}

deg_gene_sets <- list(
  gene_set_increase = orig_deg_gene_sets[["gene_set_3"]],
  gene_set_decrease = orig_deg_gene_sets[["gene_set_-3"]]
)
stopifnot(length(deg_gene_sets[["gene_set_increase"]]) == 6)
stopifnot(length(deg_gene_sets[["gene_set_decrease"]]) == 8)

lvs <- list()
for (cidx in seq_len(ncol(z))) {
  data <- z[, cidx]
  names(data) <- z_gene_names
  lvs[[colnames(z)[cidx]]] <- data
}
message("LVs to score: ", length(lvs))

set.seed(0)
started <- Sys.time()
results <- list()
for (lv in names(lvs)) {
  repetitions <- list()
  for (i in seq_len(n_reps)) {
    rep_res <- fgsea(pathways = deg_gene_sets, stats = lvs[[lv]], scoreType = "pos", eps = 0.0)[order(pval), ]
    rep_res[, "lv"] <- lv
    rep_res[, "rep_idx"] <- i
    repetitions[[i]] <- rep_res
  }
  results[[lv]] <- do.call(rbind, repetitions)
}
df <- do.call(rbind, results)
df <- df %>% mutate(leadingEdge = map_chr(leadingEdge, toString))

write.csv(df, out_path, row.names = FALSE)
message(sprintf(
  "Wrote %s (%d rows, %d LVs) in %.1f min",
  out_path, nrow(df), length(lvs), as.numeric(difftime(Sys.time(), started, units = "mins"))
))
