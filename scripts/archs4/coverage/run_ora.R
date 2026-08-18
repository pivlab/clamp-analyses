#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(data.table)
  library(jsonlite)
  library(readxl)
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

read_gmt <- function(path) {
  records <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  valid <- lengths(records) >= 3L
  records <- records[valid]
  unique(data.table(
    term = rep(vapply(records, `[[`, character(1), 1L), lengths(records) - 2L),
    gene = unlist(lapply(records, function(record) record[-c(1L, 2L)]), use.names = FALSE)
  ))
}

args <- parse_cli()
z_path <- required_arg(args, "z")
out_dir <- required_arg(args, "out_dir")
dataset <- required_arg(args, "dataset")
fraction <- as.integer(args$fraction %||% 100L)
seed <- as.integer(args$seed %||% 1L)
model <- required_arg(args, "model")
if (!is.null(args$n_samples)) {
  n_samples <- as.integer(args$n_samples)
} else if (!is.null(args$model_manifest)) {
  n_samples <- as.integer(jsonlite::read_json(args$model_manifest)$samples)
} else if (!is.null(args$subsample_info)) {
  n_samples <- as.integer(readRDS(args$subsample_info)$n_samples)
} else if (!is.null(args$samples_rds)) {
  n_samples <- length(readRDS(args$samples_rds))
} else {
  stop("Provide --n-samples, --model-manifest, --subsample-info, or --samples-rds")
}
database <- required_arg(args, "database")
database_label <- args$database_label %||% database
database_type <- required_arg(args, "database_type")
database_path <- required_arg(args, "database_path")
top_pct <- as.numeric(args$top_pct %||% 0.01)
min_size <- as.integer(args$min_size %||% 10L)
max_size <- as.integer(args$max_size %||% 500L)
pvalue_cutoff <- as.numeric(args$pvalue_cutoff %||% 0.05)
qvalue_cutoff <- as.numeric(args$qvalue_cutoff %||% 0.2)
fdr <- as.numeric(args$fdr %||% 0.05)

if (dir.exists(out_dir)) stop("Refusing to overwrite ORA directory: ", out_dir)
dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)
tmp_dir <- paste0(out_dir, ".tmp.", Sys.getpid())
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
dir.create(tmp_dir, recursive = TRUE)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

if (database_type == "gmt") {
  term2gene <- read_gmt(database_path)
  exclude_regex <- args$exclude_term_regex %||% ""
  if (nzchar(exclude_regex)) term2gene <- term2gene[!grepl(exclude_regex, term)]
} else if (database_type == "cellmarker_xlsx") {
  sheet <- args$sheet %||% "human"
  term_column <- args$term_column %||% "cell_name"
  gene_column <- args$gene_column %||% "Symbol"
  marker <- as.data.table(readxl::read_excel(database_path, sheet = sheet))
  if (!all(c(term_column, gene_column) %in% names(marker))) {
    stop("CellMarker columns are missing: ", term_column, ", ", gene_column)
  }
  term2gene <- unique(marker[
    !is.na(get(term_column)) & !is.na(get(gene_column)),
    .(term = as.character(get(term_column)), gene = as.character(get(gene_column)))
  ])
} else {
  stop("Unsupported database type: ", database_type)
}

z <- fread(z_path, check.names = FALSE)
if (ncol(z) < 2L) stop("Z must contain a gene column and at least one LV")
genes <- as.character(z[[1L]])
if (anyDuplicated(genes)) stop("Z contains duplicate genes")
loadings <- as.matrix(z[, -1L])
storage.mode(loadings) <- "double"
rownames(loadings) <- genes
lv_names <- colnames(loadings)
universe <- genes
query_target <- max(1L, ceiling(length(universe) * top_pct))

eligible <- term2gene[gene %chin% universe, .(set_size = uniqueN(gene)), by = term]
eligible <- eligible[set_size >= min_size & set_size <= max_size]
eligible_terms <- eligible$term
eligible_term2gene <- term2gene[term %chin% eligible_terms & gene %chin% universe]
effective_universe <- unique(eligible_term2gene$gene)
if (!length(eligible_terms)) stop("No eligible terms after universe and size filtering")

writeLines(universe, file.path(tmp_dir, "universe.txt"))
fwrite(eligible, file.path(tmp_dir, "eligible_terms.csv"))

results_path <- file.path(tmp_dir, "enrichment.csv.gz")
connection <- gzfile(results_path, open = "wt")
on.exit(if (isOpen(connection)) close(connection), add = TRUE)
wrote_header <- FALSE
result_rows <- 0L
tested_lvs <- 0L
recovered <- character()
started <- Sys.time()

for (lv_index in seq_len(ncol(loadings))) {
  values <- loadings[, lv_index]
  positive <- which(is.finite(values) & values > 0)
  if (!length(positive)) next
  take <- min(query_target, length(positive))
  selected <- positive[order(values[positive], decreasing = TRUE)[seq_len(take)]]
  query <- genes[selected]
  if (!all(query %chin% universe)) stop("ORA query escaped its explicit universe")

  enrichment <- suppressMessages(clusterProfiler::enricher(
    gene = query,
    universe = universe,
    TERM2GENE = eligible_term2gene,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    minGSSize = min_size,
    maxGSSize = max_size,
    qvalueCutoff = qvalue_cutoff
  ))
  tested_lvs <- tested_lvs + 1L
  if (is.null(enrichment)) next
  frame <- as.data.table(as.data.frame(enrichment))
  if (!nrow(frame)) next
  frame[, `:=`(
    LV = lv_names[[lv_index]],
    LV_index = lv_index,
    query_size = length(query)
  )]
  setcolorder(frame, c("LV", "LV_index", "query_size", setdiff(names(frame), c("LV", "LV_index", "query_size"))))
  write.table(
    frame,
    connection,
    sep = ",",
    row.names = FALSE,
    col.names = !wrote_header,
    append = wrote_header,
    quote = TRUE,
    na = ""
  )
  wrote_header <- TRUE
  result_rows <- result_rows + nrow(frame)
  recovered <- union(recovered, frame[p.adjust < fdr, ID])
}
close(connection)

if (!wrote_header) {
  empty <- data.table(
    LV = character(), LV_index = integer(), query_size = integer(),
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer()
  )
  fwrite(empty, results_path, compress = "gzip")
}

finished <- Sys.time()
summary <- data.table(
  dataset = dataset,
  fraction = fraction,
  seed = seed,
  model = model,
  n_samples = n_samples,
  database = database,
  database_label = database_label,
  universe_size = length(universe),
  effective_annotated_universe_size = length(effective_universe),
  query_target = query_target,
  latent_variables = ncol(loadings),
  tested_lvs = tested_lvs,
  eligible_pathways = length(eligible_terms),
  recovered_pathways = length(recovered),
  recovered_percent = 100 * length(recovered) / length(eligible_terms),
  fdr = fdr,
  result_rows = result_rows,
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs"))
)
fwrite(summary, file.path(tmp_dir, "summary.csv"))
writeLines(sort(recovered), file.path(tmp_dir, "recovered_pathways.txt"))

manifest <- list(
  schema_version = 1L,
  status = "complete",
  method = "clusterProfiler::enricher",
  backend = if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "slurm" else "local",
  dataset = dataset,
  fraction = fraction,
  seed = seed,
  model = model,
  n_samples = n_samples,
  database = database,
  database_label = database_label,
  database_path = normalizePath(database_path),
  database_md5 = unname(tools::md5sum(database_path)),
  z_path = normalizePath(z_path),
  z_md5 = unname(tools::md5sum(z_path)),
  universe_size = length(universe),
  effective_annotated_universe_size = length(effective_universe),
  query_target = query_target,
  eligible_pathways = length(eligible_terms),
  recovered_pathways = length(recovered),
  parameters = list(
    top_pct = top_pct,
    minGSSize = min_size,
    maxGSSize = max_size,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    recovery_fdr = fdr,
    across_lv_correction = FALSE
  ),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  package_versions = list(
    R = as.character(getRversion()),
    clusterProfiler = as.character(packageVersion("clusterProfiler"))
  )
)
jsonlite::write_json(manifest, file.path(tmp_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("complete", file.path(tmp_dir, "complete"))

if (!file.rename(tmp_dir, out_dir)) stop("Could not atomically publish ORA directory")
published <- TRUE
message("Published ORA: ", out_dir)
