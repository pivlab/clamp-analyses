#!/usr/bin/env Rscript
# Score how much of a dataset's known biology each model recovers, from the
# full member list of each mechanism's best-matching ORA gene set.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
models <- fread(required_arg(args, "models"))
ora <- fread(required_arg(args, "ora"))
mechanisms <- fromJSON(required_arg(args, "mechanisms"), simplifyVector = FALSE)
out_dir <- required_arg(args, "out_dir")
top_pct <- as.numeric(args$top_pct %||% 0.01)
min_count <- as.integer(args$min_count %||% 3L)
gene_set_gmt <- required_arg(args, "gene_set_gmt")
hallmark_gmt <- required_arg(args, "hallmark_gmt")
cellmarker_file <- required_arg(args, "cellmarker_file")
cellmarker_sheet <- args$cellmarker_sheet %||% "human"
cellmarker_term_column <- args$cellmarker_term_column %||% "cell_name"
cellmarker_gene_column <- args$cellmarker_gene_column %||% "Symbol"
fdr <- as.numeric(args$fdr %||% 0.05)
lv_fdr <- as.numeric(args$lv_fdr %||% 0.05)
lv_logfc <- as.numeric(args$lv_logfc %||% 0.05)

read_gmt_sets <- function(path) {
  x <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(x, function(row) unique(row[-c(1L, 2L)]))
  names(out) <- vapply(x, `[[`, "", 1L)
  out
}
read_cellmarker_sets <- function(path, sheet, term_col, gene_col) {
  x <- data.table::as.data.table(readxl::read_excel(path, sheet = sheet))
  if (!all(c(term_col, gene_col) %in% names(x))) {
    stop("CellMarker columns are missing: ", term_col, ", ", gene_col)
  }
  x <- x[!is.na(get(term_col)) & !is.na(get(gene_col)),
         .(term = as.character(get(term_col)), gene = as.character(get(gene_col)))]
  split(x$gene, x$term)
}

gene_sets <- list(canonical = read_gmt_sets(gene_set_gmt),
                  hallmark = read_gmt_sets(hallmark_gmt),
                  cellmarker = read_cellmarker_sets(cellmarker_file, cellmarker_sheet,
                                                     cellmarker_term_column, cellmarker_gene_column))

members_for <- function(database, term) {
  unique(as.character(gene_sets[[database]][[term]] %||% character()))
}

started <- Sys.time()
tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

hits <- rbindlist(lapply(seq_len(nrow(ora)), function(i) {
  f <- file.path(ora$dir[i], "enrichment.csv.gz")
  if (!file.exists(f)) {
    warning("Missing ORA output: ", f); return(NULL)
  }
  d <- fread(f)
  if (!nrow(d)) return(NULL)
  d[, `:=`(
    model = ora$model[i], database = ora$database[i],
    query_n = as.integer(tstrsplit(GeneRatio, "/")[[2]]),
    set_size = as.integer(tstrsplit(BgRatio, "/")[[1]]),
    bg_n = as.integer(tstrsplit(BgRatio, "/")[[2]]))]
  d[, fold := (Count / query_n) / (set_size / bg_n)]
  d[]
}), fill = TRUE)
if (is.null(hits) || !nrow(hits)) stop("No ORA results found")
setnames(hits, "ID", "term", skip_absent = TRUE)
hits <- hits[p.adjust < fdr & Count >= min_count]
fwrite(hits[order(model, database, LV, p.adjust)], file.path(tmp_dir, "lv_pathway_hits.csv"))

significant_for <- function(model_row, spec, z) {
  if (is.null(spec$contrast) || !nzchar(spec$contrast)) {
    sig <- readLines(model_row$significant_lvs, warn = FALSE)
    return(intersect(sig[nzchar(sig)], colnames(z)))
  }
  stats <- fread(model_row$lv_stats)
  intersect(unique(stats[contrast == spec$contrast & adj.P.Val < lv_fdr &
                         abs_logFC >= lv_logfc, LV]), colnames(z))
}

model_cache <- lapply(seq_len(nrow(models)), function(i) {
  list(model = models$model[i], row = models[i], z = read_matrix_csv(models$z[i]))
})
names(model_cache) <- models$model

candidate_hits <- rbindlist(lapply(model_cache, function(entry) {
  rbindlist(lapply(names(mechanisms), function(nm) {
    spec <- mechanisms[[nm]]
    category <- spec$category %||% "molecular_mechanism"
    sig <- significant_for(entry$row, spec, entry$z)
    pat <- paste(unlist(spec$terms %||% character()), collapse = "|")
    if (!length(sig) || !nzchar(pat)) return(NULL)
    d <- hits[model == entry$model & LV %chin% sig &
                grepl(pat, term, ignore.case = TRUE)][order(p.adjust)]
    if (identical(category, "molecular_mechanism")) {
      d <- d[database == (spec$database %||% "canonical")]
    }
    if (!nrow(d)) return(NULL)
    d[, `:=`(mechanism = nm,
             testable = isTRUE(spec$testable %||% TRUE),
             category = category,
             target_contrast = spec$contrast %||% NA_character_)]
    d
  }), fill = TRUE)
}), fill = TRUE)

evidence <- rbindlist(lapply(model_cache, function(entry) {
  rbindlist(lapply(names(mechanisms), function(nm) {
    spec <- mechanisms[[nm]]
    testable <- isTRUE(spec$testable %||% TRUE)
    sig <- significant_for(entry$row, spec, entry$z)
    pat <- paste(unlist(spec$terms %||% character()), collapse = "|")
    ph <- candidate_hits[model == entry$model & mechanism == nm][order(p.adjust)]
    best <- if (nrow(ph)) ph[1] else NULL
    members <- if (nrow(ph)) intersect(members_for(best$database, best$term), rownames(entry$z)) else character()
    n_top <- max(1L, ceiling(nrow(entry$z) * top_pct))
    top_genes <- if (nrow(ph)) rownames(entry$z)[order(entry$z[, best$LV], decreasing = TRUE)[seq_len(n_top)]] else character()
    data.table(
      model = entry$model, mechanism = nm,
      category = spec$category %||% "molecular_mechanism",
      source = spec$source %||% NA_character_, doi = spec$doi %||% NA_character_,
      target_contrast = spec$contrast %||% NA_character_, testable = testable,
      n_target_lvs = length(sig), n_lvs = uniqueN(ph$LV), n_terms = nrow(ph),
      best_lv = if (nrow(ph)) best$LV else NA_character_,
      best_term = if (nrow(ph)) best$term else NA_character_,
      best_database = if (nrow(ph)) best$database else NA_character_,
      best_q = if (nrow(ph)) best$p.adjust else NA_real_,
      best_fold = if (nrow(ph)) best$fold else NA_real_,
      best_gene_ratio = if (nrow(ph)) best$GeneRatio else NA_character_,
      n_gene_set_in_universe = length(members),
      n_gene_set_in_top_loading = sum(members %chin% top_genes),
      gene_set_fraction_in_top_loading = if (length(members)) sum(members %chin% top_genes) / length(members) else NA_real_,
      pathway_evidence = if (testable) nrow(ph) > 0 else NA,
      recovered = if (testable) nrow(ph) > 0 else NA)
  }), fill = TRUE)
}), fill = TRUE)
fwrite(evidence, file.path(tmp_dir, "mechanism_recovery.csv"))

if (!nrow(candidate_hits[testable == TRUE])) {
  empty_comparisons <- data.table(
    comparison_id = character(), mechanism = character(), category = character(),
    database = character(), gene_set = character(), model = character(), LV = character(),
    recovered = logical(), target_contrast = character(), ora_fdr = numeric(), n_gene_set_in_universe = integer(),
    n_gene_set_in_top_loading = integer(), top_pct = numeric(), p_value = numeric(),
    arch_fraction = numeric(), local_fraction = numeric(), p_adj = numeric(),
    comparison_note = character())
  empty_loadings <- data.table(
    comparison_id = character(), model = character(), mechanism = character(), category = character(),
    gene_set = character(), database = character(), LV = character(), rank = integer(),
    gene = character(), loading = numeric(), is_gene_set = logical(),
    in_top_loading_set = logical(), top_pct = numeric())
  fwrite(empty_comparisons, file.path(tmp_dir, "gene_set_comparisons.csv"))
  fwrite(empty_loadings, file.path(tmp_dir, "gene_set_loading_points.csv"))
  fwrite(data.table(model = character(), mechanism = character(), gene_set = character(),
                    database = character(), LV = character(), gene = character(), rank = integer(),
                    in_top_loading_set = logical()),
         file.path(tmp_dir, "gene_set_assignment.csv"))
  write_manifest(tmp_dir, method = "mechanism recovery (ORA gene sets)",
                 inputs = list(models = required_arg(args, "models"), ora = required_arg(args, "ora")),
                 parameters = list(top_pct = top_pct, min_count = min_count, fdr = fdr,
                                   lv_fdr = lv_fdr, lv_logfc = lv_logfc,
                                   gene_set_gmt = gene_set_gmt, hallmark_gmt = hallmark_gmt,
                                   cellmarker_file = cellmarker_file,
                                   mechanisms = names(mechanisms)),
                 extra = list(n_models = nrow(models), n_mechanisms = length(mechanisms),
                              n_pathway_hits = nrow(hits)), started = started)
  finish_publish(tmp_dir, out_dir)
  published <- TRUE
  quit(save = "no", status = 0)
}

target_defs <- unique(candidate_hits[testable == TRUE,
  .(mechanism, category, database, target_contrast)])
target_defs[, comparison_id := paste(mechanism, database, sep = "::")]

comparison_rows <- rbindlist(lapply(seq_len(nrow(target_defs)), function(i) {
  target <- target_defs[i]
  rbindlist(lapply(model_cache, function(entry) {
    ph <- candidate_hits[model == entry$model & mechanism == target$mechanism &
                           database == target$database][order(p.adjust)]
    recovered <- nrow(ph) > 0L
    best <- if (recovered) ph[1] else NULL
    members <- if (recovered) intersect(members_for(target$database, best$term), rownames(entry$z)) else character()
    n_top <- max(1L, ceiling(nrow(entry$z) * top_pct))
    top_genes <- if (recovered) rownames(entry$z)[order(entry$z[, best$LV], decreasing = TRUE)[seq_len(n_top)]] else character()
    data.table(comparison_id = target$comparison_id, mechanism = target$mechanism,
               category = target$category, database = target$database,
               gene_set = if (recovered) best$term else NA_character_,
               model = entry$model, LV = if (recovered) best$LV else NA_character_,
               recovered = recovered, target_contrast = target$target_contrast,
               ora_fdr = if (recovered) best$p.adjust else NA_real_,
               n_gene_set_in_universe = length(members),
               n_gene_set_in_top_loading = sum(members %chin% top_genes),
               top_pct = top_pct)
  }), fill = TRUE)
}), fill = TRUE)

comparison_tests <- comparison_rows[, {
  a <- .SD[model == "ARCHS4"]
  b <- .SD[model == "local"]
  p <- NA_real_
  arch_fraction <- NA_real_
  local_fraction <- NA_real_
  if (nrow(a) == 1L && nrow(b) == 1L && isTRUE(a$recovered) && isTRUE(b$recovered) &&
      identical(a$gene_set, b$gene_set) &&
      a$n_gene_set_in_universe > 0L && b$n_gene_set_in_universe > 0L) {
    arch_fraction <- a$n_gene_set_in_top_loading / a$n_gene_set_in_universe
    local_fraction <- b$n_gene_set_in_top_loading / b$n_gene_set_in_universe
    p <- fisher.test(matrix(c(a$n_gene_set_in_top_loading,
                              a$n_gene_set_in_universe - a$n_gene_set_in_top_loading,
                              b$n_gene_set_in_top_loading,
                              b$n_gene_set_in_universe - b$n_gene_set_in_top_loading),
                            nrow = 2, byrow = TRUE), alternative = "greater")$p.value
  }
  .(p_value = p, arch_fraction = arch_fraction, local_fraction = local_fraction)
}, by = comparison_id]
comparison_tests[, p_adj := p.adjust(p_value, method = "BH")]
comparison_rows <- merge(comparison_rows, comparison_tests, by = "comparison_id", all.x = TRUE)
comparison_rows[, comparison_note := ""]
comparison_rows[model == "local" & !recovered, comparison_note := "Not recovered"]
comparison_rows[recovered == TRUE & is.na(p_adj), comparison_note := "Different matched pathway"]
comparison_rows[recovered == TRUE & !is.na(p_adj), comparison_note := paste0(
  fifelse(arch_fraction > local_fraction, "ARCHS4 > local",
          fifelse(arch_fraction < local_fraction, "local > ARCHS4", "equal fraction")),
  " (FDR ", formatC(p_adj, format = "e", digits = 1), ")")]
fwrite(comparison_rows, file.path(tmp_dir, "gene_set_comparisons.csv"))

loading_points <- rbindlist(lapply(model_cache, function(entry) {
  rbindlist(lapply(seq_len(nrow(comparison_rows[model == entry$model & recovered == TRUE])), function(i) {
    e <- comparison_rows[model == entry$model & recovered == TRUE][i]
    members <- intersect(members_for(e$database, e$gene_set), rownames(entry$z))
    ord <- order(entry$z[, e$LV], decreasing = TRUE)
    data.table(comparison_id = e$comparison_id, model = entry$model,
               mechanism = e$mechanism, category = e$category, gene_set = e$gene_set,
               database = e$database, LV = e$LV, rank = seq_along(ord),
               gene = rownames(entry$z)[ord], loading = entry$z[ord, e$LV],
               is_gene_set = rownames(entry$z)[ord] %chin% members,
               in_top_loading_set = seq_along(ord) <= ceiling(nrow(entry$z) * top_pct),
               top_pct = top_pct)
  }), fill = TRUE)
}), fill = TRUE)
if (!nrow(loading_points)) loading_points <- data.table(
  comparison_id = character(), model = character(), mechanism = character(), category = character(),
  gene_set = character(), database = character(), LV = character(), rank = integer(),
  gene = character(), loading = numeric(), is_gene_set = logical(),
  in_top_loading_set = logical(), top_pct = numeric())
fwrite(loading_points, file.path(tmp_dir, "gene_set_loading_points.csv"))

for (mdl in models$model) {
  s <- evidence[model == mdl]
  tested <- s[testable == TRUE]
  message(mdl, ": ", sum(tested$recovered, na.rm = TRUE), "/", nrow(tested),
          " testable targets recovered",
          if (any(tested$recovered == FALSE, na.rm = TRUE))
            paste0(" | missing: ", paste(tested[recovered == FALSE]$mechanism, collapse = ", ")) else "")
}

write_manifest(tmp_dir, method = "mechanism recovery (ORA gene sets)",
               inputs = list(models = required_arg(args, "models"), ora = required_arg(args, "ora")),
               parameters = list(top_pct = top_pct, min_count = min_count, fdr = fdr,
                                 lv_fdr = lv_fdr, lv_logfc = lv_logfc,
                                 gene_set_gmt = gene_set_gmt, hallmark_gmt = hallmark_gmt,
                                 cellmarker_file = cellmarker_file,
                                 mechanisms = names(mechanisms)),
               extra = list(n_models = nrow(models), n_mechanisms = length(mechanisms),
                            n_pathway_hits = nrow(hits)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
