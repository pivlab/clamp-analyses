#!/usr/bin/env Rscript
# Gather the per-dataset outputs into the tables the report notebooks read.
# Notebooks do no computation: everything they plot is produced here.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
prod_root <- required_arg(args, "prod_root")
datasets <- strsplit(required_arg(args, "datasets"), ",")[[1]]
groups <- strsplit(required_arg(args, "groups"), ",")[[1]]
out_dir <- required_arg(args, "out_dir")
stopifnot(length(datasets) == length(groups))
names(groups) <- datasets
models <- c("ARCHS4", "local")

started <- Sys.time()
tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

gather <- function(rel, per_model = TRUE) {
  rows <- list()
  for (d in datasets) {
    targets <- if (per_model) models else NA_character_
    for (m in targets) {
      f <- file.path(prod_root, d, if (is.na(m)) rel else sprintf(rel, m))
      if (!file.exists(f)) next
      x <- fread(f)
      if (!nrow(x)) next
      x[, dataset := d][, group := groups[[d]]]
      if (!is.na(m) && !"model" %in% names(x)) x[, model := m]
      rows[[length(rows) + 1L]] <- x
    }
  }
  if (!length(rows)) return(data.table())
  rbindlist(rows, fill = TRUE)
}

lv_stats <- gather("lv_stats/%s/lv_stats.csv")
sig_lvs <- gather("lv_stats/%s/significant_lvs.csv")
mech <- gather("mechanism_recovery/mechanism_recovery.csv", per_model = FALSE)
gene_set_members <- gather("mechanism_recovery/gene_set_assignment.csv", per_model = FALSE)
gene_set_comparisons <- gather("mechanism_recovery/gene_set_comparisons.csv", per_model = FALSE)
gene_set_loadings <- gather("mechanism_recovery/gene_set_loading_points.csv", per_model = FALSE)
hits <- gather("mechanism_recovery/lv_pathway_hits.csv", per_model = FALSE)

fwrite(lv_stats, file.path(tmp_dir, "lv_stats_long.csv"))
fwrite(sig_lvs, file.path(tmp_dir, "significant_lvs_long.csv"))
fwrite(mech, file.path(tmp_dir, "mechanism_recovery_long.csv"))
fwrite(gene_set_members, file.path(tmp_dir, "gene_set_assignment_long.csv"))
fwrite(gene_set_comparisons, file.path(tmp_dir, "gene_set_comparisons_long.csv"))
fwrite(gene_set_loadings, file.path(tmp_dir, "gene_set_loading_points_long.csv"))
fwrite(hits, file.path(tmp_dir, "lv_pathway_hits_long.csv"))

if (nrow(mech)) {
  fwrite(dcast(mech, dataset + group + mechanism ~ model, value.var = "recovered"),
         file.path(tmp_dir, "mechanism_recovery_wide.csv"))
  score <- mech[testable == TRUE,
                .(n_mechanisms = .N, n_recovered = sum(recovered, na.rm = TRUE)),
                by = .(dataset, group, model)]
  fwrite(dcast(score, dataset + group ~ model,
               value.var = c("n_mechanisms", "n_recovered")),
         file.path(tmp_dir, "mechanism_score_by_model.csv"))
}
if (nrow(sig_lvs) && nrow(hits)) {
  best <- hits[order(dataset, model, LV, p.adjust),
               .SD[1], by = .(dataset, model, LV),
               .SDcols = c("term", "database", "GeneRatio", "p.adjust", "fold")]
  setnames(best, c("term", "p.adjust"), c("best_term", "best_q"))
  panel <- merge(sig_lvs, best, by = c("dataset", "model", "LV"), all.x = TRUE)
  fwrite(panel[order(dataset, model, -max_abs_logFC)], file.path(tmp_dir, "panel_ready.csv"))
}

summary_dt <- data.table(dataset = datasets, group = unname(groups[datasets]))
summary_dt <- merge(summary_dt,
  sig_lvs[, .(n_significant_lvs = .N), by = .(dataset, model)], by = "dataset", all.x = TRUE)
fwrite(summary_dt, file.path(tmp_dir, "dataset_summary.csv"))

message("Aggregated ", length(datasets), " datasets | lv_stats rows: ", nrow(lv_stats),
        " | mechanisms: ", nrow(mech))

write_manifest(tmp_dir, method = "aggregate projections",
               parameters = list(datasets = datasets, groups = unname(groups)),
               extra = list(n_lv_stats = nrow(lv_stats), n_mechanism_rows = nrow(mech),
                            n_pathway_hits = nrow(hits)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
