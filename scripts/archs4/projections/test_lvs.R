#!/usr/bin/env Rscript
# Per-LV differential activity between conditions (limma lmFit/eBayes/topTable),
# for one model of one dataset, over every contrast declared for the dataset.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(jsonlite)
  library(limma)
})

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]])))
source(file.path(dirname(script_dir), "common.R"))

args <- parse_cli()
b_path <- required_arg(args, "b")
samples_path <- required_arg(args, "samples")
out_dir <- required_arg(args, "out_dir")
contrasts_spec <- fromJSON(required_arg(args, "contrasts"), simplifyVector = FALSE)
summary_path <- args$summary %||% ""
lv_fdr <- as.numeric(args$lv_fdr %||% 0.05)
lv_logfc <- as.numeric(args$lv_logfc %||% 0.05)
pathway_fdr <- as.numeric(args$pathway_fdr %||% 0.05)
pathway_auc <- as.numeric(args$pathway_auc %||% 0.70)

started <- Sys.time()
B <- read_matrix_csv(b_path)
samples <- fread(samples_path)
if (!"sample" %in% names(samples)) stop("samples.csv needs a `sample` column")

samples <- samples[match(colnames(B), samples$sample)]
if (anyNA(samples$sample)) stop("Some columns of B have no row in samples.csv")
message("B: ", nrow(B), " LVs x ", ncol(B), " samples")

annot <- NULL
if (nzchar(summary_path) && file.exists(summary_path)) {
  s <- fread(summary_path)
  lv_col <- intersect(c("LV", "LV_index"), names(s))[1]
  if (!is.na(lv_col)) {
    setnames(s, lv_col, "LV")
    s[, LV := ifelse(grepl("^LV", LV), LV, paste0("LV", LV))]
    annot <- s[FDR < pathway_fdr & AUC > pathway_auc][order(LV, -AUC, FDR),
      .(prior_pathways = paste(pathway, collapse = " | "),
        top_AUC = max(AUC), n_paths = .N), by = LV]
  }
}

# limma_groups: ~0 + factor, makeContrasts. limma_coef: ~ block + factor, topTable(coef=).
fit_contrast <- function(spec) {
  fac <- as.character(samples[[spec$factor]])
  if (identical(spec$type, "limma_groups")) {
    levels_ <- unlist(spec$levels)
    grp <- factor(fac, levels = levels_)
    if (anyNA(grp)) stop("Samples outside declared levels for contrast ", spec$name)
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- levels(grp)
    fit <- lmFit(B, design)
    cm <- makeContrasts(contrasts = spec$contrast, levels = design)
    fit <- eBayes(contrasts.fit(fit, cm))
    topTable(fit, n = Inf, adjust.method = "BH", sort.by = "none")
  } else if (identical(spec$type, "limma_coef")) {
    grp <- relevel(factor(fac), ref = spec$reference)
    design <- if (!is.null(spec$block)) {
      blk <- factor(as.character(samples[[spec$block]]))
      model.matrix(~ blk + grp)
    } else {
      model.matrix(~ grp)
    }
    colnames(design) <- make.names(sub("^grp", spec$factor, colnames(design)))
    coef_name <- make.names(spec$coef)
    if (!coef_name %in% colnames(design)) {
      stop("Coefficient ", spec$coef, " not in design (have: ",
           paste(colnames(design), collapse = ", "), ")")
    }
    fit <- eBayes(lmFit(B, design))
    topTable(fit, coef = coef_name, n = Inf, adjust.method = "BH", sort.by = "none")
  } else {
    stop("Unknown contrast type: ", spec$type)
  }
}

tmp_dir <- begin_publish(out_dir)
published <- FALSE
on.exit(if (!published && dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE), add = TRUE)

stats <- rbindlist(lapply(contrasts_spec, function(spec) {
  tt <- fit_contrast(spec)
  dt <- as.data.table(tt, keep.rownames = "LV")
  dt[, `:=`(contrast = spec$name, abs_logFC = abs(logFC))]
  message("  ", spec$name, ": ",
          sum(dt$adj.P.Val < lv_fdr & dt$abs_logFC >= lv_logfc), " significant LVs")
  dt
}), fill = TRUE)

if (!is.null(annot)) {
  stats <- merge(stats, annot, by = "LV", all.x = TRUE)
  stats[is.na(prior_pathways), prior_pathways := "unannotated"]
  stats[is.na(top_AUC), top_AUC := 0]
  stats[is.na(n_paths), n_paths := 0L]
} else {
  stats[, `:=`(prior_pathways = "unannotated", top_AUC = 0, n_paths = 0L)]
}
setcolorder(stats, c("contrast", "LV", "logFC", "abs_logFC", "AveExpr", "t",
                     "P.Value", "adj.P.Val"))
setorder(stats, contrast, adj.P.Val, -abs_logFC)
fwrite(stats, file.path(tmp_dir, "lv_stats.csv"))

sig <- stats[adj.P.Val < lv_fdr & abs_logFC >= lv_logfc]
sig_lvs <- unique(sig$LV)
writeLines(sig_lvs, file.path(tmp_dir, "significant_lvs.txt"))

per_lv <- sig[, .(n_sig_contrasts = .N,
                  contrasts = paste(contrast, collapse = ","),
                  max_abs_logFC = max(abs_logFC), min_FDR = min(adj.P.Val),
                  direction = paste(sort(unique(ifelse(logFC > 0, "up", "down"))), collapse = "/"),
                  prior_pathways = first(prior_pathways)),
              by = LV][order(-max_abs_logFC)]
fwrite(per_lv, file.path(tmp_dir, "significant_lvs.csv"))
message("Significant LVs: ", length(sig_lvs), " of ", nrow(B))

write_manifest(tmp_dir, method = "limma::lmFit/eBayes/topTable",
               inputs = c(list(b = b_path, samples = samples_path),
                          if (nzchar(summary_path)) list(summary = summary_path) else list()),
               parameters = list(lv_fdr = lv_fdr, lv_logfc = lv_logfc,
                                 pathway_fdr = pathway_fdr, pathway_auc = pathway_auc,
                                 contrasts = vapply(contrasts_spec, function(s) s$name, "")),
               extra = list(n_lvs = nrow(B), n_samples = ncol(B),
                            n_significant_lvs = length(sig_lvs)),
               started = started)

finish_publish(tmp_dir, out_dir)
published <- TRUE
