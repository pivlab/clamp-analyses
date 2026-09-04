# Figure helpers for the projection report notebooks. S

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

PROJ_MODEL_COLOURS <- c(ARCHS4 = "#332288", local = "#DDCC77")
PROJ_DIVERGING <- list(low = "#1a9850", mid = "white", high = "#d73027")

proj_mechanism_label <- function(x) {
  out <- gsub("_", " ", x)
  out[out == "ISG core"] <- "Interferon alpha response"
  out
}

PROJ_CYTOKINE_MECHANISM_ORDER <- proj_mechanism_label(c(
  "ER_stress_UPR", "immune_crosstalk", "ISG_core", "MHC_class_I",
  "proteasome_ubiquitination", "splicing_RBP"))

proj_local_model_label <- function(dataset) {
  labels <- c(
    cyt_ifna_GSE133218 = "Cyt model",
    mono_lps_GSE193336 = "Monocyte model",
    placenta_EMTAB6701 = "Placenta model"
  )
  unname(labels[[dataset]] %||% "Local model")
}

proj_pretty_term <- function(x) {
  acr <- c("cd\\d+", "ifn[ag]?", "il\\d*", "dna", "rna", "mhc", "er", "upr", "nfkb",
           "jak", "stat", "irf\\d*", "traf\\d*", "atf\\d*", "tnf[a]?", "\\d+s", "hla")
  x <- gsub("^(REACTOME|KEGG_MEDICUS|KEGG|WP|PID|BIOCARTA|GOCC|GOBP|HALLMARK|HALL)_", "", x)
  x <- gsub("\\b([a-z])", "\\U\\1", tolower(gsub("_", " ", x)), perl = TRUE)
  for (a in acr) x <- gsub(paste0("\\b(", a, ")\\b"), "\\U\\1", x, perl = TRUE, ignore.case = TRUE)
  trimws(gsub("\\s+", " ", x))
}

proj_theme <- function(base = 12) {
  theme_bw(base_size = base) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 9, lineheight = 1.08, hjust = 1),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.title = element_text(size = 11), legend.text = element_text(size = 10))
}

proj_expected_mechanism_panel <- function(lv_stats, comparisons, dsel, msel) {
  include_cell_types <- any(comparisons[dataset == dsel & category == "cell_type"]$category == "cell_type")
  s <- comparisons[dataset == dsel & model == msel & recovered == TRUE &
                     (category == "molecular_mechanism" |
                      (include_cell_types & category == "cell_type"))]
  if (!nrow(s)) return(NULL)
  if (include_cell_types) {
    s[, cellmarker_rank := fifelse(database == "cellmarker", 0L, 1L)]
    setorder(s, mechanism, category, cellmarker_rank, ora_fdr)
    s <- s[, .SD[1L], by = .(mechanism, category)]
    s[, cellmarker_rank := NULL]
  }
  s[, paper_mechanism := fifelse(category == "cell_type",
                                 paste0("cell type: ", gsub("_", " ", mechanism)),
                                 proj_mechanism_label(mechanism))]
  stats <- copy(lv_stats[dataset == dsel & model == msel])
  d <- rbindlist(lapply(seq_len(nrow(s)), function(i) {
    target <- s[i]
    keep <- if (is.na(target$target_contrast) || !nzchar(target$target_contrast)) {
      candidates <- stats[LV == target$LV]
      candidates[which.max(abs(logFC))]
    } else if (identical(target$target_contrast, "pooled")) {
      stats[contrast != "pooled" & LV == target$LV]
    } else {
      stats[contrast == target$target_contrast & LV == target$LV]
    }
    if (!nrow(keep)) return(NULL)
    keep[, `:=`(comparison_id = target$comparison_id,
                gene_set = target$gene_set,
                category = target$category,
                paper_mechanism = target$paper_mechanism,
                ora_fdr = target$ora_fdr)]
    keep
  }), fill = TRUE)
  if (!nrow(d)) return(NULL)
  is_cytokine <- identical(dsel, "cyt_ifna_GSE133218")
  contrasts <- unique(d$contrast)
  if (is_cytokine) {
    contrasts <- intersect(c("tp_2h", "tp_8h", "tp_24h"), contrasts)
    d[, contrast := factor(contrast, levels = contrasts)]
    d[, contrast_x := as.integer(contrast)]
  } else {
    contrasts <- "selected"
    d[, contrast_x := 1L]
  }
  d[, row_label := paste(LV, proj_pretty_term(gene_set), sep = "\n")]
  d[, row_key := paste(comparison_id, LV, gene_set, paper_mechanism, sep = "||")]
  row_order <- if (identical(dsel, "cyt_ifna_GSE133218")) {
    unique(d[order(match(paper_mechanism, PROJ_CYTOKINE_MECHANISM_ORDER), ora_fdr)]$row_key)
  } else if (include_cell_types) {
    unique(d[order(factor(category, levels = c("molecular_mechanism", "cell_type")), ora_fdr)]$row_key)
  } else {
    unique(d[order(ora_fdr)]$row_key)
  }
  d[, row_key := factor(row_key, levels = rev(row_order))]
  row_labels <- unique(d[, .(row_key, row_label)])
  row_label_map <- setNames(row_labels$row_label, as.character(row_labels$row_key))
  annotation <- unique(d[, .(row_key, paper_mechanism)])
  annotation[, paper_mechanism_label := vapply(
    paper_mechanism,
    function(x) paste(strwrap(x, width = 27L), collapse = "\n"),
    character(1)
  )]
  lim <- max(abs(d$logFC), na.rm = TRUE)

  p_dot <- ggplot(d, aes(contrast_x, row_key, fill = logFC, size = -log10(pmax(ora_fdr, 1e-300)))) +
    geom_point(shape = 21, colour = "grey40", stroke = 0.3) +
    scale_x_continuous(breaks = if (is_cytokine) seq_along(contrasts) else NULL,
                       labels = if (is_cytokine) contrasts else NULL,
                       limits = c(0.55, length(contrasts) + 0.45),
                       expand = expansion(mult = 0)) +
    scale_fill_gradient2(low = PROJ_DIVERGING$low, mid = PROJ_DIVERGING$mid,
                         high = PROJ_DIVERGING$high, midpoint = 0,
                         limits = c(-lim, lim), name = "LV log2FC") +
    scale_size_continuous(range = c(1.5, 8), name = "-log10 ORA FDR") +
    guides(fill = guide_colourbar(barheight = grid::unit(28, "pt"),
                                  barwidth = grid::unit(6, "pt")),
           size = guide_legend(keyheight = grid::unit(8, "pt"),
                               keywidth = grid::unit(8, "pt"))) +
    scale_y_discrete(labels = row_label_map) +
    labs(x = NULL, y = NULL,
         title = ifelse(msel == "ARCHS4", "ARCHS4 projection",
                        proj_local_model_label(dsel))) +
    proj_theme() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5),
          plot.title.position = "panel",
          axis.text.x = element_text(angle = if (length(contrasts) > 4) 45 else 0,
                                     hjust = if (length(contrasts) > 4) 1 else 0.5),
          legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(size = 7.5, lineheight = 0.9),
          legend.text = element_text(size = 7),
          legend.key.height = grid::unit(8, "pt"),
          legend.key.width = grid::unit(8, "pt"),
          legend.spacing.y = grid::unit(1, "pt"),
          legend.box.spacing = grid::unit(2, "pt"))

  p_annotation <- ggplot(annotation, aes(x = 0, y = row_key, label = paper_mechanism_label)) +
    geom_text(hjust = 0, size = 3.7, lineheight = 0.98, colour = "grey25") +
    scale_y_discrete(limits = levels(d$row_key), drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(title = if (msel == "local") {
      if (include_cell_types) "Expected mechanism / cell type" else "Expected mechanism"
    } else NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 12),
          plot.margin = margin(5.5, 12, 5.5, 0))

  list(dot = p_dot, annotation = p_annotation)
}

proj_expected_mechanism_heatmap <- function(lv_stats, comparisons, dsel) {
  models <- c("local", "ARCHS4")
  panels <- lapply(models, function(m) proj_expected_mechanism_panel(lv_stats, comparisons, dsel, m))
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]
  models <- models[keep]
  if (!length(panels)) return(NULL)
  n_rows <- vapply(models, function(m) uniqueN(comparisons[
    dataset == dsel & model == m & recovered == TRUE &
      (category == "molecular_mechanism" |
       category == "cell_type"),
    .(mechanism, LV)]), integer(1))
  heights <- pmax(1, n_rows + 1L)
  model_column <- patchwork::wrap_plots(lapply(panels, `[[`, "dot"), ncol = 1) +
    patchwork::plot_layout(heights = heights)
  mechanism_column <- patchwork::wrap_plots(lapply(panels, `[[`, "annotation"), ncol = 1) +
    patchwork::plot_layout(heights = heights)
  patchwork::wrap_plots(list(model_column, mechanism_column), nrow = 1, widths = c(1, 0.9)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right", legend.box = "vertical",
          legend.title = element_text(size = 7.5, lineheight = 0.9),
          legend.text = element_text(size = 7),
          legend.key.height = grid::unit(8, "pt"),
          legend.key.width = grid::unit(8, "pt"),
          legend.spacing.y = grid::unit(1, "pt"),
          legend.box.spacing = grid::unit(2, "pt"))
}
