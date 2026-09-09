# Figure helpers for the projection report notebooks.
#
# Sourced by notebooks only: every function here takes the aggregated CSVs and
# returns a ggplot. Nothing in this file computes statistics, so a notebook can
# never silently become the place an analysis lives.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

PROJ_MODEL_COLOURS <- c(ARCHS4 = "#332288", local = "#DDCC77")
PROJ_DIVERGING <- list(low = "#1a9850", mid = "white", high = "#d73027")
PROJ_CYTOKINE_MECHANISM_ORDER <- gsub("_", " ", c(
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

proj_wrap <- function(x, width = 36) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), "",
         USE.NAMES = FALSE)
}

# MSigDB terms are ALL_CAPS_UNDERSCORE; make them readable without mangling the
# acronyms that carry the meaning.
proj_pretty_term <- function(x) {
  acr <- c("cd\\d+", "ifn[ag]?", "il\\d*", "dna", "rna", "mhc", "er", "upr", "nfkb",
           "jak", "stat", "irf\\d*", "traf\\d*", "atf\\d*", "tnf[a]?", "\\d+s", "hla")
  x <- gsub("^(REACTOME|KEGG_MEDICUS|KEGG|WP|PID|BIOCARTA|GOCC|GOBP|HALLMARK|HALL)_", "", x)
  x <- gsub("\\b([a-z])", "\\U\\1", tolower(gsub("_", " ", x)), perl = TRUE)
  for (a in acr) x <- gsub(paste0("\\b(", a, ")\\b"), "\\U\\1", x, perl = TRUE, ignore.case = TRUE)
  trimws(gsub("\\s+", " ", x))
}

proj_theme <- function(base = 7) {
  theme_bw(base_size = base) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 3.8, lineheight = 0.82, hjust = 1),
          plot.title = element_text(size = 7, hjust = 0.5),
          plot.subtitle = element_text(size = 5.5, colour = "grey25", lineheight = 1.1),
          legend.title = element_text(size = 5.5), legend.text = element_text(size = 5))
}

# Per-timepoint LV activity: colour = log2FC, size = -log10 FDR. Rows are the
# significant LVs of one model, labelled by their strongest pathway hit.
proj_lv_dotplot <- function(lv_stats, panel_ready, dsel, msel,
                            contrasts = NULL, max_lvs = 12, title = NULL) {
  sig <- panel_ready[dataset == dsel & model == msel][order(-max_abs_logFC)]
  if (!nrow(sig)) return(NULL)
  sig <- head(sig, max_lvs)
  d <- lv_stats[dataset == dsel & model == msel & LV %chin% sig$LV]
  if (!is.null(contrasts)) d <- d[contrast %chin% contrasts]
  lab <- sig[, .(LV, label = paste0(LV, "\n",
                 proj_wrap(ifelse(is.na(best_term), "no significant pathway",
                                  proj_pretty_term(best_term)))))]
  d <- merge(d, lab, by = "LV")
  d[, label := factor(label, levels = rev(lab$label[match(sig$LV, lab$LV)]))]
  d[, contrast := factor(contrast, levels = unique(lv_stats$contrast))]
  d[, negLogFDR := -log10(pmax(adj.P.Val, 1e-300))]
  lim <- max(abs(d$logFC), na.rm = TRUE)
  ggplot(d, aes(contrast, label, fill = logFC, size = negLogFDR)) +
    geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
    scale_fill_gradient2(low = PROJ_DIVERGING$low, mid = PROJ_DIVERGING$mid,
                         high = PROJ_DIVERGING$high, midpoint = 0,
                         limits = c(-lim, lim), name = "log2FC") +
    scale_size_continuous(range = c(1.2, 9), name = "-log10 FDR") +
    labs(x = NULL, y = NULL, title = title %||% paste(dsel, msel)) +
    proj_theme()
}

# The headline comparison.  Tiles use the selected ORA gene set: colour is the
# share of its full membership among the selected LV's top-loading genes.
proj_mechanism_heatmap <- function(mech, comparisons, dsel, title = NULL) {
  d <- mech[dataset == dsel & category %in% c("molecular_mechanism", "cell_type")]
  if (!nrow(d)) return(NULL)
  # Use exactly the same paired canonical gene set as the expected-mechanism
  # and Figure-2E panels.  The older mechanism summary used each model's own
  # best text-matched ORA term, which could inflate a model's recovery count.
  # Molecular programmes use their declared pathway collection (canonical C2
  # by default and Hallmark for IFN-alpha). Expected cell types use direct
  # CellMarker evidence when available.
  cs_raw <- comparisons[dataset == dsel &
                          (category == "molecular_mechanism" | category == "cell_type")]
  if (!nrow(cs_raw)) return(NULL)
  cs_raw[, cellmarker_rank := fifelse(category == "cell_type" & database == "cellmarker", 0L, 1L)]
  setorder(cs_raw, model, mechanism, category, cellmarker_rank, ora_fdr)
  cs_raw <- cs_raw[, .SD[1L], by = .(model, mechanism, category)]
  cs <- cs_raw[,
                    .(comparison_recovered = any(recovered),
                      comparison_n_lvs = uniqueN(LV[recovered]),
                      comparison_n_in_top = max(n_gene_set_in_top_loading, na.rm = TRUE),
                      comparison_n_universe = max(n_gene_set_in_universe, na.rm = TRUE)),
                    by = .(model, mechanism, category)]
  d <- merge(d, cs, by = c("model", "mechanism", "category"), all.x = TRUE, sort = FALSE)
  d[testable == TRUE, `:=`(
    recovered = fifelse(is.na(comparison_recovered), FALSE, comparison_recovered),
    n_lvs = fifelse(is.na(comparison_n_lvs), 0L, comparison_n_lvs),
    n_gene_set_in_top_loading = fifelse(is.na(comparison_n_in_top), 0L, comparison_n_in_top),
    n_gene_set_in_universe = fifelse(is.na(comparison_n_universe), 0L, comparison_n_universe))]
  d[testable == TRUE,
    gene_set_fraction_in_top_loading := fifelse(n_gene_set_in_universe > 0,
                                                n_gene_set_in_top_loading / n_gene_set_in_universe,
                                                0)]
  d[is.na(gene_set_fraction_in_top_loading), gene_set_fraction_in_top_loading := 0]
  d[, lab := fifelse(!testable,
                     "Not testable\n(no canonical set)",
                     fifelse(recovered,
                             sprintf("%d LV%s\n%d/%d set genes", n_lvs, ifelse(n_lvs == 1, "", "s"),
                                     n_gene_set_in_top_loading, n_gene_set_in_universe),
                             "Not recovered"))]
  d[, mechanism_label := fifelse(category == "cell_type",
                                 paste0("cell type: ", gsub("_", " ", mechanism)),
                                 gsub("_", " ", mechanism))]
  mechanism_order <- if (identical(dsel, "cyt_ifna_GSE133218")) {
    PROJ_CYTOKINE_MECHANISM_ORDER
  } else {
    c(sort(unique(d[category == "molecular_mechanism"]$mechanism_label)),
      sort(unique(d[category == "cell_type"]$mechanism_label)))
  }
  d[, mechanism_label := factor(mechanism_label, levels = rev(mechanism_order))]
  ggplot(d, aes(model, mechanism_label, fill = gene_set_fraction_in_top_loading)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), size = 2.5, lineheight = 0.95,
              colour = ifelse(d$gene_set_fraction_in_top_loading > 0.6, "white", "grey15")) +
    scale_fill_gradient(low = "#f7f7f7", high = "#332288", limits = c(0, 1),
                        name = "gene-set\nfraction") +
    labs(x = NULL, y = NULL, title = title %||% dsel,
         subtitle = "tile = expected programme or cell-type gene-set members in selected significant LVs' top-loading genes") +
    proj_theme()
}

# One bar per dataset per model: mechanisms recovered out of those defined.
proj_mechanism_scoreboard <- function(mech) {
  s <- mech[testable == TRUE,
            .(n = .N, recovered = sum(recovered, na.rm = TRUE)), by = .(dataset, group, model)]
  s[, frac := recovered / n]
  ggplot(s, aes(reorder(dataset, frac), frac, fill = model)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%d/%d", recovered, n)),
              position = position_dodge(width = 0.75), hjust = -0.2, size = 2.8) +
    scale_fill_manual(values = PROJ_MODEL_COLOURS) +
    scale_y_continuous(limits = c(0, 1.15), labels = scales::percent) +
    coord_flip() +
    labs(x = NULL, y = "mechanisms recovered", title = "Projection vs local model") +
    proj_theme()
}

# Expected-mechanism heatmap/dotplot. Rows are only LVs whose ORA hit matches a
# pathway or cell marker anticipated from the paper; the annotation at right
# states the corresponding paper mechanism.
proj_expected_mechanism_panel <- function(lv_stats, comparisons, dsel, msel) {
  include_cell_types <- any(comparisons[dataset == dsel & category == "cell_type"]$category == "cell_type")
  s <- comparisons[dataset == dsel & model == msel & recovered == TRUE &
                     (category == "molecular_mechanism" |
                      (include_cell_types & category == "cell_type"))]
  if (!nrow(s)) return(NULL)
  # One row is one paper mechanism / expected gene set, not every available
  # cell-type contrast. Pooled cytokine targets retain their time-course rows;
  # every other target is shown only in its own one-versus-rest contrast.
  # A cell type can have both a canonical-pathway and a CellMarker hit.  Use
  # the marker hit when available; it is the direct cell-type evidence.
  if (include_cell_types) {
    s[, cellmarker_rank := fifelse(database == "cellmarker", 0L, 1L)]
    setorder(s, mechanism, category, cellmarker_rank, ora_fdr)
    s <- s[, .SD[1L], by = .(mechanism, category)]
    s[, cellmarker_rank := NULL]
  }
  s[, paper_mechanism := fifelse(category == "cell_type",
                                 paste0("cell type: ", gsub("_", " ", mechanism)),
                                 gsub("_", " ", mechanism))]
  stats <- copy(lv_stats[dataset == dsel & model == msel])
  d <- rbindlist(lapply(seq_len(nrow(s)), function(i) {
    target <- s[i]
    keep <- if (is.na(target$target_contrast) || !nzchar(target$target_contrast)) {
      # Non-cytokine datasets have no pre-specified time/contrast.  One
      # selected LV must therefore contribute one representative contrast,
      # rather than every cell-type one-versus-rest comparison.
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
    # The other datasets are one selected contrast per recovered LV.  They are
    # deliberately rendered as one model column, not a sparse grid made from
    # the union of unrelated cell-type contrasts.
    contrasts <- "selected"
    d[, contrast_x := 1L]
  }
  # Mechanisms can legitimately share a selected LV and ORA term (notably
  # maternal-immune and dNK mechanisms in placenta).  Keep a unique plotting
  # key so their paper labels cannot be drawn on top of one another.
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
  # The expected (paper-defined) mechanism is a true second table column, not
  # text drawn over the plotting margin.  This keeps it readable for the
  # one-contrast datasets, where an artificial x-coordinate otherwise creates
  # a very wide empty panel and clips the annotation against the legends.
  annotation[, paper_mechanism_label := vapply(
    paper_mechanism,
    function(x) paste(strwrap(x, width = 20L), collapse = "\n"),
    character(1)
  )]
  lv_pathway <- unique(d[, .(row_key, LV, gene_set)])
  lv_pathway[, lv_pathway_label := paste0(
    LV, "\n",
    vapply(gene_set, function(x) {
      term <- proj_pretty_term(x)
      if (nchar(term) > 18L) paste0(substr(term, 1L, 15L), "...") else term
    }, character(1))
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
    scale_size_continuous(range = c(1.8, 5.0), name = "-log10 ORA FDR") +
    # The adjacent expected-mechanism column is the row key for this compact
    # panel. Omitting the duplicate LV/pathway strings leaves space for the
    # dot heatmap itself at figure scale.
    scale_y_discrete(labels = NULL) +
    labs(x = NULL, y = NULL,
         title = ifelse(msel == "ARCHS4", "ARCHS4 projection",
                        proj_local_model_label(dsel))) +
    proj_theme() +
    theme(plot.margin = margin(2, 2, 2, 2),
          plot.title.position = "panel",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = if (length(contrasts) > 4) 45 else 0,
                                     hjust = if (length(contrasts) > 4) 1 else 0.5),
          legend.position = "right",
          legend.key.size = unit(2.5, "mm"),
          legend.box = "vertical")

  p_annotation <- ggplot(annotation, aes(x = 0, y = row_key, label = paper_mechanism_label)) +
    geom_text(hjust = 0, size = 1.7, lineheight = 0.78, colour = "grey25") +
    scale_y_discrete(limits = levels(d$row_key), drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(title = if (msel == "local") "Mechanism" else NULL) +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 7),
          plot.margin = margin(2, 6, 2, 0))

  p_lv_pathway <- ggplot(lv_pathway, aes(x = 0, y = row_key, label = lv_pathway_label)) +
    geom_text(hjust = 0, size = 1.35, lineheight = 0.78, colour = "grey25") +
    scale_y_discrete(limits = levels(d$row_key), drop = FALSE) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    labs(title = if (msel == "local") "Pathway" else NULL) +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 7),
          plot.margin = margin(2, 6, 2, 0))

  list(dot = p_dot, lv_pathway = p_lv_pathway, annotation = p_annotation)
}

# One output per dataset: retain the model-specific colour and size scales from
# the reference panels, then stack the local and projected-model panels.
proj_expected_mechanism_heatmap <- function(lv_stats, comparisons, dsel) {
  models <- c("local", "ARCHS4")
  panels <- lapply(models, function(m) proj_expected_mechanism_panel(lv_stats, comparisons, dsel, m))
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]
  models <- models[keep]
  if (!length(panels)) return(NULL)
  # A model with two recovered rows should not receive the same vertical area
  # as a model with five. The +1 preserves room for titles, axes and legends.
  n_rows <- vapply(models, function(m) uniqueN(comparisons[
    dataset == dsel & model == m & recovered == TRUE &
      (category == "molecular_mechanism" |
       category == "cell_type"),
    .(mechanism, LV)]), integer(1))
  heights <- pmax(1, n_rows + 1L)
  model_column <- patchwork::wrap_plots(lapply(panels, `[[`, "dot"), ncol = 1) +
    patchwork::plot_layout(heights = heights)
  lv_column <- patchwork::wrap_plots(lapply(panels, `[[`, "lv_pathway"), ncol = 1) +
    patchwork::plot_layout(heights = heights)
  mechanism_column <- patchwork::wrap_plots(lapply(panels, `[[`, "annotation"), ncol = 1) +
    patchwork::plot_layout(heights = heights)
  patchwork::wrap_plots(list(model_column, lv_column, mechanism_column), nrow = 1,
                        widths = c(1.45, 0.65, 0.65)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right", legend.box = "horizontal",
          legend.key.height = unit(2, "mm"), legend.key.width = unit(2, "mm"))
}

# Same content as proj_expected_mechanism_heatmap(), composed with cowplot
# instead of patchwork's guides="collect". patchwork reserves a fixed-size
# (absolute-unit) guide gutter that does not shrink with its container, so at
# very small allocated widths it can crowd out the dot/annotation columns
# entirely. cowplot::plot_grid uses purely relative units throughout, which
# degrades gracefully when nested inside another plot_grid() at a tight size.
proj_expected_mechanism_heatmap_cp <- function(lv_stats, comparisons, dsel,
                                                legend_position = c("right", "bottom"),
                                                legend_size = 0.28) {
  legend_position <- match.arg(legend_position)
  models <- c("local", "ARCHS4")
  panels <- lapply(models, function(m) proj_expected_mechanism_panel(lv_stats, comparisons, dsel, m))
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]
  models <- models[keep]
  if (!length(panels)) return(NULL)
  n_rows <- vapply(models, function(m) uniqueN(comparisons[
    dataset == dsel & model == m & recovered == TRUE &
      (category == "molecular_mechanism" | category == "cell_type"),
    .(mechanism, LV)]), integer(1))
  heights <- pmax(1, n_rows + 1L)

  shared_legend <- cowplot::get_legend(
    panels[[1]]$dot + theme(
      legend.position = legend_position,
      legend.box = if (legend_position == "bottom") "horizontal" else "vertical",
      legend.key.size = unit(2.5, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5.5)
    )
  )
  dots_nolegend <- lapply(panels, function(p) p$dot + theme(legend.position = "none"))
  model_column <- cowplot::plot_grid(plotlist = dots_nolegend, ncol = 1, rel_heights = heights)
  lv_column <- cowplot::plot_grid(plotlist = lapply(panels, `[[`, "lv_pathway"),
                                   ncol = 1, rel_heights = heights)
  mech_column  <- cowplot::plot_grid(plotlist = lapply(panels, `[[`, "annotation"),
                                     ncol = 1, rel_heights = heights)
  body <- cowplot::plot_grid(model_column, lv_column, mech_column, nrow = 1,
                             rel_widths = c(1.45, 0.65, 0.65))
  if (legend_position == "bottom") {
    cowplot::plot_grid(body, shared_legend, ncol = 1,
                       rel_heights = c(1 - legend_size, legend_size))
  } else {
    cowplot::plot_grid(body, shared_legend, nrow = 1,
                       rel_widths = c(1 - legend_size, legend_size))
  }
}

# Figure 2E-style view: every gene is placed at its rank in the selected
# significant LV; all members of the selected ORA gene set are coloured.
proj_gene_set_pair_plot <- function(comparisons, loadings, dsel, comparison_key,
                                    facet_nrow = 1L, compact_x_label = FALSE,
                                    title = NULL,
                                    strip_model_only = FALSE,
                                    summary_position = c("top", "bottom")) {
  summary_position <- match.arg(summary_position)
  s <- comparisons[dataset == dsel & get("comparison_id") == comparison_key]
  if (!nrow(s)) return(NULL)
  s[, model := factor(model, levels = c("ARCHS4", "local"))]
  d <- loadings[dataset == dsel & get("comparison_id") == comparison_key]
  # Figure 2E view: show the ranked leading 1% only, rather than the full
  # loading tail. The recovery labels retain the full expected-set denominator.
  d <- d[in_top_loading_set == TRUE]
  d[, model := factor(model, levels = c("ARCHS4", "local"))]
  if (!nrow(d[model == "ARCHS4"])) return(NULL)

  labels <- d[is_gene_set == TRUE][order(model, rank), head(.SD, 5L), by = model]
  s[, label := fifelse(recovered,
                        sprintf("%d/%d in top %.0f%%", n_gene_set_in_top_loading,
                                n_gene_set_in_universe, 100 * top_pct),
                        "Not recovered")]
  s[nzchar(comparison_note), label := paste(label, comparison_note, sep = "\n")]
  missing <- s[recovered == FALSE]
  strip_labels <- if (strip_model_only) {
    setNames(ifelse(as.character(s$model) == "ARCHS4", "ARCHS4 projection",
                    proj_local_model_label(dsel)), as.character(s$model))
  } else {
    setNames(paste0(as.character(s$model), "\n", fifelse(s$recovered,
                                                           proj_wrap(proj_pretty_term(s$gene_set), width = 20),
                                                           "Not recovered")),
             as.character(s$model))
  }

  ggplot(d, aes(rank, loading)) +
    geom_point(data = d[is_gene_set == FALSE], colour = "grey78", size = 0.12,
               alpha = 0.65) +
    geom_point(data = d[is_gene_set == TRUE], aes(colour = model), size = 0.5,
               alpha = 0.95) +
    ggrepel::geom_text_repel(data = labels, aes(label = gene, colour = model),
                            size = 1.25, fontface = "italic", seed = 42,
                            max.overlaps = Inf, min.segment.length = 0,
                            segment.size = 0.15, segment.colour = "grey55",
                            box.padding = 0.12, point.padding = 0.06,
                            show.legend = FALSE) +
    geom_text(data = s[recovered == TRUE],
              aes(x = Inf, y = if (summary_position == "top") Inf else -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05,
              vjust = if (summary_position == "top") 1.15 else -0.25, size = 1.5,
              colour = "grey20") +
    geom_text(data = missing,
              aes(x = Inf, y = if (summary_position == "top") Inf else -Inf, label = label),
              inherit.aes = FALSE, hjust = 1.05,
              vjust = if (summary_position == "top") 1.15 else -0.25, size = 1.7,
              colour = "grey30") +
    facet_wrap(~model, nrow = facet_nrow, scales = "free_y", drop = FALSE,
               labeller = labeller(model = as_labeller(strip_labels))) +
    scale_colour_manual(values = PROJ_MODEL_COLOURS) +
    scale_x_continuous(expand = expansion(mult = c(0.015, 0.08))) +
    labs(x = if (compact_x_label) "Gene rank" else "Top 1% genes ranked by descending loading",
         y = "Gene loading",
         title = title %||% gsub("_", " ", s$mechanism[1])) +
    theme_classic(base_size = 6) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 5, lineheight = 0.85),
          plot.title = element_text(size = 7, hjust = 0.5),
          panel.grid = element_blank(),
          legend.position = "none")
}
