source(here::here("scripts/panels/supplementary/style.R"))
suppressPackageStartupMessages(library(ggrepel))

# These inputs are regenerated from the current main-repository notebooks by
# export_projection_data.R; they retain the curated mechanisms and contrasts.
out <- here("output/99_panels/supp5/source_data")
groups <- c("00_cytokines", "01_monocyte", "02_placenta")
S <- lapply(groups, function(g) fread(file.path(out, paste0(g, "_comparisons.csv"))))
L <- lapply(groups, function(g) fread(file.path(out, paste0(g, "_gene_loadings.csv"))))
stats <- fread(file.path(out, "lv_stats_long.csv"))
source(here("../clamp-analyses/scripts/archs4/projections/plots.R"))

# A: the four original notebook data sets, restored as conventional vertical plots.
method_colours <- c("CLAMPfull"="#3e348b", "ARCHS4 projection"="#923155",
                    "RNA-Seq"="#9c9f36", "ARCHS4 projection null"="#999999")
make_benchmark <- function(file, title, metric) {
  d <- fread(file.path(out, file))
  if ("ARI" %in% names(d)) setnames(d, "ARI", "value")
  m <- d[, .(mean_value=mean(value), max_value=max(value)), by=method][order(-mean_value)]
  # Match the paired-comparison treatment used in Fig. 2A.  Each benchmark
  # row is aligned across methods before testing, then shown directly above
  # the corresponding boxplots.
  methods <- as.character(m$method)
  pairs <- if ("RNA-Seq" %in% methods) {
    list(c("CLAMPfull", "ARCHS4 projection"), c("CLAMPfull", "RNA-Seq"),
         c("ARCHS4 projection", "ARCHS4 projection null"))
  } else {
    list(c("CLAMPfull", "ARCHS4 projection"),
         c("ARCHS4 projection", "ARCHS4 projection null"))
  }
  comparisons <- rbindlist(lapply(pairs, function(pair) {
    a <- d[method == pair[1], value]
    b <- d[method == pair[2], value]
    data.table(left=pair[1], right=pair[2], p=wilcox.test(a, b, paired=TRUE)$p.value)
  }))
  comparisons[, label:=paste0("q = ", formatC(p.adjust(p, "BH"), format="e", digits=1))]
  d[, method:=factor(method, levels=m$method)]
  m[, method:=factor(method, levels=levels(d$method))]
  ymax <- max(d$value, na.rm=TRUE)
  yrange <- diff(range(d$value, na.rm=TRUE)); if (!is.finite(yrange) || yrange == 0) yrange <- 1
  comparisons[, `:=`(x=match(left, levels(d$method)), xend=match(right, levels(d$method)),
                     y=ymax + yrange*(.24 + .28*(.I-1)), label_y=ymax + yrange*(.30 + .28*(.I-1)))]
  ggplot(d, aes(method, value, fill=method)) +
    geom_boxplot(width=.56, outlier.shape=NA, linewidth=.3, alpha=.8) +
    geom_point(position=position_jitter(width=.08, height=0), shape=21, size=.8,
               fill="white", stroke=.2) +
    geom_point(data=m, aes(x=method, y=mean_value), shape=23, size=1.6, fill="white",
               colour="black", inherit.aes=FALSE) +
    geom_text(data=m, aes(x=method, y=max_value, label=sprintf("%.3f", mean_value)),
              vjust=-.65, size=5/.pt, inherit.aes=FALSE) +
    geom_segment(data=comparisons, aes(x=x, xend=xend, y=y, yend=y), inherit.aes=FALSE, linewidth=.2) +
    geom_segment(data=comparisons, aes(x=x, xend=x, y=y-yrange*.025, yend=y), inherit.aes=FALSE, linewidth=.2) +
    geom_segment(data=comparisons, aes(x=xend, xend=xend, y=y-yrange*.025, yend=y), inherit.aes=FALSE, linewidth=.2) +
    geom_text(data=comparisons, aes(x=(x+xend)/2, y=label_y, label=label),
              inherit.aes=FALSE, size=5/.pt) +
    scale_fill_manual(values=method_colours, guide="none") +
    labs(title=title, x=NULL, y=metric) + nm_theme() +
    theme(axis.text.x=element_text(angle=42, hjust=1, size=5),
          plot.title=element_text(size=5.5, hjust=.5),
          plot.margin=margin(2,2,2,2,"mm")) +
    coord_cartesian(ylim=c(min(d$value, na.rm=TRUE), ymax + yrange*(.24 + .28*nrow(comparisons) + .12)), clip="off")
}
A <- plot_grid(
  make_benchmark("04_gtex_ari_k_matched.csv", "GTEx, K = 578", "ARI"),
  make_benchmark("04_gtex_ari_full_lvs.csv", "GTEx, 1,728 LVs", "ARI"),
  make_benchmark("03_pseudobulk_recovery_k_matched.csv", "Pseudobulk, K-matched", "Max Pearson r"),
  make_benchmark("03_pseudobulk_recovery_full_lvs.csv", "Pseudobulk, 1,728 LVs", "Max Pearson r"),
  ncol=4)

# B-D: the Figure 3 projection-panel treatment. Local and ARCHS4 use separate
# log2FC scales because the model ranges differ; point size is −log10(ORA FDR).
effect_rows <- lapply(seq_along(S), function(i) {
  s <- copy(S[[i]]); s[, ora_fdr:=fdr]
  rbindlist(lapply(c("local", "ARCHS4"), function(m) {
    p <- proj_expected_mechanism_panel(stats, s, unique(s$dataset), m)
    if (is.null(p)) return(NULL)
    d <- copy(p$dot$data)
    d[, mechanism:=comparison_id]
    d
  }), fill=TRUE)
})
dot_group <- function(i) {
  s <- copy(S[[i]])
  h <- effect_rows[[i]]
  order <- unique(s$comparison_id)
  levels_y <- rev(order)
  h[, mechanism:=factor(mechanism, levels=levels_y)]
  size_limits <- range(-log10(pmax(h$ora_fdr, 1e-300)))
  # This is intentionally the Figure 3 G mechanism-recovery construction:
  # ARCHS4 at left, local at right, wide ARCHS4 field, independent fill keys
  # and one shared ORA-FDR size key.
  dot_one <- function(model_name, show_y) {
    d <- h[model==model_name]
    miss <- s[model==model_name & !recovered]
    miss[, mechanism:=factor(comparison_id, levels=levels_y)]
    x_missing <- if (i==1L) "tp_8h" else "selected"
    lim <- max(abs(scales::breaks_pretty(n=3)(range(c(-abs(d$logFC), abs(d$logFC))))))
    ggplot(d, aes(contrast, mechanism)) +
      geom_point(aes(size=-log10(pmax(ora_fdr,1e-300)), fill=logFC), shape=21, stroke=.2) +
      geom_text(data=miss, aes(x=x_missing, y=mechanism), label="Not recovered",
                inherit.aes=FALSE, size=5/.pt, colour="grey35", fontface="italic") +
      scale_size_continuous(range=c(.7,2.4), limits=size_limits, name="−log10 ORA FDR",
                            breaks=scales::breaks_pretty(n=3)) +
      scale_fill_gradient2(low="#1a9850", mid="white", high="#d73027", midpoint=0,
                           limits=c(-lim,lim), name="LV log2FC") +
      scale_x_discrete(labels=function(x) sub("tp_", "", x), drop=FALSE) +
      scale_y_discrete(drop=FALSE) + labs(title=if(model_name=="local") "Local model" else "ARCHS4 projection", x=NULL, y=NULL) +
      nm_theme() + theme(axis.text.y=if(show_y) element_text(size=5) else element_blank(),
                         axis.ticks.y=if(show_y) element_line() else element_blank(),
                         axis.text.x=element_text(size=5, angle=45, hjust=1), plot.title=element_text(size=5.5,hjust=.5),
                         legend.position="bottom", legend.box="vertical", legend.key.width=unit(3,"mm"),legend.key.height=unit(1.2,"mm"),
                         plot.margin=margin(1,1,1,1,"mm"))
  }
  arch <- dot_one("ARCHS4", TRUE); local <- dot_one("local", FALSE)
  size_legend <- get_legend(arch + guides(fill="none") + theme(legend.position="right"))
  fill_legend_arch <- get_legend(arch + guides(size="none", fill=guide_colourbar(title.position="right", barwidth=unit(1.2,"mm"),barheight=unit(10,"mm"))) + theme(legend.position="right"))
  fill_legend_local <- get_legend(local + guides(size="none", fill=guide_colourbar(title.position="right", barwidth=unit(1.2,"mm"),barheight=unit(10,"mm"))) + theme(legend.position="right"))
  local <- local + theme(legend.position="none"); arch <- arch + theme(legend.position="none")
  dots <- plot_grid(arch, local, nrow=1, rel_widths=c(1.9,1), align="h", axis="tb")
  legend_column <- plot_grid(fill_legend_arch, fill_legend_local, size_legend, ncol=1, rel_heights=c(.34,.34,.32))
  plot_grid(dots, legend_column, nrow=1, rel_widths=c(1,.28))
}
# Cytokine recovery is already shown in Fig. 3 and is therefore omitted here.
notebook_summary_table <- function(i) {
  d <- copy(S[[i]])
  d[, fraction:=fifelse(recovered & n_gene_set_in_universe > 0,
                        n_gene_set_in_top_loading / n_gene_set_in_universe, NA_real_)]
  d[, label:=fifelse(recovered, sprintf("%d/%d", n_gene_set_in_top_loading, n_gene_set_in_universe), "—")]
  d[, model_label:=factor(fifelse(model == "ARCHS4", "ARCHS4", "Local"), levels=c("ARCHS4", "Local"))]
  d[, mechanism:=factor(comparison_id, levels=rev(unique(comparison_id)))]
  ggplot(d, aes(model_label, mechanism, fill=fraction)) +
    geom_tile(colour="white", linewidth=.2) +
    geom_text(aes(label=label), size=5/.pt) +
    scale_fill_gradient(low="white", high="#1B9E77", limits=c(0,1), na.value="grey92",
                        name="Fraction", breaks=c(0,.5,1)) +
    labs(x=NULL, y=NULL, title="Notebook summary") + nm_theme() +
    theme(axis.text.x=element_text(size=5, angle=45, hjust=1), axis.text.y=element_text(size=5),
          plot.title=element_text(size=5.5, hjust=.5), legend.position="bottom",
          legend.key.width=unit(8,"mm"), legend.key.height=unit(1.0,"mm"),
          plot.margin=margin(1,1,1,1,"mm"))
}
B <- nm_tag(plot_grid(notebook_summary_table(2), dot_group(2), nrow=1, rel_widths=c(1.15,2.25)), "B", "Monocyte mechanism recovery")
C <- nm_tag(plot_grid(notebook_summary_table(3), dot_group(3), nrow=1, rel_widths=c(1.15,2.25)), "C", "Placenta mechanism recovery")

# Select two examples per model group. Prefer ARCHS4-only recoveries, then the
# largest ARCHS4 advantage in recovered pathway members.
example_candidates <- rbindlist(lapply(seq_along(S), function(g) {
  wide <- dcast(S[[g]], comparison_id + category ~ model,
    value.var=c("recovered", "LV", "gene_set", "fdr", "n_gene_set_in_top_loading", "n_gene_set_in_universe"))
  wide[, advantage:=n_gene_set_in_top_loading_ARCHS4 - n_gene_set_in_top_loading_local]
  only_arch <- wide[recovered_ARCHS4 == TRUE & recovered_local == FALSE][order(-n_gene_set_in_top_loading_ARCHS4)]
  # Show one clear ARCHS4-only recovery and one fair head-to-head example:
  # both models recover it, but ARCHS4 has more top-loading pathway genes and
  # a significant ARCHS4 FDR.
  both <- wide[recovered_ARCHS4 == TRUE & recovered_local == TRUE &
                 advantage > 0 & fdr_ARCHS4 <= .05][order(-advantage, fdr_ARCHS4)]
  chosen <- rbind(only_arch[seq_len(min(1L, nrow(only_arch)))],
                  both[seq_len(min(1L, nrow(both)))])
  if (nrow(chosen) < 2L) {
    fallback <- wide[recovered_ARCHS4 == TRUE & recovered_local == TRUE][order(-advantage, fdr_ARCHS4)]
    chosen <- unique(rbind(chosen, fallback))[seq_len(min(2L, nrow(unique(rbind(chosen, fallback)))))]
  }
  chosen[, group:=g]
  chosen
}), fill=TRUE)
# Two examples are sufficient at final size: one ARCHS4-only recovery, and one
# paired recovery in which ARCHS4 retains more pathway genes with significant FDR.
examples <- example_candidates
stopifnot(examples[, .N, by=group][, all(N == 2L)])

short_fdr <- function(x) ifelse(is.na(x), "—", formatC(x, format="e", digits=1))
example_model_plot <- function(x, model_name) {
  g <- x$group; mechanism_name <- x$comparison_id
  s <- copy(S[[g]][comparison_id == mechanism_name])
  d <- copy(L[[g]][comparison_id == mechanism_name & model == model_name])
  st <- s[model == model_name]
  model_label <- if (model_name == "ARCHS4") "ARCHS4 projection" else "Local model"
  status <- if (st$recovered[1]) sprintf("%s; %d/%d genes; FDR %s", st$LV[1], st$n_gene_set_in_top_loading[1], st$n_gene_set_in_universe[1], short_fdr(st$fdr[1])) else "Not recovered"
  gene_labels <- d[is_gene_set == TRUE][order(rank)][seq_len(min(3L, .N))]
  ggplot(d, aes(rank, loading)) +
    geom_point(colour="grey78", size=.18, alpha=.35) +
    geom_point(data=d[is_gene_set==TRUE], colour=if (model_name == "ARCHS4") "#332288" else "#DDCC77", size=.55) +
    geom_text_repel(data=gene_labels, aes(label=gene), size=5/.pt, min.segment.length=0,
                    seed=1, max.overlaps=Inf, box.padding=.08, point.padding=.03) +
    scale_y_continuous(expand=expansion(mult=c(.03,.30))) +
    labs(title=model_label, subtitle=status, x=if(model_name == "local") "Gene rank" else NULL, y="Loading") + nm_theme() +
    theme(plot.title=element_text(size=5.5, hjust=.5), plot.subtitle=element_text(size=5, hjust=.5),
          axis.text=element_text(size=5), axis.title.y=element_text(size=5), plot.margin=margin(1,2,1,2,"mm"))
}
example_pair <- function(x) {
  g <- x$group; mechanism_name <- x$comparison_id
  s <- S[[g]][comparison_id == mechanism_name & model == "ARCHS4"]
  term <- proj_pretty_term(s$gene_set[1])
  # Each example is a direct side-by-side comparison at the same scale:
  # ARCHS4 | local.
  plot_grid(ggdraw() + draw_label(term, size=6),
            plot_grid(example_model_plot(x, "ARCHS4"), example_model_plot(x, "local"), nrow=1),
            ncol=1, rel_heights=c(.12,.88))
}
E <- plot_grid(plotlist=lapply(seq_len(nrow(examples)), function(i) example_pair(examples[i])), ncol=3)

page <- plot_grid(
  nm_tag(A, "A", "Projection benchmarks: original vertical notebook plots"),
  B, C,
  nm_tag(E, "D", "Two local-versus-ARCHS4 examples each for cytokines, monocytes and placenta; labels give LV, recovered genes and FDR"),
  ncol=1, rel_heights=c(46,40,49,104))
nm_export(list(page), 5, "Projection benchmarks, mechanism recovery and local-versus-ARCHS4 gene-loading examples")
