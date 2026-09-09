h2_stats <- copy(comparisons_H[dataset == dataset_H & comparison_id == "ISG_core::hallmark"])
h2_points <- copy(loadings_H[dataset == dataset_H & comparison_id == "ISG_core::hallmark" & in_top_loading_set == TRUE])
fwrite(h2_stats,file.path(source_dir,"H2_statistics.csv"))
fwrite(h2_points,file.path(source_dir,"H2_loading_points.csv"))
h2_panels <- lapply(c("ARCHS4","local"), function(msel) {
  d <- h2_points[model == msel]
  st <- h2_stats[model == msel]
  lbl <- d[is_gene_set == TRUE][order(rank)][seq_len(5)]
  yr <- max(d$loading)
  xr <- max(d$rank)
  # Reserve two explicit rows for the gene labels above the data.
  lbl[, `:=`(label_x = xr*c(.05,.40,.75,.22,.62),
             label_y = yr*c(1.25,1.25,1.25,1.05,1.05))]
  p <- ggplot(d,aes(rank,loading)) +
    geom_point(data=d[is_gene_set == FALSE], colour="grey78",size=.15) +
    geom_point(data=d[is_gene_set == TRUE],colour=PROJ_MODEL_COLOURS[[msel]],size=.5) +
    geom_segment(data=lbl,aes(x=rank,y=loading,xend=label_x,yend=label_y),
                 colour="grey60",linewidth=.15) +
    geom_text(data=lbl,aes(x=label_x,y=label_y,label=gene),size=5/ggplot2::.pt,
              fontface="italic",colour="black",hjust=0) +
    annotate("text",x=xr,y=yr*1.52,hjust=1,size=5/ggplot2::.pt,
             label=sprintf("%d/%d in top %.0f%%", st$n_gene_set_in_top_loading,
                           st$n_gene_set_in_universe,100*st$top_pct)) +
    scale_x_continuous(expand=expansion(mult=c(.02,.16)),breaks=c(0,50,100,150)) +
    scale_y_continuous(limits=c(0,yr*1.85),breaks=scales::breaks_pretty(n=3)) +
    labs(x=if(msel == "local") "Gene rank" else NULL,y="Gene loading",
         title=if(msel == "ARCHS4") "ARCHS4 projection" else "Cyt model") + theme_nm() +
    theme(plot.title=element_text(size=5.5,face="bold",hjust=.5),plot.margin=margin(0,4,0,0))
  if(msel == "ARCHS4") p <- p + annotate("text",x=xr,y=yr*1.78,hjust=1,
                                             label=paste0("q = ",fmt_p_G(st$p_adj)),size=5/ggplot2::.pt)
  p
})
plot_H2 <- plot_grid(ggdraw()+draw_label("Interferon Alpha Response",size=6),
                     plot_grid(plotlist=h2_panels,ncol=1),ncol=1,rel_heights=c(.07,.93))
