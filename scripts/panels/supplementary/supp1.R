source(here::here("scripts/panels/supplementary/style.R"))
set.seed(123)
source(here("scripts/panels/supplementary/supp1_panels.R"))
nm_sources(1,unlist(panel_input)[!grepl("notebook",names(panel_input))])
# One portrait page; all six cohorts are retained in both heatmap blocks.
panel_B <- panel_B+theme(axis.text.y=element_text(size=5.5),axis.text.x=element_text(size=5.5))
top <- plot_grid(panel_A,plot_grid(panel_B,panel_C_timing,ncol=1),ncol=2,rel_widths=c(2,1))
purity_plots <- lapply(c(DATASETS_ROW1,DATASETS_ROW2),function(ds)
 make_purity_panel(ds)+theme(legend.position="none",plot.margin=margin(1,1,1,1,"mm")))
purity <- plot_grid(plot_grid(plotlist=purity_plots[1:3],ncol=3),plot_grid(plotlist=purity_plots[4:6],ncol=3),ncol=1,rel_heights=c(1.35,1))
purity_legend <- get_legend(make_purity_panel(DATASETS_ROW1[1])+theme(legend.position="bottom"))
purity <- plot_grid(purity,purity_legend,ncol=1,rel_heights=c(1,.1))
markers <- plot_grid(plot_grid(plotlist=lapply(marker_panels[1:3],function(p)p+theme(legend.position="none")),ncol=3),
 plot_grid(plotlist=lapply(marker_panels[4:6],function(p)p+theme(legend.position="none")),ncol=3),
 marker_legend,ncol=1,rel_heights=c(43,36,7))
page <- plot_grid(top,panel_C_ortho,nm_tag(purity,"E","Cell-type recovery in all six datasets"),
 nm_tag(markers,"F","Marker enrichment for all matched cell types"),ncol=1,rel_heights=c(48,25,82,86))
nm_export(list(page),1,"Pseudobulk benchmarks, orthogonality, recovery and marker enrichment")
