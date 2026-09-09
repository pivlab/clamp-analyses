source(here::here("scripts/panels/supplementary/style.R"))
set.seed(123)
source(here("scripts/panels/supplementary/supp2_panels.R"))
nm_sources(2,unlist(panel_input))
map_theme <- theme_void(base_family="Helvetica",base_size=5.5)+
 theme(panel.border=element_rect(colour="grey65",fill=NA,linewidth=.2),
       plot.title=element_text(size=5.5,hjust=.5,lineheight=1),plot.margin=margin(1,1,1,1,"mm"))
feature_maps <- function(dataset_id) {
 d<-cells[dataset==dataset_id];meta<-lvs[dataset==dataset_id][order(plot_order)];lim<-square_limits(d)
 ctr<-d[,.(umap1=median(umap1),umap2=median(umap2)),by=mapped_cell_type]
 ctr[,label:=as.character(seq_len(.N))]
 annotation<-ggplot(d,aes(umap1,umap2,colour=mapped_cell_type))+
  ggrastr::rasterise(geom_point(size=.16,alpha=.78),dpi=300)+
  ggrepel::geom_text_repel(data=ctr,aes(label=label),colour="black",size=5/.pt,seed=123,
    box.padding=.4,point.padding=.2,force=5,max.time=3,max.iter=100000,max.overlaps=Inf,min.segment.length=0,segment.size=.15)+
  scale_colour_manual(values=CELL_COLORS)+coord_equal(xlim=lim$x,ylim=lim$y,expand=FALSE)+
  labs(title="Cell-type annotation")+map_theme+theme(legend.position="none")
 features<-lapply(seq_len(nrow(meta)),function(i) {
  m<-meta[i];score_col<-paste0("activity_",i)
  ggplot(d,aes(umap1,umap2,colour=.data[[score_col]]))+
   ggrastr::rasterise(geom_point(size=.16),dpi=300)+
   scale_colour_gradientn(colours=c("#D8E2EF","white",CELL_COLORS[[m$cell_type]]))+
   coord_equal(xlim=lim$x,ylim=lim$y,expand=FALSE)+
   labs(title=paste(paste(strwrap(pretty_cell_type(m$cell_type),width=18),collapse="\n"),m$LV,sep=" · "))+map_theme+theme(legend.position="none")
 })
 spare <- if(dataset_id %in% c("Heart_Datar2026","Lung_Sikkema2023","PBMC_1k1k"))2L else 1L
 nrkey <- ceiling(nrow(ctr)/spare)
 ctr[,`:=`(keycol=(seq_len(.N)-1)%/%nrkey,keyrow=(seq_len(.N)-1)%%nrkey+1)]
 legend <- ggplot(ctr,aes(x=keycol,y=-keyrow))+
  geom_point(aes(colour=mapped_cell_type),size=.7)+
  geom_text(aes(x=keycol+.04,label=paste(label,pretty_cell_type(mapped_cell_type))),hjust=0,size=5/.pt,colour="black")+
  scale_colour_manual(values=CELL_COLORS,guide="none")+
  coord_cartesian(xlim=c(-.04,spare-.08),ylim=c(-nrkey-.3,-.7),expand=FALSE,clip="off")+
  theme_void()+theme(plot.margin=margin(1,1,1,1,"mm"))
 result <- c(list(annotation),features)
 attr(result,"key") <- legend
 result
}
feature_block<-function(ds,ncol,label) {
 maps<-feature_maps(ds);n<-length(maps);nr<-ceiling(n/ncol);last<-(nr-1)*ncol+1
 rows<-lapply(seq_len(nr-1),function(r)plot_grid(plotlist=maps[((r-1)*ncol+1):(r*ncol)],ncol=ncol))
 rows<-c(rows,list(plot_grid(plot_grid(plotlist=maps[last:n],ncol=n-last+1),attr(maps,"key"),ncol=2,
  rel_widths=c(n-last+1,ncol-(n-last+1)))))
 nm_tag(plot_grid(plotlist=rows,ncol=1),label,DATASET_LABELS[[ds]])
}
panel_A <- panel_A+labs(y="Top-1% annotated-cell\npurity (%)")+theme(axis.title.y=element_text(size=5.5))
# Uniform six-column map grid, cohort-specific blocks, and no omitted LVs.
page <- plot_grid(nm_tag(panel_A,"A","Donor-bulk projection purity"),
 feature_block("Heart_Datar2026",6,"B"),
 plot_grid(feature_block("Brain_Mathys2023",4,"C"),feature_block("Brain_Xiong2023",4,"D"),ncol=2),
 feature_block("Lung_Sikkema2023",6,"E"),
 feature_block("PBMC_1k1k",6,"F"),
 feature_block("PBMC_Perez2022",6,"G"),
 ncol=1,rel_heights=c(30,58,44,43,43,23))
nm_export(list(page),2,"Donor-bulk recovery and all 47 cell-type latent-variable maps")
