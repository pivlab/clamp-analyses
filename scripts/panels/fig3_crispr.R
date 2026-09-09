# Rebuild the two main-branch 02_network notebook plots from its source tables.
crispr_dir <- file.path(DATA_DIR, "04_crispercas/02_network/figures")
crispr_tables <- lapply(c("community_genes", "community_edges", "trait_supported_gene_edges",
                        "inter_lv_shared_gene_edges", "network_lv_summary",
                        "network_lv_significant_traits", "network_lv_significant_pathways",
                        "poster_style_twas_matrix"), function(n) {
  d <- fread(file.path(crispr_dir, paste0(n, ".csv")))
  fwrite(d, file.path(source_dir, paste0("I_", n, ".csv")))
  d
})
names(crispr_tables) <- c("genes", "edges", "trait_edges", "shared", "summary", "traits", "pathways", "twas")
# Twelve module instances, four per column, matching the supplied sketch.
icol <- c("Alzheimer’s"="#8870AA", Diabetes="#477BA8", "High cholesterol"="#378A80",
          "Fat mass"="#C8668C", "Biliary diseases"="#869746", "Coronary diseases"="#C18539")
column_groups <- list(c("Alzheimer’s","Fat mass"),c("Diabetes","Biliary diseases"),
                      c("High cholesterol","Coronary diseases"))
centers_I <- rbindlist(lapply(seq_along(column_groups),function(j) {
  groups <- column_groups[[j]]
  m <- rbindlist(lapply(groups,function(g) copy(crispr_tables$summary[trait_group == g])))
  stopifnot(nrow(m)==4L)
  m[, `:=`(cx=12+(j-1)*43,cy=c(61,43,25,7),module_id=paste(trait_group,LV))]
  m
}))
centers_I[, title := paste0(LV," - ",pathway_cell_type)]
anchors_I <- centers_I[, .(ax=unique(cx)+17,ay=mean(cy)),by=trait_group]
anchors_I[, label := vapply(trait_group,function(s) paste(strwrap(s,width=11),collapse="\n"),"")]
anchors_I[trait_group %in% c("High cholesterol","Coronary diseases"), label := trait_group]
# Stop each connector outside the trait label, including its padding.
anchors_I[, edge_x := ax - vapply(strsplit(label,"\n"), function(lines) max(nchar(lines)), 1L)*.95/2 - 1]
nodes_I <- merge(crispr_tables$genes,centers_I[,.(trait_group,LV,cx,cy,module_id)],by=c("trait_group","LV"),sort=FALSE)
nodes_I[, angle := pi/2+2*pi*(seq_len(.N)-1)/.N,by=module_id]
nodes_I[, `:=`(x=cx+3.2*cos(angle),y=cy+3.2*sin(angle))]
nodes_I[, `:=`(tx=cx+4*cos(angle),ty=cy+4*sin(angle),
               hj=ifelse(cos(angle)>.1,0,ifelse(cos(angle)< -.1,1,.5)))]
circles_I <- centers_I[,.(x=cx+5.5*cos(seq(0,2*pi,length.out=100)),
                          y=cy+5.5*sin(seq(0,2*pi,length.out=100))),by=.(trait_group,module_id)]
path_edges_I <- merge(crispr_tables$edges,nodes_I[,.(trait_group,LV,gene1=gene,x,y)],by=c("trait_group","LV","gene1"))
path_edges_I <- merge(path_edges_I,nodes_I[,.(trait_group,LV,gene2=gene,xend=x,yend=y)],by=c("trait_group","LV","gene2"))
trait_edges_I <- merge(crispr_tables$trait_edges,nodes_I,by=c("trait_group","LV","gene"))
trait_edges_I <- merge(trait_edges_I,anchors_I,by="trait_group")
fwrite(nodes_I,file.path(source_dir,"I_network_node_positions.csv"))
fwrite(centers_I,file.path(source_dir,"I_network_module_positions.csv"))
plot_I1 <- ggplot() +
  geom_polygon(data=circles_I,aes(x,y,group=module_id,fill=trait_group),alpha=.13) +
  geom_segment(data=path_edges_I,aes(x,y,xend=xend,yend=yend),colour="grey65",linewidth=.15,linetype="dashed") +
  geom_segment(data=trait_edges_I,aes(x,y,xend=cx+8,yend=cy,colour=trait_group),linewidth=.18,alpha=.6) +
  geom_segment(data=unique(trait_edges_I[,.(trait_group,cx,cy,ax,ay,edge_x)]),
               aes(x=cx+8,y=cy,xend=edge_x,yend=ay,colour=trait_group),linewidth=.22) +
  geom_point(data=nodes_I,aes(x,y,shape=recovered,colour=trait_group),size=1.1) +
  geom_text(data=nodes_I,aes(tx,ty,label=gene,hjust=hj),size=5/ggplot2::.pt,fontface="italic") +
  geom_text(data=centers_I,aes(x=cx-11,y=cy+6.7,label=title),hjust=0,vjust=.5,
            size=5/ggplot2::.pt,lineheight=1) +
  geom_label(data=anchors_I,aes(ax,ay,label=label,fill=trait_group),
             size=5/ggplot2::.pt,linewidth=0,label.padding=unit(.3,"mm")) +
  scale_shape_manual(values=c("TRUE"=8,"FALSE"=16),guide="none") +
  scale_colour_manual(values=icol,guide="none") +
  scale_fill_manual(values=scales::alpha(icol,.3),guide="none") +
  coord_fixed(xlim=c(0,129),ylim=c(0,70),expand=FALSE,clip="off") +
  theme_void()+theme(plot.margin=margin(0,0,0,0))
twas_I <- copy(crispr_tables$twas)
twas_I[, gene := factor(gene, levels=rev(unique(gene)))]
twas_I[, trait := factor(trait, levels=unique(trait))]
missed <- c("DGAT2", "SOX9", "PCYT2", "PTEN", "HILPDA")
plot_I2 <- ggplot(twas_I, aes(trait,gene)) +
  geom_point(aes(size=neg_log10_p_capped60, fill=neg_log10_p_capped60), shape=21, stroke=.2, colour="#879B95") +
  scale_size_area(max_size=3, guide="none") +
  scale_fill_gradient(low="#DCEFEA",high="#1B4D43", limits=c(0,60), name="S-MultiXcan\n−log10(p)\n(capped at 60)") +
  scale_x_discrete(position="top") +
  labs(x=NULL,y=NULL) + theme_nm() +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text.x=element_text(angle=55,hjust=0,size=5),
        axis.text.y=element_text(face=ifelse(levels(twas_I$gene) %in% missed,"bold.italic","italic"),size=5),
        panel.grid.major=element_line(colour="grey92",linewidth=.2),
        legend.position="bottom",legend.key.width=unit(5,"mm"),legend.key.height=unit(1.5,"mm"))
plot_I <- plot_grid(plot_I1, plot_I2, nrow=1, rel_widths=c(2.6,1))
