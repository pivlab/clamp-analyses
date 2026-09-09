source(here::here("scripts/panels/supplementary/style.R"))
base <- here("output/03_model_biology/02_archs4/04_crispercas")
paths <- c(fgsea=file.path(base,"00_gene_enrichment_CRISPRCas9/fgsea_all_lvs.csv"),
           significant_lvs=file.path(base,"00_gene_enrichment_CRISPRCas9/lv_sig_sub_ordered.csv"),
           ora_table=file.path(base,"01_biology_LVs/ora_table.csv"),
           traits_table=file.path(base,"01_biology_LVs/traits_table.csv"),
           lincs_table=file.path(base,"01_biology_LVs/lincs_table.csv"))
out <- nm_sources(7,paths)
d <- fread(paths[1]);sig<-fread(paths[2]);ora<-fread(paths[3]);traits<-fread(paths[4]);lincs<-fread(paths[5])
# Exactly the main notebook's minimum padj and maximum NES across repeats.
best <- d[,.(best_padj=min(padj),best_NES=max(NES)),by=.(lv,pathway)]
stopifnot(setequal(best[best_padj<.05,lv],sig$lv))
fwrite(best,file.path(out,"best_adjusted_p_by_lv_and_gene_set.csv"))
A <- ggplot(best,aes(-log10(pmax(best_padj,.Machine$double.xmin)),fill=pathway))+
 geom_histogram(bins=60,alpha=.6,position="identity",linewidth=.2)+
 geom_vline(xintercept=-log10(.05),linetype="dashed",colour="grey30",linewidth=.3)+
 scale_fill_manual(values=c("#0072B2","#D55E00"),name="CRISPR gene set")+
 labs(x="−log10(best adjusted p-value across 10 repeats)",y="LV count",
      subtitle=sprintf("%d significant LVs of %d; threshold = 0.05",nrow(sig),uniqueN(d$lv)))+nm_theme()+theme(legend.position="bottom")
count_terms <- function(x) vapply(x,function(s) if(is.na(s)||!nzchar(s)) 0L else length(strsplit(s,"; ",fixed=TRUE)[[1]]),1L)
summary <- copy(ora)
for(nm in setdiff(names(summary),"LV"))set(summary,j=nm,value=count_terms(summary[[nm]]))
summary <- merge(summary,traits,by="LV",all=TRUE);summary[,traits:=count_terms(traits)]
summary <- merge(summary,lincs,by="LV",all=TRUE);summary[,lincs_top_drugs:=count_terms(lincs_top_drugs)]
counts <- melt(summary,id.vars="LV",variable.name="collection",value.name="n_terms")
fwrite(counts,file.path(out,"annotation_counts_by_lv.csv"))
collection_names <- c(bp="GO:BP",canonical="Canonical pathways",cellmarker="CellMarker",archs4_cell_lines="ARCHS4 cell lines",azimuth="Azimuth",traits="Traits",lincs_top_drugs="LINCS drugs")
counts[,collection:=factor(collection,levels=names(collection_names),labels=collection_names)]
B <- ggplot(counts,aes(collection,n_terms))+geom_boxplot(fill="#0072B2",alpha=.35,width=.6,outlier.shape=NA,linewidth=.25)+
 geom_point(position=position_jitter(width=.12),size=.6,alpha=.4)+
 scale_y_continuous(trans="log1p",breaks=c(0,1,10,100,1000))+
 labs(x=NULL,y="Annotations per significant LV (log1p scale)")+nm_theme()+theme(axis.text.x=element_text(angle=35,hjust=1))
ids <- sig$lv
heat <- function(selected) {
 dd<-copy(counts[LV%in%selected]);dd[,LV:=factor(LV,levels=rev(selected))]
 ggplot(dd,aes(collection,LV,fill=log1p(n_terms)))+geom_tile(colour="white",linewidth=.2)+
  geom_text(aes(label=n_terms,colour=n_terms>30),size=5/.pt,show.legend=FALSE)+
  scale_colour_manual(values=c("FALSE"="black","TRUE"="white"))+
  scale_fill_gradient(low="white",high="#0072B2",name="Annotation count",breaks=log1p(c(0,10,100,1000)),labels=c(0,10,100,1000))+
  labs(x=NULL,y=NULL)+nm_theme()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=5),
                                     axis.text.x=element_text(angle=40,hjust=1,size=5),legend.position="bottom")
}
pages<-list(plot_grid(nm_tag(A,"A","CRISPR enrichment across all latent variables"),nm_tag(B,"B","Biological annotations of significant LVs"),ncol=2),
            plot_grid(nm_tag(heat(head(ids,45)),"C","Significant LVs, ranks 1–45"),nm_tag(heat(tail(ids,-45)),"D","Significant LVs, ranks 46–89"),ncol=2))
# Retain full network-support tables from main; network and TWAS visuals are in Figure 3I.
for(nm in c("network_lv_summary","network_lv_significant_traits","network_lv_significant_pathways","community_genes","community_edges","trait_supported_gene_edges","inter_lv_shared_gene_edges","poster_style_twas_matrix"))
 file.copy(file.path(base,"02_network/figures",paste0(nm,".csv")),file.path(out,paste0(nm,".csv")),overwrite=TRUE)
nm_export(list(plot_grid(plotlist=pages,ncol=1,rel_heights=c(65,176))),7,"CRISPR-Cas9 enrichment and all 89 significant latent variables")
