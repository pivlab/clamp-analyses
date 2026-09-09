source(here::here("scripts/panels/supplementary/style.R"))
NM_HEIGHT <- 170
set.seed(123)
out <- here("output/99_panels/supp6/source_data")
# Python export is kept separate because the upstream score artifact is a pandas pickle.
stopifnot(file.exists(file.path(out,"prediction_scores.csv")))
scores <- fread(file.path(out,"prediction_scores.csv"));roc <- fread(file.path(out,"roc_coordinates.csv"));perf <- fread(file.path(out,"performance_summary.csv"))
base <- here("output/03_model_biology/02_archs4/02_drug_diseases_canonical/three_model_comparison/data")
paths <- setNames(file.path(base,paste0(c("per_tissue_metrics","max_aggregate_reference","paired_tissue_tests","ordering_stability"),".csv")),
                  c("per_tissue_metrics","max_aggregate_reference","paired_tissue_tests","ordering_stability"))
nm_sources(6,paths)
tissue <- fread(paths[1]);tests <- fread(paths[3]);agg <- fread(paths[2])
order <- c("ARCHS4","recount2","GTEx","Gene-based")
colours <- c(ARCHS4="#1B4D43",recount2="#2E9B8E",GTEx="#E8836B",`Gene-based`="#888888")
scores[,model:=factor(model,levels=order)];roc[,model:=factor(model,levels=order)]
A <- ggplot(scores,aes(score,colour=model))+
  geom_histogram(aes(y=after_stat(density)),bins=45,position="identity",fill=NA,linewidth=.35)+
  scale_colour_manual(values=colours,name="Method")+labs(x="Aggregated drug–disease score",y="Density")+nm_theme()+theme(legend.position="bottom")
roc_labels <- setNames(sprintf("%s (AUC=%.3f)",perf$model,perf$AUROC),perf$model)
B <- ggplot(roc,aes(fpr,tpr,colour=model))+geom_path(linewidth=.55)+
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60",linewidth=.3)+
  scale_colour_manual(values=colours,labels=roc_labels,name="Method")+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  labs(x="False positive rate",y="True positive rate")+nm_theme()+
  theme(legend.position="bottom")+guides(colour=guide_legend(ncol=1))
perf[,model:=factor(model,levels=rev(order))]
C <- ggplot(perf,aes(AUROC,model,colour=model))+geom_vline(xintercept=.5,linetype="dashed",colour="grey60",linewidth=.3)+
  geom_point(size=2)+geom_text(aes(label=sprintf("%.3f",AUROC)),hjust=-.4,size=5.5/.pt,colour="black")+
  scale_colour_manual(values=colours,guide="none")+scale_x_continuous(limits=c(.48,.69),breaks=seq(.5,.65,.05))+
  labs(x="Max-aggregated AUROC",y=NULL,subtitle="Dashed line: random chance (AUROC = 0.5)")+nm_theme()
tests[,comparison:=paste(group1,"vs",group2)]
tests[,label:=paste0("q = ",nm_q(q_value_bh))]
D <- ggplot(tests,aes(-log10(q_value_bh),factor(comparison,levels=rev(comparison))))+
  geom_vline(xintercept=-log10(.05),linetype="dashed",colour="grey60",linewidth=.3)+
  geom_point(size=1.5,colour="#1B4D43")+geom_text(aes(label=label),hjust=-.15,size=5/.pt)+
  scale_x_continuous(expand=expansion(mult=c(.03,.5)))+
  labs(x="−log10(BH-adjusted paired-test p-value)",y=NULL,subtitle="49 tissues; Wilcoxon signed-rank tests")+nm_theme()
# AUROC comparison is in main Figure 3H; preserve every remaining metric in source data.
pages <- list(plot_grid(nm_tag(A,"A","Prediction score distributions"),nm_tag(B,"B","ROC curves"),
                        nm_tag(C,"C","Aggregate performance"),nm_tag(D,"D","All paired method comparisons"),ncol=2))
titles <- "Drug–disease prediction and complete method comparisons"
if("auprc" %in% names(tissue)) {
  tissue[,method:=factor(method,levels=order)]
  E <- ggplot(tissue,aes(method,auprc,fill=method))+geom_boxplot(width=.55,outlier.shape=NA,alpha=.5,linewidth=.3)+
    geom_point(position=position_jitter(width=.08),size=.7,alpha=.6)+
    geom_point(data=agg,aes(method,auprc),shape=23,fill="white",size=2)+
    scale_fill_manual(values=colours,guide="none")+labs(x=NULL,y="AUPRC per GTEx tissue",subtitle="Diamonds: max-aggregate reference")+nm_theme()
  F <- ggplot(tissue,aes(auroc,auprc,colour=method))+geom_point(size=1,alpha=.75)+
    scale_colour_manual(values=colours,name="Method")+labs(x="AUROC per tissue",y="AUPRC per tissue")+nm_theme()+theme(legend.position="bottom")
  # Full 49-tissue values remain readable as two horizontal panels.
  tissue[,tissue:=factor(tissue,levels=rev(sort(unique(tissue))))]
  L <- ggplot(tissue,aes(auroc,tissue,colour=method))+geom_point(position=position_dodge(width=.55),size=.9)+
    scale_colour_manual(values=colours,name="Method")+labs(x="AUROC",y=NULL)+nm_theme()+theme(axis.text.y=element_text(size=5),legend.position="bottom")
  M <- ggplot(tissue,aes(auprc,tissue,colour=method))+geom_point(position=position_dodge(width=.55),size=.9)+
    scale_colour_manual(values=colours,name="Method")+labs(x="AUPRC",y=NULL)+nm_theme()+theme(axis.text.y=element_text(size=5),legend.position="bottom")
  pages <- c(pages,list(plot_grid(nm_tag(E,"E","Precision–recall performance"),nm_tag(F,"F","Per-tissue metric agreement"),ncol=2),
                       plot_grid(nm_tag(L,"G","All 49 tissues: AUROC"),nm_tag(M,"H","All 49 tissues: AUPRC"),ncol=2)))
  titles <- c(titles,"Complementary precision–recall metrics","Tissue-specific performance")
}
nm_export(pages,6,titles)
