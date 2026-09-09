source(here::here("scripts/panels/supplementary/style.R"))
set.seed(123)
base <- here("output/03_model_biology/02_archs4")
paths <- c(pathway_coverage=file.path(base,"00_coverage/coverage_long.csv"),
           trait_coverage=file.path(base,"00_coverage/coverage_trait_recovery.csv"),
           trait_full_compendia=file.path(base,"00_coverage/final_trait_recovery.csv"),
           pathway_saturation=file.path(base,"01_saturation/saturation_long.csv"),
           trait_saturation=file.path(base,"01_saturation/saturation_trait_recovery_k1728.csv"))
out <- nm_sources(4,paths)
cov <- fread(paths[1])[database!="reactome"]
traits <- fread(paths[2]); finals <- fread(paths[3]); sat <- fread(paths[4])[database!="reactome"]; tsat <- fread(paths[5])
cols <- unlist(yaml::read_yaml(here("config.yaml"))$MODEL_COLORS)[c("CLAMPbase","CLAMPfull")]
models <- c("CLAMPbase","CLAMPfull")
paths_cov <- cov[dataset=="archs4" & fraction<100,.(recovered=sum(recovered_pathways),eligible=sum(eligible_pathways)),by=.(model,fraction,seed)]
traits[,model:=sub("CLAMPfull_bp","CLAMPfull",model)]
trait_cov <- traits[fraction<100,.(model,fraction,seed,recovered=recovered_traits,eligible=eligible_traits)]
psat <- sat[dataset=="archs4",.(recovered=sum(recovered_pathways),eligible=sum(eligible_pathways)),by=.(model,fraction,k,seed)]
tsat[,model:=sub("CLAMPfull_bp","CLAMPfull",model)]
tsat <- tsat[,.(model,fraction,k,seed,recovered=recovered_traits,eligible=eligible_traits)]
recovery_plot <- function(d,model_name,metric,saturation=FALSE,paired=TRUE,all_tests=FALSE) {
  d <- copy(d[model==model_name]);d[,fraction_label:=factor(paste0(fraction,"%"),levels=paste0(sort(unique(fraction)),"%"))]
  totals <- d[,.(recovered=median(recovered),eligible=median(eligible),minimum=min(recovered)),by=fraction_label]
  span <- max(1,diff(range(d$recovered)))
  totals[,label:=sprintf("%s\n(%.0f%%)",format(round(recovered),big.mark=",",trim=TRUE),100*recovered/eligible)]
  refs <- if(all_tests) setdiff(sort(unique(d$fraction)),75) else c(25,50)
  tests <- rbindlist(lapply(refs,function(f) {
    ref <- d[fraction==75][order(seed)]; cand <- d[fraction==f][order(seed)]
    if(paired) {common<-intersect(ref$seed,cand$seed);ref<-ref[seed%in%common];cand<-cand[seed%in%common]}
    if(nrow(ref)<2 || nrow(cand)<2)return(NULL)
    data.table(fraction=f,p=t.test(ref$recovered,cand$recovered,paired=paired,var.equal=!paired,alternative="greater")$p.value,
               delta=mean(ref$recovered)-mean(cand$recovered))
  }))
  tests[,q:=p.adjust(p,"BH")]
  tests[,`:=`(x=match(paste0(fraction,"%"),levels(d$fraction_label)),xend=match("75%",levels(d$fraction_label)),
               y=max(d$recovered)+span*(.16+.20*seq_len(.N)),label=sprintf("q=%s; %+.0f",nm_q(q),delta))]
  fwrite(tests,file.path(out,paste0(model_name,"_",metric,if(saturation)"_saturation" else "_coverage","_tests.csv")))
  p <- ggplot(d,aes(fraction_label,recovered))+
    geom_boxplot(fill=cols[[model_name]],alpha=.6,width=.6,outlier.shape=NA,linewidth=.25)+
    geom_point(position=position_jitter(width=.08,height=0),size=.7)+
    geom_text(data=totals,aes(y=minimum-span*.10,label=label),vjust=1,size=5.5/.pt,lineheight=.95)+
    geom_segment(data=tests,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,linewidth=.25)+
    geom_text(data=tests,aes(x=(x+xend)/2,y=y+span*.025,label=label),inherit.aes=FALSE,vjust=0,size=5.5/.pt)+
    scale_y_continuous(expand=expansion(mult=c(.20,.12)))+
    labs(x="Studies used",y=paste("Recovered",metric),title=model_name)+nm_theme()
  if(!saturation)p<-p+stat_summary(aes(group=1),fun=median,geom="line",linetype="dashed",linewidth=.3)
  p
}
cross <- cov[fraction==100 & ((model=="CLAMPfull"&seed==1) | (model=="CLAMPbase"&((dataset=="archs4"&seed==1)|(dataset!="archs4"&seed==0)))),.(recovered=sum(recovered_pathways)),by=.(dataset,model)]
finals[,model:=sub("CLAMPfull_bp","CLAMPfull",model)]
cross_plot <- function(d,m,metric) {
  d<-copy(d[model==m]);d[,dataset:=factor(dataset,levels=c("archs4","gtex","recount2"),labels=c("ARCHS4","GTEx","recount2"))]
  ggplot(d,aes(dataset,recovered,colour=dataset))+geom_point(size=2)+
    geom_text(aes(label=format(recovered,big.mark=",",trim=TRUE)),vjust=-1,size=5.5/.pt,colour="black")+
    scale_colour_manual(values=c(ARCHS4="#0072B2",GTEx="#E69F00",recount2="#009E73"),guide="none")+
    scale_y_continuous(expand=expansion(mult=c(.1,.25)))+labs(x=NULL,y=paste("Recovered",metric),title=m)+nm_theme()
}
saturation_plot <- function(m) {
 d<-copy(psat[model==m]);d[,`:=`(K=factor(k),fraction_label=factor(paste0(fraction,"%"),levels=paste0(sort(unique(fraction)),"%")))]
 ggplot(d,aes(K,recovered,fill=fraction_label,colour=fraction_label))+
   geom_boxplot(position=position_dodge(.85),alpha=.35,outlier.shape=NA,linewidth=.25)+
   geom_point(position=position_jitterdodge(jitter.width=.1,dodge.width=.85),size=.5,show.legend=FALSE)+
   stat_summary(aes(group=fraction_label),fun=mean,geom="line",position=position_dodge(.85),linetype="dashed",linewidth=.3)+
   scale_fill_viridis_d(end=.9,direction=-1)+scale_colour_viridis_d(end=.9,direction=-1)+
   labs(x="Latent variables (K)",y="Recovered pathways",title=m,fill="Studies used",colour="Studies used")+nm_theme()+
   theme(legend.position="bottom")+guides(fill=guide_legend(nrow=1),colour=guide_legend(nrow=1))
}
# Complete remaining analyses; main-Figure-3 counterparts are recorded in the inventory.
pages <- list(
 plot_grid(nm_tag(recovery_plot(paths_cov,"CLAMPfull","pathways"),"A","Pathway coverage"),
           nm_tag(recovery_plot(trait_cov,"CLAMPfull","traits"),"B","Trait coverage"),
           nm_tag(cross_plot(cross,"CLAMPfull","pathways"),"C","Full-compendium pathways"),
           nm_tag(cross_plot(finals[,.(dataset,model,recovered=recovered_traits)],"CLAMPbase","traits"),"D","Full-compendium traits"),ncol=2),
 plot_grid(nm_tag(saturation_plot("CLAMPfull"),"E","Pathway saturation"),
           nm_tag(recovery_plot(psat[k==max(k)],"CLAMPfull","pathways",TRUE,TRUE,FALSE),"F","Full-rank pathway recovery"),
           nm_tag(recovery_plot(psat[k==max(k)],"CLAMPbase","pathways",TRUE,TRUE,FALSE),"G","CLAMPbase full-rank recovery"),
           nm_tag(recovery_plot(tsat,"CLAMPfull","traits",TRUE),"H","Full-rank trait recovery"),
           nm_tag(recovery_plot(tsat,"CLAMPbase","traits",TRUE),"I","Full-rank trait recovery"),ncol=2))
nm_export(list(plot_grid(plotlist=pages,ncol=1,rel_heights=c(2,3))),4,"Complete pathway and trait coverage and saturation")
