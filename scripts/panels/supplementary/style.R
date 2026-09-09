suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot); library(patchwork)
  library(here); library(grid); library(yaml)
})
NM_WIDTH <- 180
NM_HEIGHT <- 247
NM_FONT <- 5.5
nm_theme <- function() {
  theme_classic(base_size=NM_FONT,base_family="Helvetica") +
    theme(axis.title=element_text(size=6),axis.text=element_text(size=5.5,colour="black"),
          plot.title=element_text(size=6.5,face="bold",hjust=0),
          plot.subtitle=element_text(size=5.5),strip.text=element_text(size=6,face="bold"),
          legend.title=element_text(size=5.5),legend.text=element_text(size=5.5),
          legend.key.size=unit(2.5,"mm"),axis.ticks.length=unit(1,"mm"),
          plot.margin=margin(2,2,2,2,"mm"))
}
nm_tag <- function(p,label,title=NULL) {
  heading <- paste(label,if(!is.null(title)) title else "")
  plot_grid(ggdraw()+draw_label(heading,x=.01,hjust=0,size=6.5,fontface="bold",fontfamily="Helvetica"),
            p,ncol=1,rel_heights=c(.08,.92))
}
nm_fonts <- function(g) {
  if(inherits(g,"text")) {
    fs <- g$gp$fontsize; if(is.null(fs)) fs <- NM_FONT
    g$gp$fontsize <- pmin(7,pmax(5,fs));g$gp$fontfamily <- "Helvetica"
  }
  if(!is.null(g$grobs)) g$grobs <- lapply(g$grobs,nm_fonts)
  if(!is.null(g$children)) for(i in seq_along(g$children)) g$children[[i]] <- nm_fonts(g$children[[i]])
  g
}
nm_export <- function(pages,number,titles) {
  out <- here(sprintf("output/99_panels/supp%d",number))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  stopifnot(length(pages)==1L, length(titles)==1L)
  grobs <- vector("list",length(pages))
  for(i in seq_along(pages)) {
    heading <- sprintf("Supplementary Figure %d%s | %s",number,
                       if(length(pages)>1) sprintf(" (%d/%d)",i,length(pages)) else "",titles[i])
    page <- plot_grid(ggdraw()+draw_label(heading,x=.012,hjust=0,size=7,fontface="bold",fontfamily="Helvetica"),
                       pages[[i]],ncol=1,rel_heights=c(.018,.982))
    # Resolve text on the actual output-sized device before enforcing the font floor.
    grDevices::pdf(NULL,width=NM_WIDTH/25.4,height=NM_HEIGHT/25.4)
    g <- nm_fonts(grid::forceGrob(ggplotGrob(page)))
    grDevices::dev.off()
    grobs[[i]] <- g
    # Every supplementary figure is a single Nature Methods page. Export only
    # the submission artefacts; do not leave intermediate page files behind.
    stem <- sprintf("supp%d",number)
    ggsave(file.path(out,paste0(stem,".svg")),g,width=NM_WIDTH,height=NM_HEIGHT,units="mm",device=svglite::svglite,bg="white")
    ggsave(file.path(out,paste0(stem,".png")),g,width=NM_WIDTH,height=NM_HEIGHT,units="mm",dpi=300,device=ragg::agg_png,bg="white")
  }
  raw <- tempfile(fileext=".pdf")
  cairo_pdf(raw,width=NM_WIDTH/25.4,height=NM_HEIGHT/25.4,onefile=TRUE,family="Helvetica")
  for(g in grobs) {grid.newpage();grid.draw(g)}
  dev.off()
  dest <- file.path(out,sprintf("supp%d.pdf",number))
  if(nzchar(Sys.which("gs"))) {
    rc <- system2("gs",c("-q","-dNOPAUSE","-dBATCH","-sDEVICE=pdfwrite",
                          sprintf("-dDEVICEWIDTHPOINTS=%.6f",NM_WIDTH/25.4*72),
                          sprintf("-dDEVICEHEIGHTPOINTS=%.6f",NM_HEIGHT/25.4*72),
                          "-dFIXEDMEDIA","-dPDFFitPage",paste0("-sOutputFile=",dest),raw))
    stopifnot(rc==0L)
  } else file.copy(raw,dest,overwrite=TRUE)
  unlink(raw)
  cat(sprintf("Supplementary Figure %d: %d pages, %d x %d mm\n",number,length(pages),NM_WIDTH,NM_HEIGHT))
  invisible(grobs)
}
nm_sources <- function(number,paths) {
  out <- here(sprintf("output/99_panels/supp%d/source_data",number));dir.create(out,recursive=TRUE,showWarnings=FALSE)
  stopifnot(all(file.exists(paths)))
  dest <- file.path(out,paste0(names(paths),ifelse(grepl("\\.gz$",paths),".csv.gz",".csv")))
  for(i in seq_along(paths)) file.copy(paths[i],dest[i],overwrite=TRUE)
  fwrite(data.table(source=names(paths),path=paths,md5=unname(tools::md5sum(paths))),file.path(out,"input_manifest.csv"))
  out
}
nm_q <- function(p) ifelse(p<.001,format(p,scientific=TRUE,digits=2),sprintf("%.3f",p))
