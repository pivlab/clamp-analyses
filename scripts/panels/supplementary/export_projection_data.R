# Run from the main repository to reproduce its curated selections and tests.
suppressPackageStartupMessages({library(data.table);library(here);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
script_root <- normalizePath(args[1])
source(here("scripts/archs4/projections/plots.R"))
agg_dir <- here("output/03_model_biology/02_archs4/03_projections/aggregate")
prod_root <- here("output/01_model_building/02_archs4/03_projections")
out <- here("output/99_panels/supp5/source_data")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(group in c("00_cytokines","01_monocyte","02_placenta")) {
 source(file.path(script_root,paste0("supp5_",group,"_data.R")))
 for(n in c("manual_recovery","gene_loadings","comparisons"))
  fwrite(get(n),file.path(out,paste0(group,"_",n,".csv")))
 rm(z_mats);gc()
}
file.copy(file.path(agg_dir,"lv_stats_long.csv"),file.path(out,"lv_stats_long.csv"),overwrite=TRUE)

for(group in c("03_pseudobulk_recovery","04_gtex_ari")) source(file.path(script_root,paste0("supp5_",group,"_data.R")))
