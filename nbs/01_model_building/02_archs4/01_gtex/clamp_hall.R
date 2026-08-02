#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(CLAMP)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
gtex_genes <- readRDS(required_arg(args, "genes"))
samples <- readRDS(required_arg(args, "samples"))
gtex_fbm_filt <- readRDS(required_arg(args, "fbm_filt"))
gtex_svdRes <- readRDS(required_arg(args, "svd_res"))
gtex_baseRes <- readRDS(required_arg(args, "base_rds"))
CLAMP_K_gtex <- readRDS(required_arg(args, "k"))
pathways_dir <- required_arg(args, "pathways_dir")
out_dir <- required_arg(args, "out_dir")
multiplier <- as.integer(args$multiplier %||% 100L)
max_iter <- as.integer(args$max_iter %||% 5000L)

message("Genes: ", length(gtex_genes))
message("Samples: ", length(samples))
message("CLAMP K: ", CLAMP_K_gtex)

# Load and match all pathways prior
hall_gmt <- CLAMP:::read_gmt(file.path(pathways_dir, "h.all.v2026.1.Hs.symbols.gmt"))
reactome_gmt <- CLAMP:::read_gmt(file.path(pathways_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"))
gocc_gmt <- CLAMP:::read_gmt(file.path(pathways_dir, "c5.go.cc.v2026.1.Hs.symbols.gmt"))
c8_gmt <- CLAMP:::read_gmt(file.path(pathways_dir, "c8.all.v2026.1.Hs.symbols.gmt"))

names(hall_gmt) <- paste0("HALL_", names(hall_gmt))
names(reactome_gmt) <- paste0("REACTOME_", names(reactome_gmt))
names(gocc_gmt) <- paste0("GOCC_", names(gocc_gmt))
names(c8_gmt) <- paste0("C8_", names(c8_gmt))

all_pathways_list <- list(
  HALL = hall_gmt,
  REACTOME = reactome_gmt,
  GOCC = gocc_gmt,
  C8 = c8_gmt
)

all_pathways_pathMat <- gmtListToSparseMat(all_pathways_list)
all_pathways_matched <- getMatchedPathwayMat(all_pathways_pathMat, gtex_genes)
message("Loaded and matched all pathways matrix against GTEx genes")

# Run CLAMPfull with all pathways prior
message("Running CLAMPfull with all pathways prior on GTEx (all samples)...")

gtex_fullRes_hall <- CLAMPfull(
  Y = gtex_fbm_filt,
  svdres = gtex_svdRes,
  priorMat = all_pathways_matched,
  clamp.base.result = gtex_baseRes,
  use_cpp = TRUE,
  trace = TRUE,
  multiplier = multiplier,
  max.iter = max_iter,
  clamp_k = CLAMP_K_gtex
)

gtex_fullRes_hall$Z <- data.frame(gtex_fullRes_hall$Z)
rownames(gtex_fullRes_hall$Z) <- gtex_genes

gtex_fullRes_hall$B <- data.frame(gtex_fullRes_hall$B)
colnames(gtex_fullRes_hall$B) <- samples

gtex_fullRes_hall$summary <- gtex_fullRes_hall$summary %>%
  dplyr::rename(LV = LV_index) %>%
  dplyr::mutate(LV = paste0("LV", LV))

message("CLAMPfull completed")

# Save results
dst_dir <- file.path(out_dir, "CLAMPfull_hall")
dir.create(dst_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(gtex_fullRes_hall, file = file.path(out_dir, "CLAMPfull_hall.rds"))

write.csv(gtex_fullRes_hall$B, file.path(dst_dir, "B.csv"))
write.csv(gtex_fullRes_hall$Z, file.path(dst_dir, "Z.csv"))
write.csv(gtex_fullRes_hall$summary, file.path(dst_dir, "summary.csv"))

message("Results saved to: ", dst_dir)
