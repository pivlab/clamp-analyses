#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(matrixStats)
  library(factoextra)
  library(cluster)
  library(GenomicSuperSignature)
  library(msigdbr)
  library(fgsea)
  library(dplyr)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))

args <- parse_cli()
gtex_data <- readRDS(required_arg(args, "df_gtex_fbm_filt"))
n <- readRDS(required_arg(args, "k"))
out_dir <- required_arg(args, "out_dir")
study <- args$dataset %||% "GTEx"
d <- as.integer(args$d_cluster %||% 4L)
seed <- as.integer(args$seed %||% 123L)
set.seed(seed)

# PCA
pca_res <- prcomp(t(as.matrix(gtex_data))) # x is a matrix with genes(row) x samples(column)

trainingData_PCA <- list()
trainingData_PCA[[study]] <- list()

trainingData_PCA[[study]]$rotation <- pca_res$rotation[, 1:n]
colnames(trainingData_PCA[[study]]$rotation) <- paste0(study, ".PC", 1:n)

eigs <- pca_res$sdev^2

pca_summary <- rbind(SD = sqrt(eigs),
                      Variance = eigs / sum(eigs),
                      Cumulative = cumsum(eigs) / sum(eigs))

trainingData_PCA[[study]]$variance <- pca_summary[, 1:n]
colnames(trainingData_PCA[[study]]$variance) <- paste0(study, ".PC", c(1:n))

# Hierarchical clustering
allZ <- trainingData_PCA[[study]]$rotation
storage.mode(allZ) <- "double"
all <- t(allZ)

res.dist <- factoextra::get_dist(all, method = "spearman")

# Cut the tree
k <- round(nrow(all) / d, 0)
res.hcut <- factoextra::hcut(res.dist, k = k, hc_func = "hclust",
                              hc_method = "ward.D", hc_metric = "spearman")

# Build avgLoading
trainingData_PCclusters <- buildAvgLoading(allZ, k, cluster = res.hcut$cluster)

# Silhouette Width
cl <- trainingData_PCclusters$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters$sw <- cl_silh_width # add silhouette width to the result

# Final model
trainingData_df <- DataFrame(
  PCAsummary = I(list(trainingData_PCA[[study]]$variance))
)
rownames(trainingData_df) <- study

# Construct PCAGenomicSignatures
RAVmodel <- PCAGenomicSignatures(
  assays = list(RAVindex = as.matrix(trainingData_PCclusters$avgLoading)),
  trainingData = trainingData_df
)

# Attach metadata analogous to the multi-study build
metadata(RAVmodel) <- trainingData_PCclusters[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))

geneSets(RAVmodel) <- "Custom" # label as you wish
studies(RAVmodel) <- trainingData_PCclusters$studies # PC->study map
silhouetteWidth(RAVmodel) <- trainingData_PCclusters$sw
updateNote(RAVmodel) <- paste0("Single-matrix GTEx model; PCs = ", n, ".")
metadata(RAVmodel)$version <- "0.1.0-single"

B <- assays(RAVmodel)[["RAVindex"]]

# genes x samples (numeric matrix)
expr <- as.matrix(gtex_data)
storage.mode(expr) <- "double"
rownames(expr) <- make.unique(rownames(expr))

# loadings from model: genes x RAVs
RAVindex <- assays(RAVmodel)[["RAVindex"]] |> as.matrix()
storage.mode(RAVindex) <- "double"

# align by genes (same order in both)
common <- sort(intersect(rownames(expr), rownames(RAVindex)))
expr_c <- expr[common, , drop = FALSE]
RAVindex_c <- RAVindex[common, , drop = FALSE]

# RAV x sample scores (B in LV-space)
B_RAVxSample <- crossprod(RAVindex_c, expr_c) # == t(RAVindex_c) %*% expr_c

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(B_RAVxSample, file = file.path(out_dir, "gtex_B.csv"), quote = FALSE)

message("GenomicSuperSignature saved -> ", out_dir)
