#!/usr/bin/env Rscript
# GenomicSuperSignature pseudobulk models
suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(GenomicSuperSignature)
  library(S4Vectors)
  library(SummarizedExperiment)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
source(file.path(script_dir, "common.R"))
args <- parse_cli()
dataset <- required_arg(args, "dataset")
loaded <- read_norm_and_k(required_arg(args, "norm"), required_arg(args, "k"))
d_cluster <- as.integer(if (is.null(args$d_cluster)) 4L else args$d_cluster)
set.seed(as.integer(if (is.null(args$seed)) 123L else args$seed))

pca <- prcomp(t(loaded$norm))
rotation <- pca$rotation[, seq_len(loaded$k), drop = FALSE]
colnames(rotation) <- paste0(dataset, ".PC", seq_len(loaded$k))
eigs <- pca$sdev^2
variance <- rbind(SD = sqrt(eigs), Variance = eigs / sum(eigs), Cumulative = cumsum(eigs) / sum(eigs))
variance <- variance[, seq_len(loaded$k), drop = FALSE]
colnames(variance) <- colnames(rotation)

k_clust <- max(round(ncol(rotation) / d_cluster), 2L)
distance <- as.dist(1 - cor(rotation, method = "spearman"))
clusters <- cutree(hclust(distance, method = "ward.D"), k = k_clust)
avg <- GenomicSuperSignature::buildAvgLoading(rotation, k_clust, cluster = clusters)
sil <- tryCatch(cluster::silhouette(avg$cluster, distance), error = function(e) NULL)
avg$sw <- if (is.null(sil)) rep(NA_real_, k_clust) else summary(sil)$clus.avg.widths

training <- S4Vectors::DataFrame(PCAsummary = I(list(variance))); rownames(training) <- dataset
model <- GenomicSuperSignature::PCAGenomicSignatures(
  assays = list(RAVindex = as.matrix(avg$avgLoading)), trainingData = training
)
S4Vectors::metadata(model) <- avg[c("cluster", "size", "k", "n")]
names(S4Vectors::metadata(model)$size) <- paste0("RAV", seq_len(ncol(model)))
GenomicSuperSignature::geneSets(model) <- "Custom"
GenomicSuperSignature::studies(model) <- avg$studies
GenomicSuperSignature::silhouetteWidth(model) <- avg$sw

index <- SummarizedExperiment::assays(model)[["RAVindex"]] |> as.matrix()
common <- sort(intersect(rownames(loaded$norm), rownames(index)))
B <- crossprod(index[common, , drop = FALSE], loaded$norm[common, , drop = FALSE])
rownames(B) <- paste0("RAV", seq_len(nrow(B))); colnames(B) <- colnames(loaded$norm)
Z <- index[common, , drop = FALSE]; colnames(Z) <- paste0("RAV", seq_len(ncol(Z)))
out <- required_arg(args, "out_dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(B, file.path(out, "B.csv")); write.csv(Z, file.path(out, "Z.csv"))
saveRDS(model, file.path(out, "RAVmodel.rds"))
message("GSS complete")
