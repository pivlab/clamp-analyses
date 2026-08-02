# Shared GTEx helpers used by multiple production scripts.
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    key <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) stop("Missing --", gsub("_", "-", name))
  value
}

# Build the GO:BP pathway prior used by clamp.R and plier.R.
# genes: character vector of gene symbols to match the prior against.
# gmt_path: path to a pinned local copy of the GO_Biological_Process gene-set
#   library (same source pseudobulk's rule pathway_prior downloads), read
#   with CLAMP:::read_gmt() instead of a live Enrichr call.
# Returns: list(matched = <matched pathway matrix>, chat = <getChat() result>).
build_go_bp_prior <- function(genes, gmt_path) {
  bp_gmt <- CLAMP:::read_gmt(gmt_path)
  names(bp_gmt) <- paste0("BP_", names(bp_gmt))
  gmt_list <- list(BP = bp_gmt)

  path_mat <- CLAMP::gmtListToSparseMat(gmt_list)
  matched <- CLAMP::getMatchedPathwayMat(path_mat, genes)
  chat <- CLAMP::getChat(matched)

  list(matched = matched, chat = chat)
}
