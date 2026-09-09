mono_dataset <- "mono_lps_GSE193336"
hits <- fread(file.path(agg_dir, 'lv_pathway_hits_long.csv'))[dataset == mono_dataset]
setnames(hits, 'p.adjust', 'fdr', skip_absent = TRUE)

excluded_hub_lvs <- c("LV1152", "LV1206")
hits <- hits[!LV %chin% excluded_hub_lvs]

manual_mechanisms <- list(
  "NAD+ depletion / SIRT1 inhibition"                       = list(pattern = "NAD_BIOSYNTHESIS|NAD_METABOLISM|SIRTUIN",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Increased glucose transport / glycolytic reprogramming"  = list(pattern = "GLYCOLYSIS|GLUCOSE_METABOLISM|GLUCOSE_TRANSPORT",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Reduced TCA cycle / OXPHOS"                               = list(pattern = "OXIDATIVE_PHOSPHORYLATION|TCA_CYCLE|CITRIC_ACID_CYCLE|RESPIRATORY_ELECTRON_TRANSPORT|\\bOXPHOS\\b",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Oxidative stress / ROS"                                   = list(pattern = "OXIDATIVE_STRESS|REACTIVE_OXYGEN|ROS_AND_RNS|KEAP1|NFE2L2|NRF2",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "HIF1A / NF-κB inflammatory signaling"                     = list(pattern = "HIF1|NF_KB|NFKB|TOLL_LIKE|TLR4|\\bTNF\\b|INFLAMMASOME",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Itaconate / macrophage metabolic rewiring"                = list(pattern = "ITACONATE|IRG1|ACOD1",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Glutathione / redox metabolism"                           = list(pattern = "GLUTATHIONE",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Neutrophil activation"                                    = list(pattern = "NEUTROPHIL_DEGRANULATION|NEUTROPHIL_ACTIVATION",
                                                                     database = "canonical",
                                                                     category = "molecular_mechanism", target_contrast = "pooled"),
  "Macrophage"                                                = list(pattern = "Macrophage", database = "cellmarker",
                                                                     category = "cell_type", target_contrast = NA_character_),
  "Neutrophil"                                                = list(pattern = "Neutrophil", database = "cellmarker",
                                                                     category = "cell_type", target_contrast = NA_character_)
)

mechanism_key <- vapply(manual_mechanisms, function(s) paste(s$pattern, s$database %||% ""), character(1))
search_specs <- manual_mechanisms[!duplicated(mechanism_key)]
names(search_specs) <- mechanism_key[!duplicated(mechanism_key)]

manual_recovery <- rbindlist(lapply(c("local", "ARCHS4"), function(mdl) {
  hm <- hits[model == mdl]
  candidates <- rbindlist(lapply(names(search_specs), function(key) {
    spec <- search_specs[[key]]
    d <- hm[grepl(spec$pattern, term, ignore.case = TRUE)]
    if (!is.null(spec$database)) d <- d[database %chin% spec$database]
    d <- d[order(fdr), .SD[1L], by = LV]
    if (!nrow(d)) return(NULL)
    d[, search_key := key][, .(search_key, LV, term, fdr, database)]
  }))
  setorder(candidates, fdr)
  claimed_lv <- character(0)
  claimed_key <- character(0)
  won <- rbindlist(lapply(seq_len(nrow(candidates)), function(i) {
    r <- candidates[i]
    if (r$LV %chin% claimed_lv || r$search_key %chin% claimed_key) return(NULL)
    claimed_lv <<- c(claimed_lv, r$LV)
    claimed_key <<- c(claimed_key, r$search_key)
    r
  }))
  by_key <- if (nrow(won)) split(won, won$search_key) else list()
  rbindlist(lapply(names(manual_mechanisms), function(nm) {
    key <- mechanism_key[[nm]]
    r <- by_key[[key]]
    spec <- manual_mechanisms[[nm]]
    if (is.null(r)) {
      data.table(mechanism = nm, LV = NA_character_, term = NA_character_,
                 fdr = NA_real_, database = NA_character_, model = mdl,
                 category = spec$category, target_contrast = spec$target_contrast)
    } else {
      data.table(mechanism = nm, LV = r$LV, term = r$term, fdr = r$fdr,
                 database = r$database, model = mdl,
                 category = spec$category, target_contrast = spec$target_contrast)
    }
  }))
}))
manual_recovery[, recovered := !is.na(LV)]
manual_recovery[, mechanism := factor(mechanism, levels = names(manual_mechanisms))]
setorder(manual_recovery, mechanism, model)

print(manual_recovery[, .(mechanism, model, LV, pathway = term, fdr, recovered)])
source(here("scripts", "archs4", "common.R"))

read_gmt_sets <- function(path) {
  x <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- lapply(x, function(row) unique(row[-c(1L, 2L)]))
  names(out) <- vapply(x, `[[`, "", 1L)
  out
}
read_cellmarker_sets <- function(path, sheet, term_col, gene_col) {
  x <- data.table::as.data.table(readxl::read_excel(path, sheet = sheet))
  x <- x[!is.na(get(term_col)) & !is.na(get(gene_col)),
         .(term = as.character(get(term_col)), gene = as.character(get(gene_col)))]
  split(x$gene, x$term)
}
gene_sets <- list(
  canonical  = read_gmt_sets(here("data", "pathways", "c2.cp.v2026.1.Hs.symbols.gmt")),
  hallmark   = read_gmt_sets(here("data", "pathways", "h.all.v2026.1.Hs.symbols.gmt")),
  cellmarker = read_cellmarker_sets(here("data", "pathways", "Cell_marker_Human.xlsx"),
                                     "human", "cell_name", "Symbol")
)

z_registry <- fread(file.path(prod_root, mono_dataset, "mechanism_models.tsv"))
z_mats <- list(
  ARCHS4 = read_matrix_csv(here("output", "98_final_models", "clampfull", "canonical", "archs4", "Z.csv")),
  local  = read_matrix_csv(here(z_registry[model == "local"]$z))
)

top_pct <- 0.01
recovered_rows <- manual_recovery[recovered == TRUE]

gene_loadings <- rbindlist(lapply(seq_len(nrow(recovered_rows)), function(i) {
  r <- recovered_rows[i]
  z <- z_mats[[r$model]]
  members <- intersect(gene_sets[[r$database]][[r$term]], rownames(z))
  n_top <- max(1L, ceiling(nrow(z) * top_pct))
  ord <- order(z[, r$LV], decreasing = TRUE)[seq_len(n_top)]
  data.table(dataset = mono_dataset, comparison_id = as.character(r$mechanism), model = r$model,
             rank = seq_along(ord), gene = rownames(z)[ord], loading = z[ord, r$LV],
             is_gene_set = rownames(z)[ord] %chin% members, in_top_loading_set = TRUE,
             n_gene_set_in_universe = length(members))
}))

comparisons <- copy(manual_recovery)
setnames(comparisons, "term", "gene_set")
comparisons[, `:=`(dataset = mono_dataset, comparison_id = as.character(mechanism), top_pct = top_pct)]
comparisons <- merge(
  comparisons,
  gene_loadings[, .(n_gene_set_in_top_loading = sum(is_gene_set),
                     n_gene_set_in_universe = max(n_gene_set_in_universe)),
                by = .(comparison_id, model)],
  by = c("comparison_id", "model"), all.x = TRUE)
comparisons[is.na(n_gene_set_in_top_loading), n_gene_set_in_top_loading := 0L]
comparisons[is.na(n_gene_set_in_universe), n_gene_set_in_universe := 0L]

comparison_tests <- comparisons[, {
  a <- .SD[model == "ARCHS4"]
  b <- .SD[model == "local"]
  p <- NA_real_
  arch_fraction <- NA_real_
  local_fraction <- NA_real_
  if (nrow(a) == 1L && nrow(b) == 1L && isTRUE(a$recovered) && isTRUE(b$recovered) &&
      identical(a$gene_set, b$gene_set) &&
      a$n_gene_set_in_universe > 0L && b$n_gene_set_in_universe > 0L) {
    arch_fraction <- a$n_gene_set_in_top_loading / a$n_gene_set_in_universe
    local_fraction <- b$n_gene_set_in_top_loading / b$n_gene_set_in_universe
    p <- fisher.test(matrix(c(a$n_gene_set_in_top_loading,
                              a$n_gene_set_in_universe - a$n_gene_set_in_top_loading,
                              b$n_gene_set_in_top_loading,
                              b$n_gene_set_in_universe - b$n_gene_set_in_top_loading),
                            nrow = 2, byrow = TRUE), alternative = "greater")$p.value
  }
  .(p_value = p, arch_fraction = arch_fraction, local_fraction = local_fraction)
}, by = comparison_id]
comparison_tests[, p_adj := p.adjust(p_value, method = "BH")]
comparisons <- merge(comparisons, comparison_tests, by = "comparison_id", all.x = TRUE)

comparisons[, comparison_note := ""]
comparisons[model == "local" & !recovered, comparison_note := "Not recovered"]
comparisons[recovered == TRUE & is.na(p_adj), comparison_note := "Different matched pathway"]
comparisons[recovered == TRUE & !is.na(p_adj), comparison_note := paste0(
  fifelse(arch_fraction > local_fraction, "ARCHS4 > Monocyte model",
          fifelse(arch_fraction < local_fraction, "Monocyte model > ARCHS4", "equal fraction")),
  " (FDR ", formatC(p_adj, format = "e", digits = 1), ")")]
