placenta_dataset <- "placenta_EMTAB6701"

ora_dbs <- c("canonical", "hallmark", "cellmarker", "azimuth")
hits <- rbindlist(lapply(c("local", "ARCHS4"), function(mdl) rbindlist(lapply(ora_dbs, function(db) {
  f <- file.path(prod_root, placenta_dataset, "ora", mdl, db, "enrichment.csv.gz")
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  if (!nrow(d)) return(NULL)
  d[, `:=`(model = mdl, database = db)]
  d
}))))
setnames(hits, c("ID", "p.adjust"), c("term", "fdr"), skip_absent = TRUE)
hits <- hits[fdr < 0.05 & Count >= 3]

manual_mechanisms <- list(
  "EVT invasion / EMT" =
    list(pattern = "EPITHELIAL_MESENCHYMAL_TRANSITION",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "Spiral-artery remodelling" =
    list(pattern = "CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "dNK-trophoblast HLA-C / HLA-E / HLA-G interactions" =
    list(pattern = "HLAC_ALLOTYPES_INTERACTIONS_WITH_KIR_ON_DNK_CELLS",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "dNK chemokine / immunomodulatory signaling" =
    list(pattern = "CHEMOKINE_SIGNALING",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "Immune checkpoint / maternal-fetal immune tolerance" =
    list(pattern = "CO_INHIBITION_BY_PD_1|IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
         exclude_lv = c("LV1269", "LV948", "LV90"),
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "Adenosine-mediated immunoregulation" =
    list(pattern = "ADORA2B",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "dNK1 glycolytic metabolic priming" =
    list(pattern = "KEGG_MEDICUS_REFERENCE_GLYCOLYSIS",
         category = "molecular_mechanism", target_contrast = "trophoblast_vs_other"),
  "Trophoblast" =
    list(pattern = "trophoblast", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Extravillous trophoblast (EVT)" =
    list(pattern = "Extravillous Trophoblasts", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Syncytiotrophoblast (SCT)" =
    list(pattern = "Syncytiotrophoblasts And Villous Cytotrophoblasts", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Villous cytotrophoblast (VCT)" =
    list(pattern = "Syncytiotrophoblasts And Villous Cytotrophoblasts", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Decidual NK / dNK" =
    list(pattern = "UTERINE_NATURAL_KILLER", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Decidual stromal cells" =
    list(pattern = "Stromal", prefer_cellmarker = TRUE, exclude_lv = c("LV90", "LV126"),
         category = "cell_type", target_contrast = "trophoblast_vs_other"),
  "Decidual macrophages / maternal myeloid cells" =
    list(pattern = "Macrophage", prefer_cellmarker = TRUE,
         category = "cell_type", target_contrast = "trophoblast_vs_other")
)

find_best <- function(hm, spec) {
  d <- hm[grepl(spec$pattern, term, ignore.case = TRUE)]
  if (!is.null(spec$exclude_lv)) d <- d[!LV %chin% spec$exclude_lv]
  if (!nrow(d)) return(NULL)
  if (isTRUE(spec$prefer_cellmarker)) {
    d[, rnk := fifelse(database == "cellmarker", 0L, 1L)]
    setorder(d, rnk, fdr)
    d[, rnk := NULL]
  } else setorder(d, fdr)
  d[1]
}

manual_recovery <- rbindlist(lapply(c("local", "ARCHS4"), function(mdl) {
  hm <- hits[model == mdl]
  rbindlist(lapply(names(manual_mechanisms), function(nm) {
    spec <- manual_mechanisms[[nm]]
    r <- find_best(hm, spec)
    if (is.null(r)) {
      data.table(mechanism = nm, LV = NA_character_, term = NA_character_, fdr = NA_real_,
                 database = NA_character_, model = mdl, category = spec$category,
                 target_contrast = spec$target_contrast)
    } else {
      data.table(mechanism = nm, LV = r$LV, term = r$term, fdr = r$fdr,
                 database = r$database, model = mdl, category = spec$category,
                 target_contrast = spec$target_contrast)
    }
  }))
}))
manual_recovery[, LV_label := LV]

combined <- rbindlist(lapply(c("local", "ARCHS4"), function(mdl) {
  both <- manual_recovery[model == mdl & mechanism %chin% c("Extravillous trophoblast (EVT)", "Syncytiotrophoblast (SCT)") & !is.na(LV)]
  if (!nrow(both)) {
    return(data.table(mechanism = "Trophoblast differentiation into EVT / SCT", LV = NA_character_,
                       LV_label = NA_character_, term = NA_character_, fdr = NA_real_,
                       database = NA_character_, model = mdl, category = "molecular_mechanism",
                       target_contrast = "trophoblast_vs_other"))
  }
  setorder(both, fdr)
  primary <- both[1]
  data.table(mechanism = "Trophoblast differentiation into EVT / SCT", LV = primary$LV,
             LV_label = paste(unique(both$LV), collapse = " / "),
             term = paste(sort(unique(both$term)), collapse = " / "),
             fdr = primary$fdr, database = primary$database, model = mdl,
             category = "molecular_mechanism", target_contrast = "trophoblast_vs_other")
}))
manual_recovery <- rbind(manual_recovery, combined)

row_order <- c(
  "Trophoblast differentiation into EVT / SCT", "EVT invasion / EMT", "Spiral-artery remodelling",
  "dNK-trophoblast HLA-C / HLA-E / HLA-G interactions", "dNK chemokine / immunomodulatory signaling",
  "Immune checkpoint / maternal-fetal immune tolerance", "Adenosine-mediated immunoregulation",
  "dNK1 glycolytic metabolic priming", "Trophoblast", "Extravillous trophoblast (EVT)",
  "Syncytiotrophoblast (SCT)", "Villous cytotrophoblast (VCT)", "Decidual NK / dNK",
  "Decidual stromal cells", "Decidual macrophages / maternal myeloid cells")
manual_recovery[, recovered := !is.na(LV)]
manual_recovery[, mechanism := factor(mechanism, levels = row_order)]
setorder(manual_recovery, mechanism, model)

print(manual_recovery[, .(mechanism, model, LV_label, pathway = term, fdr, recovered)])

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
  azimuth    = read_gmt_sets(here("data", "pathways", "Azimuth_2023.txt")),
  cellmarker = read_cellmarker_sets(here("data", "pathways", "Cell_marker_Human.xlsx"),
                                     "human", "cell_name", "Symbol")
)

z_registry <- fread(file.path(prod_root, placenta_dataset, "mechanism_models.tsv"))
z_mats <- list(
  ARCHS4 = read_matrix_csv(here("output", "98_final_models", "clampfull", "canonical", "archs4", "Z.csv")),
  local  = read_matrix_csv(here(z_registry[model == "local"]$z))
)

top_pct <- 0.01
recovered_rows <- manual_recovery[recovered == TRUE]

gene_loadings <- rbindlist(lapply(seq_len(nrow(recovered_rows)), function(i) {
  r <- recovered_rows[i]
  z <- z_mats[[r$model]]
  members <- intersect(unique(unlist(gene_sets[[r$database]][strsplit(r$term, " / ", fixed = TRUE)[[1]]])), rownames(z))
  n_top <- max(1L, ceiling(nrow(z) * top_pct))
  ord <- order(z[, r$LV], decreasing = TRUE)[seq_len(n_top)]
  data.table(dataset = placenta_dataset, comparison_id = as.character(r$mechanism), model = r$model,
             rank = seq_along(ord), gene = rownames(z)[ord], loading = z[ord, r$LV],
             is_gene_set = rownames(z)[ord] %chin% members, in_top_loading_set = TRUE,
             n_gene_set_in_universe = length(members))
}))

comparisons <- copy(manual_recovery)
setnames(comparisons, "term", "gene_set")
comparisons[, `:=`(dataset = placenta_dataset, comparison_id = as.character(mechanism), top_pct = top_pct)]
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
  fifelse(arch_fraction > local_fraction, "ARCHS4 > Placenta model",
          fifelse(arch_fraction < local_fraction, "Placenta model > ARCHS4", "equal fraction")),
  " (FDR ", formatC(p_adj, format = "e", digits = 1), ")")]
