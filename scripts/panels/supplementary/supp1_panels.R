# Original analysis and provenance guards retained; layout/export lives in supp1.R.
panel_input <- lapply(jsonlite::fromJSON(here::here("scripts/panels/supplementary/supp1_inputs.json")),here::here)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(yaml)
  library(here)
  library(grid)
  library(ragg)
  library(svglite)
})

FONT_FAMILY <- "Helvetica"

FS_TAG          <- 7
FS_TITLE        <- 7
FS_SUBTITLE     <- 6.5
FS_AXIS_TITLE   <- 7
FS_AXIS_TEXT    <- 6.5
FS_LEGEND       <- 5.5
FS_LEGEND_TITLE <- 6.0
FS_STAT         <- 6.0
FS_HEAT_LABEL   <- 5.5
FS_HEAT_VALUE   <- 5.5
FS_CELL_VALUE   <- 5
FS_DENSE_VALUE  <- 5
FS_A_MEAN       <- 5
FS_HEAT_LABEL_H <- 5

LINE_W <- 0.3
GRID_W <- 0.2

pt_mm <- function(pt) pt / .pt

theme_nature_methods <- function(base_size = FS_AXIS_TEXT,
                                 grid = c("none", "y", "x", "both")) {
  grid <- match.arg(grid)
  th <- theme_classic(base_size = base_size, base_family = FONT_FAMILY) %+replace%
    theme(
      plot.title        = element_text(size = FS_TITLE, face = "bold", hjust = 0.5,
                                       margin = margin(b = 1)),
      plot.subtitle     = element_text(size = FS_SUBTITLE, face = "plain", hjust = 0.5,
                                       lineheight = 0.95, margin = margin(b = 1)),
      plot.tag          = element_text(size = FS_TAG, face = "bold", family = FONT_FAMILY),
      plot.tag.position = "topleft",
      plot.tag.location = "margin",
      axis.line         = element_line(linewidth = LINE_W, colour = "black"),
      axis.ticks        = element_line(linewidth = LINE_W, colour = "black"),
      axis.ticks.length = unit(0.9, "pt"),
      axis.text         = element_text(size = base_size, colour = "black"),
      axis.text.x       = element_text(margin = margin(t = 0.8)),
      axis.text.y       = element_text(hjust = 1, margin = margin(r = 0.8)),
      axis.title        = element_text(size = FS_AXIS_TITLE, colour = "black"),
      axis.title.x      = element_text(margin = margin(t = 1)),
      axis.title.y      = element_text(angle = 90, margin = margin(r = 1)),
      legend.text       = element_text(size = FS_LEGEND),
      legend.title      = element_text(size = FS_LEGEND_TITLE),
      legend.key.size   = unit(2, "mm"),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      strip.text        = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 1)),
      strip.background  = element_blank(),
      panel.background  = element_blank(),
      panel.grid        = element_blank(),
      plot.background   = element_blank(),
      plot.margin       = margin(1, 1, 1, 1, "mm")
    )
  if (grid %in% c("y", "both"))
    th <- th + theme(panel.grid.major.y = element_line(colour = "grey88", linewidth = GRID_W))
  if (grid %in% c("x", "both"))
    th <- th + theme(panel.grid.major.x = element_line(colour = "grey88", linewidth = GRID_W))
  th
}

theme_nature_heatmap <- function(x_angle = 45, label_size = FS_HEAT_LABEL) {
  theme_nature_methods() %+replace%
    theme(
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = element_text(size = label_size, angle = x_angle,
                                 hjust = if (x_angle == 0) 0.5 else 1,
                                 vjust = if (x_angle == 90) 0.5 else 1,
                                 colour = "black",
                                 margin = margin(t = 0.5), lineheight = 0.9),
      axis.text.y = element_text(size = label_size, hjust = 1, colour = "black",
                                 margin = margin(r = 0.5), lineheight = 0.9),
      plot.margin = margin(1, 1, 1, 1, "mm")
    )
}

add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag          = element_text(size = FS_TAG, face = "bold",
                                           family = FONT_FAMILY),
          plot.tag.position = "topleft",
          plot.tag.location = "margin")
}

add_overlay_tag <- function(p, tag) {
  tagged <- ggdraw(p) +
    draw_label(tag, x = 0, y = 1, hjust = 0, vjust = 1,
               size = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY)
  tagged
}

DATASET_LABELS <- c(
  Brain_Mathys2023 = "Brain Mathys",
  Brain_Xiong2023  = "Brain Xiong",
  Heart_Datar2026  = "Heart Datar",
  PBMC_1k1k        = "PBMC 1k1k",
  PBMC_Perez2022   = "PBMC Perez",
  Lung_Sikkema2023 = "Lung Sikkema"
)

TISSUE_MAP <- c(
  Brain_Mathys2023 = "Brain", Brain_Xiong2023 = "Brain",
  Heart_Datar2026  = "Heart",
  Lung_Sikkema2023 = "Lung",
  PBMC_1k1k        = "PBMC", PBMC_Perez2022 = "PBMC"
)

DATASETS_ROW1 <- c("Heart_Datar2026", "PBMC_1k1k", "Lung_Sikkema2023")
DATASETS_ROW2 <- c("Brain_Mathys2023", "Brain_Xiong2023", "PBMC_Perez2022")


TISSUE_ORDER <- c("Brain", "Heart", "Lung", "PBMC")

abbreviate_ct <- function(x) {
  x <- gsub("Oligodendrocyte [Pp]rogenitor [Cc]ells?", "OPC", x)
  x <- gsub("Oligodendrocyte [Pp]recursor [Cc]ells?", "OPC", x)
  x
}

CT_SHORT <- c(
  "Oligodendrocyte Precursor Cells" = "OPC",
  "Plasmacytoid Dendritic Cells"    = "pDC",
  "LymphaticEndothelial"            = "Lymphatic EC",
  "Alveolar epithelium"             = "Alveolar ep.",
  "Airway epithelium"               = "Airway ep.",
  "Excitatory neurons"              = "Excitatory",
  "Inhibitory neurons"              = "Inhibitory",
  "Endothelial cells"               = "Endothelial",
  "Fibroblast lineage"              = "Fibroblast",
  "Submucosal Gland"                = "Submucosal",
  "Oligodendrocytes"                = "Oligodendro.",
  "CD14+ Monocytes"                 = "CD14+ Mono",
  "CD16+ Monocytes"                 = "CD16+ Mono",
  "Dendritic cells"                 = "Dendritic",
  "Plasma B cells"                  = "Plasma B",
  "Myeloid cells"                   = "Myeloid",
  "Malignant cells"                 = "Malignant",
  "Vascular cells"                  = "Vascular",
  "Cancer cells"                    = "Cancer",
  "Mast cells"                      = "Mast",
  "Blood vessels"                   = "Blood vessel"
)
short_ct <- function(x) {
  x <- abbreviate_ct(x)
  ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
}

wrap_label <- function(x, width = 14) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

fmt_corr_cell <- function(x) sub("0.", ".", sprintf("%.2f", x), fixed = TRUE)

shorten_gtex <- function(x) {
  x <- gsub("Adipose - Subcutaneous",                    "Adipose - Subcut.", x)
  x <- gsub("Adipose - Visceral \\(Omentum\\)",             "Adipose - Visceral", x)
  x <- gsub("Artery - ",                                 "Artery - ", x)
  x <- gsub("Brain - Amygdala",                          "Brain - Amygdala", x)
  x <- gsub("Brain - Anterior cingulate cortex \\(BA24\\)", "Brain - ACC (BA24)", x)
  x <- gsub("Brain - Caudate \\(basal ganglia\\)",         "Brain - Caudate", x)
  x <- gsub("Brain - Cerebellar Hemisphere",             "Brain - Cereb. hem.", x)
  x <- gsub("Brain - Frontal Cortex \\(BA9\\)",            "Brain - FC (BA9)", x)
  x <- gsub("Brain - Nucleus accumbens \\(basal ganglia\\)", "Brain - NAc", x)
  x <- gsub("Brain - Putamen \\(basal ganglia\\)",         "Brain - Putamen", x)
  x <- gsub("Brain - Spinal cord \\(cervical c-1\\)",      "Brain - Spinal cord", x)
  x <- gsub("Brain - Substantia nigra",                  "Brain - Subst. nigra", x)
  x <- gsub("Breast - Mammary Tissue",                   "Breast - Mammary", x)
  x <- gsub("Cells - Cultured fibroblasts",              "Cells - Fibroblasts", x)
  x <- gsub("Cells - EBV-transformed lymphocytes",       "Cells - EBV lymph.", x)
  x <- gsub("Esophagus - Gastroesophageal Junction",     "Esoph. - GE junction", x)
  x <- gsub("Esophagus - Mucosa",                        "Esoph. - Mucosa", x)
  x <- gsub("Esophagus - Muscularis",                    "Esoph. - Muscularis", x)
  x <- gsub("Heart - Atrial Appendage",                  "Heart - Atrial app.", x)
  x <- gsub("Minor Salivary Gland",                      "Minor saliv. gland", x)
  x <- gsub("Skin - Not Sun Exposed \\(Suprapubic\\)",     "Skin - Not sun exp.", x)
  x <- gsub("Skin - Sun Exposed \\(Lower leg\\)",          "Skin - Sun exp.", x)
  x <- gsub("Small Intestine - Terminal Ileum",          "Sm. intestine - Ileum", x)
  trimws(x)
}

fmt_q_compact <- function(q) {
  if (is.na(q)) return('"n.a."')
  if (q >= 0.001) return(sprintf('"%.3f"', q))
  e_str    <- formatC(q, format = "e", digits = 1)
  parts    <- strsplit(e_str, "e")[[1]]
  sprintf('%s%%*%%10^{%d}', trimws(parts[1]), as.integer(parts[2]))
}

assign_bracket_tiers <- function(comp_df, xpos) {
  comp_ord <- copy(as.data.table(comp_df))
  comp_ord[, x1 := xpos[a]]
  comp_ord[, x2 := xpos[b]]
  comp_ord[, left  := pmin(x1, x2)]
  comp_ord[, right := pmax(x1, x2)]
  comp_ord[, span  := right - left]
  setorder(comp_ord, span, q)

  levels_used <- list()
  comp_ord[, tier := 0L]
  for (i in seq_len(nrow(comp_ord))) {
    left  <- comp_ord$left[i]; right <- comp_ord$right[i]; tier <- 1L
    repeat {
      current  <- if (tier <= length(levels_used)) levels_used[[tier]] else NULL
      overlaps <- !is.null(current) && any(vapply(current, function(iv) {
        !(right < iv[1] || left > iv[2])
      }, logical(1)))
      if (!overlaps) break
      tier <- tier + 1L
    }
    comp_ord$tier[i] <- tier
    prior <- if (tier <= length(levels_used)) levels_used[[tier]] else list()
    levels_used[[tier]] <- c(prior, list(c(left, right)))
  }
  comp_ord
}

add_brackets <- function(p, comp_ord, y_base, y_step, h, size = pt_mm(FS_STAT)) {
  for (i in seq_len(nrow(comp_ord))) {
    r <- comp_ord[i, ]
    y <- y_base + (r$tier - 1) * y_step
    p <- p +
      annotate("segment", x = r$left,  xend = r$left,  y = y,     yend = y + h, linewidth = 0.2) +
      annotate("segment", x = r$left,  xend = r$right, y = y + h, yend = y + h, linewidth = 0.2) +
      annotate("segment", x = r$right, xend = r$right, y = y,     yend = y + h, linewidth = 0.2) +
      annotate("text", x = (r$left + r$right) / 2, y = y + h + 0.005,
               label = fmt_q_compact(r$q), size = size, vjust = 0,
               parse = TRUE, family = FONT_FAMILY)
  }
  p
}

cfg <- yaml::read_yaml(here("config.yaml"))
MODEL_COLORS_RAW <- unlist(cfg$MODEL_COLORS)
names(MODEL_COLORS_RAW)[names(MODEL_COLORS_RAW) == "GenomicSuperSignature"] <- "GSSig"

ct_labels_df <- read.csv(here("data", "pseudobulk", "cell_type_labels.csv"), stringsAsFactors = FALSE)
CT_LABELS <- setNames(ct_labels_df$label, ct_labels_df$cell_type)
ct_label <- function(x) ifelse(x %in% names(CT_LABELS), CT_LABELS[x], x)

GREEN_SCALE      <- unlist(cfg$GREEN_SCALE)
DIVERGING_COLORS <- unlist(cfg$DIVERGING_COLORS)
DIVERGING_VALUES <- as.numeric(cfg$DIVERGING_VALUES)

TILE_TEXT <- c(`TRUE` = "white", `FALSE` = "black")

long <- fread(panel_input[["benchmark_long"]])
stopifnot(all(c("dataset", "method", "truth", "cor") %in% names(long)))

long_box <- long[truth == "v0" & !is.na(cor)]
long_box[, tissue := TISSUE_MAP[dataset]]
stopifnot(!anyNA(long_box$tissue))

METHOD_ORDER  <- long_box[, .(m = mean(cor, na.rm = TRUE)), by = method][order(-m), as.character(method)]
METHOD_COLORS <- MODEL_COLORS_RAW[METHOD_ORDER]
METHOD_COLORS[is.na(METHOD_COLORS)] <- "grey70"
names(METHOD_COLORS) <- METHOD_ORDER

long_box[, method := factor(method, levels = METHOD_ORDER)]
mean_by_tissue <- long_box[, .(mean_cor = mean(cor, na.rm = TRUE),
                               max_cor  = max(cor, na.rm = TRUE)), by = .(tissue, method)]

A_Y_MAX <- 1.16

make_tissue_box <- function(tis, show_x, show_y) {
  d <- long_box[tissue == tis]
  m <- mean_by_tissue[tissue == tis]
  p <- ggplot(d, aes(method, cor, fill = method)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, colour = "black",
                 linewidth = LINE_W, alpha = 0.85) +
    geom_jitter(width = 0.10, size = 0.3, shape = 21, fill = "white",
                colour = "#333333", stroke = 0.12, alpha = 0.7) +
    geom_point(data = m, aes(x = method, y = mean_cor), shape = 23, size = 1.1,
               fill = "white", colour = "black", stroke = 0.3, inherit.aes = FALSE) +
    geom_text(data = m, aes(x = method, y = 1.05, label = sprintf("%.3f", mean_cor)),
              size = pt_mm(FS_A_MEAN), colour = "black", inherit.aes = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                       expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(0, A_Y_MAX), clip = "on") +
    scale_fill_manual(values = METHOD_COLORS, na.value = "grey70", drop = FALSE) +
    labs(x = NULL, y = "Max Pearson r\nper cell type", title = tis) +
    theme_nature_methods(grid = "none") +
    theme(legend.position = "none")
  p <- if (show_x) {
    p + theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1,
                                         size = FS_AXIS_TEXT, margin = margin(t = 0.8)))
  } else {
    p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

plot_A_brain <- make_tissue_box("Brain", show_x = FALSE, show_y = TRUE)
plot_A_heart <- make_tissue_box("Heart", show_x = FALSE, show_y = FALSE)
plot_A_lung  <- make_tissue_box("Lung",  show_x = TRUE,  show_y = TRUE)
plot_A_pbmc  <- make_tissue_box("PBMC",  show_x = TRUE,  show_y = FALSE)

A_grid <- plot_grid(plot_A_brain, plot_A_heart, plot_A_lung, plot_A_pbmc,ncol=2)

panel_A <- add_overlay_tag(A_grid, "A")

options(repr.plot.width = 6, repr.plot.height = 3)

win_rate <- fread(panel_input[["bootstrap"]])
stopifnot(all(c("method", "win_rate") %in% names(win_rate)))

setorder(win_rate, win_rate)
win_rate[, method := factor(method, levels = method)]

panel_B <- ggplot(win_rate, aes(x = win_rate, y = method, fill = method)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * win_rate)),
            hjust = -0.2, size = pt_mm(FS_AXIS_TEXT), colour = "black") +
  scale_fill_manual(values = MODEL_COLORS_RAW, na.value = "grey70") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = function(x) paste0(100 * x, "%"),
                     expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = NULL) +
  theme_nature_methods(grid = "none") +
  theme(legend.position = "none",
        axis.line.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.text.x   = element_text(margin = margin(t = 0)),
        axis.ticks.length.x = unit(0.6, "pt"),
        plot.margin   = margin(1, 1, 0.3, 1, "mm"))

panel_B <- add_tag(panel_B, "B")

options(repr.plot.width = 2, repr.plot.height = 3)

heatmap_long <- fread(panel_input[["heatmap_long"]])
pseudobulk_recovery <- fread(panel_input[["pseudobulk_recovery"]])
pseudobulk_overall <- fread(panel_input[["pseudobulk_overall"]])
stopifnot(all(c("dataset", "row_cell_type", "col_cell_type", "pct") %in% names(heatmap_long)),
          all(c("dataset", "cell_type", "n_top", "diag_count", "recovery_pct",
                "assigned_lv", "purity_ratio", "purity_lift") %in% names(pseudobulk_recovery)),
          all(c("metric", "value") %in% names(pseudobulk_overall)),
          grepl("/03_b_matrix_singlecell/", panel_input[["heatmap_long"]], fixed=TRUE),
          grepl("/03_b_matrix_singlecell/", panel_input[["pseudobulk_recovery"]], fixed=TRUE),
          grepl("/03_b_matrix_singlecell/", panel_input[["pseudobulk_overall"]], fixed=TRUE),
          nrow(pseudobulk_recovery) == 47L)

# Hard provenance guard for Supplementary Fig. 1e. These values belong to
# the pseudobulk projection pipeline and must not be replaced by the nearby
# donor-bulk values used in Supplementary Fig. 2a.
assert_close <- function(actual, expected, tolerance, label) {
  if (length(actual) != 1L || !is.finite(actual) || abs(actual - expected) > tolerance)
    stop(sprintf("%s provenance check failed: observed %.6f, expected %.6f +/- %.6f",
                 label, actual, expected, tolerance))
}
pb_metric <- setNames(pseudobulk_overall$value, pseudobulk_overall$metric)
assert_close(pb_metric[["global_recovery_pct_pooled"]], 70.269472, 0.05, "pseudobulk pooled purity (%)")
assert_close(pb_metric[["mean_recovery_pct_by_dataset"]], 76.801048, 0.05, "pseudobulk mean dataset purity (%)")
assert_close(pb_metric[["mean_purity_ratio"]], 0.772615, 0.001, "pseudobulk mean purity ratio")
assert_close(pb_metric[["mean_purity_lift"]], 19.060873, 0.05, "pseudobulk mean purity lift")
assert_close(sum(pseudobulk_recovery$diag_count) / sum(pseudobulk_recovery$n_top) * 100,
             pb_metric[["global_recovery_pct_pooled"]], 1e-6, "pseudobulk pooled purity recomputation")
perez_b <- pseudobulk_recovery[dataset == "PBMC_Perez2022" & cell_type == "B_cell"]
stopifnot(nrow(perez_b) == 1L, perez_b$assigned_lv == "LV10")
stopifnot(setequal(unique(paste(heatmap_long$dataset, heatmap_long$col_cell_type)),
                   paste(pseudobulk_recovery$dataset, pseudobulk_recovery$cell_type)))

PURITY_LABEL_MIN <- 15

make_purity_panel <- function(ds) {
  d <- heatmap_long[dataset == ds]
  stopifnot(nrow(d) > 0)

  ct_order  <- unique(d$row_cell_type)
  ct_order  <- ct_order[order(ct_label(ct_order))]
  labs_ord  <- short_ct(ct_label(ct_order))

  d[, row_label   := factor(short_ct(ct_label(row_cell_type)), levels = labs_ord)]
  d[, col_label   := factor(short_ct(ct_label(col_cell_type)), levels = labs_ord)]
  d[, is_diagonal := row_cell_type == col_cell_type]
  d[, show_label  := is_diagonal | pct >= PURITY_LABEL_MIN]

  ggplot(d, aes(x = col_label, y = row_label, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_tile(data = d[is_diagonal == TRUE], fill = NA, colour = "black", linewidth = 0.4) +
    geom_text(data = d[show_label == TRUE],
              aes(label = sprintf("%.0f", pct), colour = pct >= 65),
              size = pt_mm(FS_CELL_VALUE), show.legend = FALSE) +
    scale_fill_gradientn(colours = GREEN_SCALE, limits = c(0, 100),
                         name = "Top 1% purity (%)",
                         guide = guide_colourbar(barwidth = unit(20, "mm"),
                                                 barheight = unit(1.8, "mm"),
                                                 title.position = "left",
                                                 title.vjust = 1)) +
    scale_colour_manual(values = TILE_TEXT, guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = DATASET_LABELS[ds]) +
    theme_nature_heatmap(x_angle = 45, label_size = FS_HEAT_LABEL_H)
}

C_row1_panels <- lapply(DATASETS_ROW1, make_purity_panel)
C_row2_panels <- lapply(DATASETS_ROW2, make_purity_panel)

C_N1 <- vapply(DATASETS_ROW1, function(ds) uniqueN(heatmap_long[dataset == ds, row_cell_type]), numeric(1))
C_N2 <- vapply(DATASETS_ROW2, function(ds) uniqueN(heatmap_long[dataset == ds, row_cell_type]), numeric(1))

C_LEG <- max(C_N2)

C_row1 <- wrap_plots(C_row1_panels, nrow = 1, widths = C_N1) &
  theme(legend.position = "none")
C_row2 <- wrap_plots(c(C_row2_panels, list(guide_area())),
                     nrow = 1, widths = c(C_N2, C_LEG),
                     guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "left")

panel_C_core <- wrap_plots(C_row1, C_row2, ncol = 1, heights = c(1.35, 1))
panel_C <- wrap_elements(full =
  ggdraw() +
    draw_plot(panel_C_core, x = 0.042, y = 0.038, width = 0.958, height = 0.962) +
    draw_label("LV assigned to cell type", x = 0.52, y = 0.002,
               hjust = 0.5, vjust = 0, size = FS_AXIS_TITLE,
               fontfamily = FONT_FAMILY) +
    draw_label("Cell type of top 1% projected cells", x = 0.006, y = 0.52,
               angle = 90, hjust = 0.5, vjust = 0, size = FS_AXIS_TITLE,
               fontfamily = FONT_FAMILY) +
    draw_label("E", x = 0, y = 1, hjust = 0, vjust = 1,
               size = FS_TAG, fontface = "bold", fontfamily = FONT_FAMILY)
)

options(repr.plot.width = 7.2, repr.plot.height = 2.5)

corr_full          <- fread(panel_input[["corr_full"]])
assignments        <- fread(panel_input[["assignments"]])
stopifnot(all(c("dataset", "LV", "cell_type", "cor") %in% names(corr_full)))
stopifnot(all(c("dataset", "LV", "cell_type") %in% names(assignments)))

CORR_LABEL_MIN <- 0.6

make_corr_panel <- function(ds) {
  d        <- corr_full[dataset == ds & !is.na(cell_type) & cell_type != "NA"]
  assigned <- assignments[dataset == ds & !is.na(cell_type) & cell_type != "NA"]
  stopifnot(nrow(d) > 0, nrow(assigned) > 0)
  d <- d[LV %in% assigned$LV & cell_type %in% assigned$cell_type]

  ord      <- order(assigned$cell_type)
  lv_order <- unique(assigned$LV[ord])
  ct_order <- unique(short_ct(ct_label(assigned$cell_type[ord])))

  assigned_key <- paste(assigned$LV, assigned$cell_type)
  d[, cell_type_label := factor(short_ct(ct_label(cell_type)), levels = ct_order)]
  d[, LV              := factor(LV, levels = lv_order)]
  d[, is_assigned     := paste(LV, cell_type) %in% assigned_key]
  d[, show_label      := is_assigned | abs(cor) >= CORR_LABEL_MIN]
  value_size <- if (ds %in% DATASETS_ROW1) FS_DENSE_VALUE else FS_CELL_VALUE

  ggplot(d, aes(x = LV, y = cell_type_label, fill = cor)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_tile(data = d[is_assigned == TRUE], fill = NA, colour = "black", linewidth = 0.4) +
    geom_text(data = d[show_label == TRUE],
              aes(label = fmt_corr_cell(cor), colour = abs(cor) >= 0.8),
              size = pt_mm(value_size), show.legend = FALSE) +
    scale_fill_gradientn(colours = DIVERGING_COLORS, values = DIVERGING_VALUES,
                         limits = c(-1, 1), name = "Pearson r",
                         guide = guide_colourbar(barwidth = unit(20, "mm"),
                                                 barheight = unit(1.8, "mm"),
                                                 title.position = "left",
                                                 title.vjust = 1)) +
    scale_colour_manual(values = TILE_TEXT, guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = DATASET_LABELS[ds]) +
    theme_nature_heatmap(x_angle = 45, label_size = FS_HEAT_LABEL_H) +
    theme(plot.title    = element_text(size = FS_SUBTITLE + 0.5, face = "bold",
                                       hjust = 0.5, margin = margin(b = 0.5)),
          plot.subtitle = element_blank())
}

D_row1_panels <- lapply(DATASETS_ROW1, make_corr_panel)
D_row2_panels <- lapply(DATASETS_ROW2, make_corr_panel)
D_row1_panels[[1]] <- add_tag(D_row1_panels[[1]], "D")

D_ncat <- function(ds) uniqueN(assignments[dataset == ds & !is.na(cell_type) & cell_type != "NA", LV])
D_N1 <- vapply(DATASETS_ROW1, D_ncat, numeric(1))
D_N2 <- vapply(DATASETS_ROW2, D_ncat, numeric(1))
D_row1 <- wrap_plots(D_row1_panels, nrow = 1, widths = D_N1)
D_row2 <- wrap_plots(D_row2_panels, nrow = 1, widths = D_N2)

panel_D <- wrap_plots(D_row1, D_row2, ncol = 1, heights = c(1.35, 1))

options(repr.plot.width = 5, repr.plot.height = 2.5)

related_corr <- fread(panel_input[["related_corr"]])
stopifnot(all(c("group_id", "LV", "cell_type", "assigned_cell_type", "cor") %in% names(related_corr)))

E_GROUP_ORDER <- sort(unique(related_corr$group_id))

pair_axis_label <- function(x) {
  x <- short_ct(ct_label(x))
  x <- sub(" cells$", "", x)
  x <- sub("^VentricularCM$", "Ventric. CM", x)
  x <- sub("^AtrialCM$", "Atrial CM", x)
  x
}

pair_y_axis_label <- function(x) {
  x <- pair_axis_label(x)
  x[x == "B"] <- "B cell"
  x[x == "T"] <- "T cell"
  x
}

make_pair_panel <- function(gid) {
  d         <- related_corr[group_id == gid]
  members   <- unique(d$assigned_cell_type)
  member_lv <- d$LV[match(members, d$assigned_cell_type)]
  d         <- d[cell_type %in% members]

  short_title <- short_ct(ct_label(members))
  short_axis  <- pair_axis_label(members)
  row_levels  <- wrap_label(pair_y_axis_label(members), width = 10)
  col_levels <- paste0(wrap_label(short_axis, width = 10), "\n", member_lv)

  d[, row_label := factor(wrap_label(pair_y_axis_label(cell_type), width = 10),
                          levels = row_levels)]
  d[, col_label := factor(paste0(wrap_label(pair_axis_label(assigned_cell_type), width = 10),
                                 "\n", LV),
                          levels = col_levels)]

  panel_title <- paste(strwrap(paste(short_title, collapse = " vs. "), width = 16), collapse = "\n")

  ggplot(d, aes(x = col_label, y = row_label, fill = cor)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", cor), colour = abs(cor) >= 0.8),
              size = pt_mm(FS_HEAT_VALUE), show.legend = FALSE) +
    scale_fill_gradientn(colours = DIVERGING_COLORS, values = DIVERGING_VALUES,
                         limits = c(-1, 1), name = "Pearson r",
                         guide = guide_colourbar(barwidth = unit(20, "mm"),
                                                 barheight = unit(1.8, "mm"),
                                                 title.position = "left",
                                                 title.vjust = 1)) +
    scale_colour_manual(values = TILE_TEXT, guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = panel_title) +
    theme_nature_heatmap(x_angle = 0, label_size = FS_HEAT_LABEL_H) +
    theme(plot.title = element_text(size = FS_SUBTITLE, face = "bold", hjust = 0.5,
                                    lineheight = 0.95, margin = margin(b = 1)))
}

E_panels <- lapply(E_GROUP_ORDER, make_pair_panel)
E_panels[[1]] <- add_tag(E_panels[[1]], "E")

panel_E <- wrap_plots(E_panels, ncol = 2)

options(repr.plot.width = 2.2, repr.plot.height = 2.5)

if (FALSE) {
gene_frac_ari         <- fread(panel_input[["gene_fraction_ari_data"]])
gene_frac_comparisons <- fread(panel_input[["gene_fraction_ari_comparisons"]])
stopifnot(all(c("fraction", "ari") %in% names(gene_frac_ari)))
stopifnot(all(c("a", "b", "q") %in% names(gene_frac_comparisons)))

FRACTION_LEVELS <- c("100%", "75%", "50%", "25%", "10%", "5%", "1%")
FRACTION_TICKS  <- sub("%$", "", FRACTION_LEVELS)
gene_frac_ari[, fraction := factor(fraction, levels = FRACTION_LEVELS)]
stopifnot(!anyNA(gene_frac_ari$fraction))

FRACTION_COLORS <- unlist(cfg$RNASEQ_FRACTION_COLORS)[FRACTION_LEVELS]

mean_df_genefrac <- gene_frac_ari[, .(mean_ari = mean(ari, na.rm = TRUE),
                                      max_ari  = max(ari, na.rm = TRUE)), by = fraction]

GENEFRAC_SHOWN <- c("75%", "50%", "10%", "1%")
comp_genefrac  <- gene_frac_comparisons[a == "100%" & b %in% GENEFRAC_SHOWN]
stopifnot(nrow(comp_genefrac) == length(GENEFRAC_SHOWN))

xpos_genefrac  <- setNames(seq_along(FRACTION_LEVELS), FRACTION_LEVELS)
comp_genefrac  <- assign_bracket_tiers(comp_genefrac, xpos_genefrac)

y_base_genefrac <- max(mean_df_genefrac$max_ari) + 0.04
y_step_genefrac <- 0.090
h_genefrac      <- 0.012
y_max_genefrac  <- y_base_genefrac + (max(comp_genefrac$tier) - 1) * y_step_genefrac +
                   h_genefrac + 0.085

panel_F <- ggplot(gene_frac_ari, aes(fraction, ari, fill = fraction)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black",
               linewidth = 0.25, alpha = 0.85) +
  geom_jitter(width = 0.10, size = 0.3, shape = 21, fill = "white",
              colour = "#333333", stroke = 0.12, alpha = 0.75) +
  geom_point(data = mean_df_genefrac, aes(x = fraction, y = mean_ari), shape = 23,
             size = 1.1, fill = "white", colour = "black", stroke = 0.3,
             inherit.aes = FALSE) +
  scale_fill_manual(values = FRACTION_COLORS, na.value = "grey70") +
  scale_x_discrete(labels = setNames(FRACTION_TICKS, FRACTION_LEVELS)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(ylim = c(0, y_max_genefrac), clip = "on") +
  labs(x = "Genes used for RNA-Seq (%)", y = "Adjusted Rand Index (ARI)") +
  theme_nature_methods(grid = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = FS_LEGEND_TITLE),
        axis.title.y = element_text(angle = 90, margin = margin(r = 0.3)),
        plot.margin = margin(1, 1, 1, 0.3, "mm"))

panel_F <- add_brackets(panel_F, comp_genefrac,
                        y_base = y_base_genefrac, y_step = y_step_genefrac, h = h_genefrac,
                        size = pt_mm(FS_HEAT_LABEL_H))
panel_F <- add_overlay_tag(panel_F, "F")

options(repr.plot.width = 3.6, repr.plot.height = 2.2)

}

if (FALSE) {
ari_data        <- fread(panel_input[["ari_data"]])
ari_comparisons <- fread(panel_input[["ari_comparisons"]])
stopifnot(all(c("method", "ari") %in% names(ari_data)))
stopifnot(all(c("a", "b", "q") %in% names(ari_comparisons)))

ari_data[, method := fifelse(method == "GenomicSuperSignature", "GSSig", method)]
ari_comparisons[, a := fifelse(a == "GenomicSuperSignature", "GSSig", a)]
ari_comparisons[, b := fifelse(b == "GenomicSuperSignature", "GSSig", b)]

ari_method_means <- ari_data[, .(m = mean(ari, na.rm = TRUE)), by = method]
G_PINNED_ORDER <- c("CLAMPfull", "PLIER", "CLAMPbase", "NMF")
G_PINNED_ORDER <- G_PINNED_ORDER[G_PINNED_ORDER %in% ari_method_means$method]
method_order_ari <- c(
  G_PINNED_ORDER,
  ari_method_means[!method %in% G_PINNED_ORDER][order(-m), as.character(method)]
)
ari_data[, method := factor(method, levels = method_order_ari)]

METHOD_COLORS_ARI <- MODEL_COLORS_RAW[method_order_ari]
METHOD_COLORS_ARI[is.na(METHOD_COLORS_ARI)] <- "grey70"
names(METHOD_COLORS_ARI) <- method_order_ari

mean_df_ari <- ari_data[, .(mean_ari = mean(ari, na.rm = TRUE),
                            max_ari  = max(ari, na.rm = TRUE)), by = method]
mean_df_ari[, method := factor(method, levels = method_order_ari)]

comp_ari <- copy(ari_comparisons)
xpos_ari <- setNames(seq_along(method_order_ari), method_order_ari)
comp_ari <- assign_bracket_tiers(comp_ari, xpos_ari)

y_base_ari <- max(mean_df_ari$max_ari) + 0.05
y_step_ari <- 0.080
h_ari      <- 0.015
y_max_ari  <- y_base_ari + (max(comp_ari$tier) - 1) * y_step_ari + h_ari + 0.085

panel_G <- ggplot(ari_data, aes(method, ari, fill = method)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black",
               linewidth = 0.25, alpha = 0.85) +
  geom_jitter(width = 0.10, size = 0.3, shape = 21, fill = "white",
              colour = "#333333", stroke = 0.12, alpha = 0.75) +
  geom_point(data = mean_df_ari, aes(x = method, y = mean_ari), shape = 23,
             size = 1.1, fill = "white", colour = "black", stroke = 0.3,
             inherit.aes = FALSE) +
  scale_fill_manual(values = METHOD_COLORS_ARI, na.value = "grey70") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1.0"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(ylim = c(0, y_max_ari), clip = "on") +
  labs(x = NULL, y = "Adjusted Rand Index (ARI)") +
  theme_nature_methods(grid = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = FS_AXIS_TEXT),
        axis.title.y = element_text(angle = 90, margin = margin(r = 0.3)),
        plot.margin = margin(1, 1, 1, 0.3, "mm"))

panel_G <- add_brackets(panel_G, comp_ari,
                        y_base = y_base_ari, y_step = y_step_ari, h = h_ari,
                        size = pt_mm(FS_HEAT_LABEL_H))
panel_G <- add_overlay_tag(panel_G, "G")

options(repr.plot.width = 3.6, repr.plot.height = 2.2)

}

if (FALSE) {
tissue_subtissue_heatmap <- fread(panel_input[["tissue_subtissue_heatmap"]])
z_matrix_subtissue       <- fread(panel_input[["z_matrix_subtissue"]])
stopifnot(all(c("Tissue", "Predicted_Tissue", "Pct") %in% names(tissue_subtissue_heatmap)))

hm <- copy(tissue_subtissue_heatmap)

hm_order  <- sort(unique(hm$Tissue))
hm_short  <- shorten_gtex(hm_order)
N_SUBTISSUES <- length(hm_order)
stopifnot(!anyDuplicated(hm_short))

hm[, row_lab := factor(shorten_gtex(Tissue), levels = rev(hm_short))]
hm[, col_lab := factor(shorten_gtex(Predicted_Tissue), levels = hm_short)]

SUBTISSUE_LABEL_MIN <- 20
hm[, is_diagonal := Tissue == Predicted_Tissue]
hm[, offdiag_rank := frank(fifelse(is_diagonal, Inf, -Pct), ties.method = "first"),
   by = Tissue]
hm[, show_label := !is_diagonal & Pct >= SUBTISSUE_LABEL_MIN & offdiag_rank == 1]
hm_labels <- hm[show_label == TRUE][order(-Pct)]

diag_cells <- hm[Tissue == Predicted_Tissue]
diag_cells <- merge(diag_cells, z_matrix_subtissue[, .(Tissue, tissue_correct)],
                    by.x = "Predicted_Tissue", by.y = "Tissue", all.x = TRUE)
stopifnot(!anyNA(diag_cells$tissue_correct))

plot_H_core <- ggplot(hm, aes(x = col_lab, y = row_lab, fill = Pct)) +
  geom_tile(colour = "#f0f0f0", linewidth = 0.1) +
  geom_tile(data = diag_cells[tissue_correct == TRUE],
            aes(x = col_lab, y = row_lab), fill = NA, colour = "black",
            linewidth = 0.35, inherit.aes = FALSE) +
  geom_tile(data = diag_cells[tissue_correct == FALSE],
            aes(x = col_lab, y = row_lab), fill = NA, colour = "red",
            linetype = "22", linewidth = 0.35, inherit.aes = FALSE) +
  geom_text(data = hm_labels,
            aes(label = sprintf("%.0f", Pct), colour = Pct >= 65),
            size = pt_mm(FS_DENSE_VALUE), check_overlap = TRUE, show.legend = FALSE) +
  scale_fill_gradientn(colours = GREEN_SCALE, limits = c(0, 100),
                       name = "% of pooled samples (column-normalised)",
                       guide = guide_colourbar(barwidth = unit(24, "mm"),
                                               barheight = unit(1.8, "mm"),
                                               title.position = "left",
                                               title.vjust = 1)) +
  scale_colour_manual(values = TILE_TEXT, guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Predicted GTEx subtissue", y = "GTEx subtissue") +
  theme_nature_heatmap(x_angle = 45, label_size = FS_HEAT_LABEL_H) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

panel_H <- add_tag(wrap_elements(full = plot_H_core), "H")

options(repr.plot.width = 7.2, repr.plot.height = 4.5)

}

FIG_W <- 183
FIG_H <- 285
stopifnot(FIG_W == 183, FIG_H == 285)

# C: pooled pseudobulk computational timing moved from Figure 2.
runtime_seed_totals <- fread(panel_input[["runtime_seed_totals"]])
stopifnot(all(c("method", "seed", "total_minutes", "n_datasets", "n_runs") %in%
              names(runtime_seed_totals)),
          nrow(runtime_seed_totals) == 30L,
          all(runtime_seed_totals$n_datasets == 6L),
          all(runtime_seed_totals$n_runs == 6L),
          all(is.finite(runtime_seed_totals$total_minutes)),
          all(runtime_seed_totals$total_minutes > 0))
timing_method_order <- runtime_seed_totals[,
  .(median_minutes = median(total_minutes)), by = method][
    order(median_minutes), as.character(method)]
runtime_seed_totals[, method := factor(method, levels = timing_method_order)]
panel_C_timing <- ggplot(runtime_seed_totals,
                         aes(method, total_minutes, fill = method)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, color = "black",
               linewidth = 0.22, alpha = 0.85) +
  geom_point(aes(shape = factor(seed)),
             position = position_jitter(width = 0.07, height = 0),
             size = 0.8, fill = "white", color = "black", stroke = 0.22) +
  scale_fill_manual(values = METHOD_COLORS, na.value = "grey70") +
  scale_shape_manual(values = c(`123` = 21, `456` = 22, `789` = 24),
                     guide = "none") +
  scale_y_log10(breaks = scales::log_breaks(n = 4),
                labels = scales::label_number(accuracy = 0.1)) +
  labs(x = NULL, y = "Total time") +
  theme_nature_methods(grid = "none") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1,
                                   size = FS_HEAT_LABEL_H),
        axis.text.y = element_text(size = FS_HEAT_LABEL_H),
        axis.title.y = element_text(size = FS_LEGEND_TITLE),
        plot.margin = margin(1, 1, 1, 1, "mm"))
panel_C_timing <- add_overlay_tag(panel_C_timing, "C")

# D: orthogonality of evaluation collections to the GO:BP training prior.
ortho <- fread(panel_input[["orthogonality_per_term"]])
ortho[, collection := factor(collection,
  levels = c("GTEx Tissues", "Azimuth 2023", "CellMarker Human",
             "Allen Brain Atlas"))]
ortho_means <- ortho[, .(mean_jaccard = mean(max_jaccard), n = .N), by = collection]
ortho_means[, label := sprintf("mean %.3f", mean_jaccard)]

panel_C_ortho <- ggplot(ortho, aes(max_jaccard)) +
  geom_density(fill = "#2166AC", colour = "#2166AC", alpha = 0.28,
               linewidth = 0.3, adjust = 1.2) +
  geom_vline(data = ortho_means, aes(xintercept = mean_jaccard),
             linetype = "dashed", linewidth = 0.3) +
  geom_text(data = ortho_means,
            aes(x = mean_jaccard, y = Inf, label = label),
            hjust = -0.06, vjust = 1.4, size = pt_mm(FS_LEGEND_TITLE),
            family = FONT_FAMILY) +
  facet_wrap(~collection, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 0.30, 0.05),
                     expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = "Maximum Jaccard similarity to any GO:BP training term",
       y = "Density") +
  theme_nature_methods(grid = "none") +
  theme(strip.text = element_text(size = FS_SUBTITLE, face = "bold"),
        axis.text = element_text(size = FS_LEGEND_TITLE),
        plot.margin = margin(1, 1, 1, 1, "mm"))
panel_C_ortho <- add_overlay_tag(panel_C_ortho, "D")

# F: former Figure 2 D marker-recovery grid.
marker_ready <- fread(panel_input[["marker"]])
stopifnot(unique(marker_ready$marker_source) == "CellMarker + Allen + Azimuth")
GRID_COLORS <- unlist(cfg$GRID_COLORS)
GRID_VALUES <- as.numeric(cfg$GRID_VALUES)

marker_lv_label <- function(x) {
  vapply(strsplit(x, " - ", fixed = TRUE), function(parts) {
    lv <- tail(parts, 1)
    gene_set <- parts[1]
    gene_set <- short_ct(gene_set)
    gene_set <- gsub("^Human ", "", gene_set)
    gene_set <- gsub("Astrocyte", "Astro.", gene_set, ignore.case = TRUE)
    gene_set <- gsub("Ciliated epithelial cell", "Ciliated ep.", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Inhibitory neuron", "Inhibitory", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Microglial cell", "Microglia", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Oligodendrocytes", "Oligodendro.", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Endothelial cell", "Endothelial", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Cardiomyocyte", "CM", gene_set, ignore.case = TRUE)
    gene_set <- gsub("Natural killer", "NK", gene_set, ignore.case = TRUE)
    gene_set <- gsub("Lymphatic endothelial", "Lymphatic EC", gene_set,
                     ignore.case = TRUE)
    gene_set <- gsub("Central memory CD4\\+ T cell", "Central mem. CD4", gene_set)
    gene_set <- gsub("Memory CD8\\+ T cell", "Memory CD8", gene_set)
    gene_set <- gsub("Natural killer T ?\\(NKT\\) cell", "NKT cell", gene_set)
    gene_set <- gsub("Type II pneumocyte", "Type II pneumo.", gene_set)
    if (identical(gene_set, lv)) lv else paste(gene_set, lv)
  }, character(1), USE.NAMES = FALSE)
}

make_marker_panel <- function(d) {
  ds <- d$dataset[1]
  d[, marker_row_label := short_ct(marker_row_label)]
  d[, lv_col_label := marker_lv_label(lv_col_label)]
  row_levels <- unique(d[order(marker_row_rank), marker_row_label])
  col_levels <- unique(d[order(lv_col_rank), lv_col_label])
  d[, marker_row_label := factor(marker_row_label, levels = row_levels)]
  d[, lv_col_label := factor(lv_col_label, levels = col_levels)]
  diagonal <- d[diag_recovered == TRUE]

  ggplot(d, aes(lv_col_label, marker_row_label)) +
    geom_point(aes(size = neg_log10_fdr, colour = row_effect)) +
    {if (nrow(diagonal)) geom_point(data = diagonal, aes(size = neg_log10_fdr),
                                    shape = 1, colour = "black", stroke = 0.5)} +
    scale_size_continuous(limits = c(0, d$Q_LIM[1]), range = c(0.25, 2.1),
                          name = expression("Combined ORA " * -log[10] ~ "(FDR)")) +
    scale_colour_gradientn(colours = GRID_COLORS, values = GRID_VALUES,
                           limits = c(-d$Z_LIM[1], d$Z_LIM[1]),
                           name = "LV effect") +
    scale_x_discrete(drop = FALSE, labels=function(x)vapply(x,function(t)paste(strwrap(t,15),collapse="\n"),character(1))) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian() +
    labs(x = NULL, y = NULL,
         title = paste0(DATASET_LABELS[ds], " (n = ", d$n_samples[1], ")")) +
    theme_nature_heatmap(x_angle = 90, label_size = FS_HEAT_LABEL_H) +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = GRID_W),
          plot.title = element_text(size = FS_SUBTITLE, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 5, lineheight = 0.88),
          axis.text.y = element_text(size = 5),
          plot.margin = margin(0.5, 0.7, 0.5, 0.7, "mm"))
}

marker_split <- split(copy(marker_ready), marker_ready$dataset)
marker_panels <- lapply(c(DATASETS_ROW1, DATASETS_ROW2),
                        function(ds) make_marker_panel(marker_split[[ds]]))
marker_legend <- get_legend(
  marker_panels[[1]] +
    guides(colour = guide_colorbar(order = 1, title.position = "top",
                                   barwidth = unit(15, "mm"),
                                   barheight = unit(1.5, "mm")),
           size = guide_legend(order = 2, title.position = "top", nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "horizontal", legend.text = element_text(size = 4.7),
          legend.title = element_text(size = 5))
)
marker_top <- plot_grid(plotlist = lapply(marker_panels[1:3],
                                          function(p) p + theme(legend.position = "none")),
                        nrow = 1)
marker_bottom <- plot_grid(plotlist = lapply(marker_panels[4:6],
                                             function(p) p + theme(legend.position = "none")),
                           nrow = 1)
marker_body <- plot_grid(marker_top, marker_bottom, ncol = 1,
                         rel_heights = c(1, 1))
panel_E_marker <- plot_grid(marker_body, marker_legend, ncol = 1,
                            rel_heights = c(1, 0.07))
panel_E_marker <- add_overlay_tag(panel_E_marker, "F")

