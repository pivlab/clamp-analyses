# Original analysis and provenance guards retained; layout/export lives in supp2.R.
panel_input <- lapply(jsonlite::fromJSON(here::here("scripts/panels/supplementary/supp2_inputs.json")),here::here)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(ggrastr)
  library(ggrepel)
  library(grid)
  library(svglite)
  library(ragg)
})


FIG_W <- 183
FIG_H <- 445
FONT_FAMILY <- "Helvetica"
FS_TAG <- 7
FS_TITLE <- 6.5
FS_AXIS_TITLE <- 6.5
FS_AXIS_TEXT <- 5
LINE_W <- 0.25

DATASET_LEVELS <- c("Heart_Datar2026", "PBMC_1k1k",
  "Brain_Mathys2023", "Brain_Xiong2023",
  "Lung_Sikkema2023", "PBMC_Perez2022")
DATASET_LABELS <- c(Heart_Datar2026="Heart: Datar", PBMC_1k1k="PBMC: 1k1k",
  Brain_Mathys2023="Brain: Mathys", Brain_Xiong2023="Brain: Xiong",
  Lung_Sikkema2023="Lung: Sikkema", PBMC_Perez2022="PBMC: Perez")
DATASET_COLORS <- c(Brain_Mathys2023="#5B8FF9", Brain_Xiong2023="#9270CA",
  Heart_Datar2026="#E8684A", PBMC_1k1k="#5AD8A6",
  PBMC_Perez2022="#F6BD16", Lung_Sikkema2023="#6DC8EC")
KNOWN_CELL_COLORS <- c(B_cell="#4477AA", Myeloid="#EE6677",
  NK="#228833", T_cell="#CCBB44", CD4_T="#66CCEE",
  CD8_T="#AA3377", CD14_Mono="#EE7733", CD16_Mono="#0077BB",
  DC="#BBBBBB", Plasma_B="#EE3377", gd_T="#009988")

purity <- fread(panel_input[["purity"]])
cells <- fread(panel_input[["umap_cells"]])
lvs <- fread(panel_input[["umap_lvs"]])
donor_bulk_recovery <- fread(panel_input[["donor_bulk_recovery"]])
donor_bulk_overall <- fread(panel_input[["donor_bulk_overall"]])
stopifnot(grepl("/06_donor_bulk_recovery/", panel_input[["purity"]], fixed=TRUE),
          grepl("/06_donor_bulk_recovery/", panel_input[["umap_lvs"]], fixed=TRUE),
          grepl("/06_donor_bulk_recovery/", panel_input[["donor_bulk_recovery"]], fixed=TRUE),
          grepl("/06_donor_bulk_recovery/", panel_input[["donor_bulk_overall"]], fixed=TRUE),
          nrow(donor_bulk_recovery) == 47L)

# Hard provenance guard for Supplementary Fig. 2a. These are donor-bulk
# projection results, distinct from the pseudobulk values in Supp. Fig. 1e.
assert_close <- function(actual, expected, tolerance, label) {
  if (length(actual) != 1L || !is.finite(actual) || abs(actual - expected) > tolerance)
    stop(sprintf("%s provenance check failed: observed %.6f, expected %.6f +/- %.6f",
                 label, actual, expected, tolerance))
}
assert_close(donor_bulk_overall$global_recovery_pct_pooled, 70.618975, 0.05, "donor-bulk pooled purity (%)")
assert_close(donor_bulk_overall$mean_recovery_pct_by_dataset, 77.455512, 0.05, "donor-bulk mean dataset purity (%)")
assert_close(donor_bulk_overall$mean_purity_ratio, 0.781929, 0.001, "donor-bulk mean purity ratio")
assert_close(donor_bulk_overall$mean_purity_lift, 19.595991, 0.05, "donor-bulk mean purity lift")
assert_close(sum(donor_bulk_recovery$diag_count) / sum(donor_bulk_recovery$n_top) * 100,
             donor_bulk_overall$global_recovery_pct_pooled, 1e-6, "donor-bulk pooled purity recomputation")
perez_b <- donor_bulk_recovery[dataset == "PBMC_Perez2022" & cell_type == "B_cell"]
stopifnot(nrow(perez_b) == 1L, perez_b$assigned_lv == "LV32",
          setequal(paste(purity$dataset, purity$cell_type),
                   paste(donor_bulk_recovery$dataset, donor_bulk_recovery$cell_type)))
activity_columns <- grep("^activity_[0-9]+$", names(cells), value=TRUE)
activity_columns <- activity_columns[order(as.integer(sub("activity_", "", activity_columns)))]
stopifnot(setequal(unique(cells$dataset), DATASET_LEVELS),
          setequal(unique(lvs$dataset), DATASET_LEVELS),
          length(activity_columns) == max(lvs$plot_order),
          lvs[order(plot_order), all(plot_order == seq_len(.N)), by=dataset]$V1,
          cells[, all(uniqueN(cell_index) == .N), by=dataset]$V1,
          all(is.finite(as.matrix(cells[, .(umap1, umap2)]))))
for (dataset_id in DATASET_LEVELS) {
  selected_columns <- paste0("activity_", seq_len(lvs[dataset == dataset_id, .N]))
  stopifnot(all(is.finite(as.matrix(cells[dataset == dataset_id,
    selected_columns, with=FALSE]))))
}

all_cell_types <- sort(unique(cells$mapped_cell_type))
CELL_COLORS <- setNames(grDevices::hcl.colors(length(all_cell_types), "Dark 3"),
                        all_cell_types)
known <- intersect(names(KNOWN_CELL_COLORS), names(CELL_COLORS))
CELL_COLORS[known] <- KNOWN_CELL_COLORS[known]
pretty_cell_type <- function(x) gsub("_", " ", x, fixed=TRUE)

theme_nm <- function() {
  theme_classic(base_size=FS_AXIS_TEXT, base_family=FONT_FAMILY) %+replace%
    theme(axis.line=element_line(linewidth=LINE_W, colour="black"),
      axis.ticks=element_line(linewidth=LINE_W, colour="black"),
      axis.text=element_text(size=FS_AXIS_TEXT, colour="black"),
      axis.title=element_text(size=FS_AXIS_TITLE, colour="black"),
      panel.grid=element_blank(), plot.margin=margin(1,1,1,1,"mm"))
}
add_tag <- function(p, label) {
  ggdraw() + draw_plot(p) +
    draw_label(label, x=0, y=1, hjust=0, vjust=1, size=FS_TAG,
               fontface="bold", fontfamily=FONT_FAMILY)
}

# A: top-1% annotated-cell purity across all matched LVs.
purity[, dataset := factor(dataset, levels=DATASET_LEVELS)]
purity[, dataset_label := factor(DATASET_LABELS[as.character(dataset)],
                                 levels=DATASET_LABELS[DATASET_LEVELS])]
panel_A <- ggplot(purity, aes(dataset_label, recovery_pct, fill=dataset)) +
  geom_boxplot(width=0.58, outlier.shape=NA, linewidth=0.28, alpha=0.92) +
  geom_jitter(width=0.10, size=0.70, shape=21, fill="#333333",
              colour="#333333", stroke=0, alpha=0.75) +
  scale_fill_manual(values=DATASET_COLORS, drop=FALSE) +
  scale_y_continuous(limits=c(0, 100), breaks=seq(0,100,25),
                     expand=expansion(mult=c(0.01,0.04))) +
  labs(x=NULL, y="Top-1% annotated-cell purity (%)") +
  theme_nm() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=25, hjust=1, size=5.5))

map_theme <- theme_void(base_family=FONT_FAMILY, base_size=4.2) +
  theme(panel.border=element_rect(colour="black", fill=NA, linewidth=0.25),
        plot.title=element_text(size=4.0, hjust=0.5, lineheight=0.92,
                                margin=margin(0,0,0.35,0,"mm")),
        plot.margin=margin(0.35,0.35,0.35,0.35,"mm"))

square_limits <- function(d, pad=0.04) {
  xr <- range(d$umap1, finite=TRUE)
  yr <- range(d$umap2, finite=TRUE)
  span <- max(diff(xr), diff(yr)) * (1 + 2 * pad)
  list(x=mean(xr) + c(-0.5, 0.5) * span,
       y=mean(yr) + c(-0.5, 0.5) * span)
}

feature_sheet <- function(dataset_id) {
  d <- cells[dataset == dataset_id]
  meta <- lvs[dataset == dataset_id][order(plot_order)]
  stopifnot(nrow(meta) > 0L, nrow(d) > 0L)
  limits <- square_limits(d)
  centroids <- d[, .(umap1=median(umap1), umap2=median(umap2)),
                 by=mapped_cell_type]
  centroids[, label := pretty_cell_type(mapped_cell_type)]
  annotation <- ggplot(d, aes(umap1, umap2, colour=mapped_cell_type)) +
    ggrastr::rasterise(geom_point(size=0.16, alpha=0.78), dpi=300) +
    ggrepel::geom_text_repel(data=centroids, aes(label=label),
      colour="black", size=1.20, fontface="bold", family=FONT_FAMILY,
      seed=123, box.padding=0.10, point.padding=0.05, min.segment.length=0,
      segment.size=0.16, max.overlaps=Inf, show.legend=FALSE) +
    scale_colour_manual(values=CELL_COLORS) +
    coord_equal(xlim=limits$x, ylim=limits$y, expand=FALSE) +
    labs(title="Cell-type annotation") + map_theme +
    theme(legend.position="none",
          plot.title=element_text(size=4.5, face="bold"))
  features <- lapply(seq_len(nrow(meta)), function(i) {
    m <- meta[plot_order == i]
    high <- CELL_COLORS[[m$cell_type]]
    score_col <- paste0("activity_", i)
    ggplot(d, aes(umap1, umap2, colour=.data[[score_col]])) +
      ggrastr::rasterise(geom_point(size=0.16), dpi=300) +
      scale_colour_gradientn(colours=c("#D8E2EF", "white", high)) +
      coord_equal(xlim=limits$x, ylim=limits$y, expand=FALSE) +
      labs(title=sprintf("%s · %s", pretty_cell_type(m$cell_type), m$LV)) +
      map_theme + theme(legend.position="none")
  })
  feature_grid <- plot_grid(plotlist=features, ncol=5, align="hv")
  map_row <- plot_grid(annotation, feature_grid, nrow=1,
                       rel_widths=c(0.30, 0.70))
  header <- ggdraw() + draw_label(DATASET_LABELS[[dataset_id]],
    x=0.5, y=0.45, hjust=0.5, size=FS_TITLE, fontface="bold",
    fontfamily=FONT_FAMILY)
  plot_grid(header, map_row, ncol=1, rel_heights=c(0.08, 0.92))
}

