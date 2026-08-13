#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (i in seq(1L, length(args), by = 2L)) {
    if (i == length(args)) stop("Missing value for ", args[[i]])
    out[[gsub("-", "_", substring(args[[i]], 3L), fixed = TRUE)]] <- args[[i + 1L]]
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", gsub("_", "-", name))
  value
}

args <- parse_cli()
ora_root <- required_arg(args, "ora_root")
coverage_out <- required_arg(args, "coverage_out")
cross_out <- required_arg(args, "cross_out")
panel_out <- required_arg(args, "panel_out")
figure_prefix <- required_arg(args, "figure_prefix")

# Only published ORA directories count.  run_ora.R builds a complete result set,
# `complete` marker included, inside `<out>.tmp.<pid>` and publishes it with an
# atomic rename, so an interrupted run leaves behind a directory that looks
# finished.  Matching the canonical <dataset>/rs<fraction>/seed<seed>/<model>/
# <database> layout skips those instead of double-counting them.
summary_files <- list.files(ora_root, pattern = "^summary\\.csv$", recursive = TRUE, full.names = TRUE)
canonical <- "/[^/]+/rs[0-9]+/seed[0-9]+/[^./]+/[^./]+/summary\\.csv$"
summary_files <- grep(canonical, summary_files, value = TRUE)
if (!length(summary_files)) stop("No published ORA summaries found below ", ora_root)
complete_files <- file.path(dirname(summary_files), "complete")
if (any(!file.exists(complete_files))) {
  stop("Incomplete ORA directories: ", paste(dirname(summary_files)[!file.exists(complete_files)], collapse = ", "))
}

coverage <- rbindlist(lapply(summary_files, fread), fill = TRUE)
required <- c(
  "dataset", "fraction", "seed", "model", "database", "database_label",
  "n_samples", "eligible_pathways", "recovered_pathways", "recovered_percent", "fdr"
)
if (!all(required %in% names(coverage))) {
  stop("ORA summaries are missing: ", paste(setdiff(required, names(coverage)), collapse = ", "))
}
if (anyDuplicated(coverage[, .(dataset, fraction, seed, model, database)])) {
  stop("Duplicate model/database ORA summaries found")
}
setorder(coverage, dataset, database, fraction, model, seed)

cross <- coverage[fraction == 100]
panel <- coverage[, .(
  mean_recovered = mean(recovered_pathways),
  sd_recovered = if (.N > 1L) sd(recovered_pathways) else NA_real_,
  mean_percent = mean(recovered_percent),
  eligible_pathways = unique(eligible_pathways),
  n_observations = .N,
  median_samples = as.integer(median(n_samples))
), by = .(dataset, fraction, model, database, database_label, fdr)]

dir.create(dirname(coverage_out), recursive = TRUE, showWarnings = FALSE)
fwrite(coverage, coverage_out)
fwrite(cross, cross_out)
fwrite(panel, panel_out)

model_labels <- c(CLAMPfull = "CLAMPfull (GO:BP prior)", CLAMPbase = "CLAMPbase")
model_colors <- c(CLAMPfull = "#0072B2", CLAMPbase = "#777777")
database_order <- c("reactome", "canonical", "cellmarker")
coverage[, database := factor(database, levels = database_order)]
coverage[, model := factor(model, levels = c("CLAMPfull", "CLAMPbase"))]

arch <- coverage[dataset == "archs4"]
fraction_order <- sort(unique(arch$fraction))
sample_labels <- arch[, .(n_samples = as.integer(median(n_samples))), by = fraction]
sample_labels[, label := sprintf("%s\n%s samples", paste0(fraction, "%"), format(n_samples, big.mark = ","))]
arch[, fraction_label := factor(
  fraction,
  levels = sample_labels$fraction,
  labels = sample_labels$label
)]
arch_full <- arch[model == "CLAMPfull"]
arch_base <- arch[model == "CLAMPbase"]

arch_means <- arch[, .(
  recovered_pathways = mean(recovered_pathways),
  recovered_percent = mean(recovered_percent)
), by = .(database, database_label, model, fraction_label)]
arch_full_labels <- arch_means[model == "CLAMPfull"]
arch_full_labels[, annotation := sprintf("%.0f\n%.0f%%", recovered_pathways, recovered_percent)]

p_arch <- ggplot(arch_full, aes(fraction_label, recovered_pathways, group = fraction_label)) +
  geom_boxplot(
    width = 0.48, outlier.shape = NA, alpha = 0.25, linewidth = 0.45,
    fill = model_colors[["CLAMPfull"]], colour = model_colors[["CLAMPfull"]]
  ) +
  geom_point(
    position = position_jitter(width = 0.06, height = 0),
    size = 1.5, alpha = 0.85, colour = model_colors[["CLAMPfull"]]
  ) +
  geom_line(
    data = arch_means[model == "CLAMPfull"],
    aes(group = 1),
    linetype = "dashed", linewidth = 0.45, colour = "#222222"
  ) +
  geom_point(
    data = arch_base,
    aes(group = NULL),
    shape = 95, size = 5, colour = model_colors[["CLAMPbase"]], alpha = 0.75
  ) +
  geom_text(
    data = arch_full_labels,
    aes(label = annotation, group = NULL),
    vjust = -0.7, colour = "black", size = 2.8, lineheight = 0.9
  ) +
  facet_wrap(~ database_label, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "Training compendium coverage", y = "Recovered pathways") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_line(colour = "#E3E3E3", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

cross_full <- cross[model == "CLAMPfull"]
cross_base <- cross[model == "CLAMPbase"]
dataset_order <- c("archs4", "gtex", "recount2")
dataset_names <- c(archs4 = "ARCHS4", gtex = "GTEx", recount2 = "recount2")
dataset_samples <- cross_full[, .(n_samples = as.integer(median(n_samples))), by = dataset]
dataset_samples[, label := sprintf("%s\n%s samples", dataset_names[dataset], format(n_samples, big.mark = ","))]
cross_full[, dataset_label := factor(
  dataset,
  levels = dataset_order,
  labels = dataset_samples[match(dataset_order, dataset), label]
)]
cross_base[, dataset_label := factor(
  dataset,
  levels = dataset_order,
  labels = dataset_samples[match(dataset_order, dataset), label]
)]
cross_means <- cross_full[, .(
  recovered_pathways = mean(recovered_pathways),
  recovered_percent = mean(recovered_percent)
), by = .(database, database_label, dataset_label)]
cross_means[, annotation := sprintf("%.0f\n%.0f%%", recovered_pathways, recovered_percent)]

p_cross <- ggplot(cross_full, aes(dataset_label, recovered_pathways, group = dataset_label)) +
  geom_boxplot(width = 0.48, fill = "#56B4E9", colour = "#0072B2", outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.06, height = 0, size = 1.7, colour = "#0072B2") +
  geom_point(
    data = cross_base,
    shape = 95, size = 5, colour = model_colors[["CLAMPbase"]], alpha = 0.75
  ) +
  geom_line(
    data = cross_means,
    aes(group = 1),
    linetype = "dashed", linewidth = 0.45, colour = "#333333"
  ) +
  geom_text(
    data = cross_means,
    aes(label = annotation),
    vjust = -0.7, colour = "black", size = 2.8, lineheight = 0.9
  ) +
  facet_wrap(~ database_label, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = NULL, y = "Recovered pathways") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_line(colour = "#E3E3E3", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(figure_prefix, "_archs4.png"), p_arch, width = 11, height = 5.8, dpi = 300)
ggsave(paste0(figure_prefix, "_archs4.pdf"), p_arch, width = 11, height = 5.8)
ggsave(paste0(figure_prefix, "_cross_dataset.png"), p_cross, width = 9, height = 5.5, dpi = 300)
ggsave(paste0(figure_prefix, "_cross_dataset.pdf"), p_cross, width = 9, height = 5.5)

print(p_arch)
print(p_cross)
