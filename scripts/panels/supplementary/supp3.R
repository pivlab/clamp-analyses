source(here::here("scripts/panels/supplementary/style.R"))
set.seed(123)
source(here("scripts/panels/supplementary/supp3_panels.R"))
nm_sources(3,unlist(panel_input))
# Preserve all numeric labels selected by the source report; give the matrix a full page.
for(i in seq_along(panel_C$layers)) if(inherits(panel_C$layers[[i]]$geom,"GeomText")) {
 panel_C$layers[[i]]$geom_params$check_overlap <- FALSE
panel_C$layers[[i]]$aes_params$size <- 5/.pt
}
# The complete GTEx timing report records both elapsed time and peak resident
# memory for every seed. They are shown separately to preserve their scales.
runtime <- fread(panel_input[["gtex_runtime_per_fit"]])
runtime[, method := factor(method, levels = runtime[, .(median_minutes = median(elapsed_minutes)), by = method][order(median_minutes), method])]
runtime_theme <- nm_theme() + theme(legend.position = "none", axis.text.x = element_text(angle = 42, hjust = 1, size = 5))
panel_D <- ggplot(runtime, aes(method, elapsed_minutes, fill = method)) +
  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .25, alpha = .75) +
  geom_point(position = position_jitter(width = .08, height = 0), size = .65) +
  scale_fill_manual(values = MODEL_COLORS_RAW, na.value = "grey70") +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  labs(x = NULL, y = "Elapsed time (min)") + runtime_theme
panel_E <- ggplot(runtime, aes(method, peak_rss_gb, fill = method)) +
  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .25, alpha = .75) +
  geom_point(position = position_jitter(width = .08, height = 0), size = .65) +
  scale_fill_manual(values = MODEL_COLORS_RAW, na.value = "grey70") +
  labs(x = NULL, y = "Peak RAM (GB)") + runtime_theme
top <- plot_grid(
  nm_tag(panel_A,"A","Gene-subsampling robustness"),
  nm_tag(panel_B,"B","Tissue-clustering benchmark"),
  nm_tag(panel_D,"C","GTEx computational time"),
  nm_tag(panel_E,"D","GTEx peak memory"), ncol=2)
page <- plot_grid(top,
  nm_tag(panel_C,"E","Subtissue concordance with loading-matrix confirmation"),
  ncol=1,rel_heights=c(82,159))
nm_export(list(page),3,"GTEx robustness, runtime, memory and subtissue concordance")
