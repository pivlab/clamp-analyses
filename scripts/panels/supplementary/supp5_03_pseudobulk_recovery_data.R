library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

result_dir <- here::here('output/03_model_biology/02_archs4/03_projections/04_pseudobulk_recovery_archs4')
projection <- read_csv(file.path(result_dir, 'pb_canon_dim.csv'), show_col_types = FALSE)
clamp_and_null <- read_csv(file.path(result_dir, 'pb_dimcontrol.csv'), show_col_types = FALSE)
full_null <- read_csv(here::here('output/03_model_biology/02_archs4/03_projections/pseudobulk_recovery/pseudobulk_recovery_long.csv'), show_col_types = FALSE) %>% filter(method == 'ARCHS4_null_projection')

per_dataset <- projection %>%
  select(dataset, cell_type, canon_k) %>%
  left_join(clamp_and_null %>% select(dataset, cell_type, clampfull_k, null_k), by = c('dataset', 'cell_type')) %>%
  group_by(dataset) %>%
  summarise(across(c(canon_k, clampfull_k, null_k), ~mean(.x, na.rm = TRUE)), .groups = 'drop')
per_dataset
per_dataset %>% summarise(across(c(canon_k, clampfull_k, null_k), mean))

method_colours <- c('CLAMPfull' = '#3e348b', 'ARCHS4 projection' = '#923155', 'ARCHS4 projection null' = '#bdbdbd')
comparisons <- tibble(a = c('CLAMPfull', 'ARCHS4 projection'), b = c('ARCHS4 projection', 'ARCHS4 projection null'))

format_p <- function(p) if (p < 0.001) format(p, scientific = TRUE, digits = 2) else sprintf('%.3f', p)

add_comparisons <- function(p, d, y_base, y_step, bracket_height) {
  for (i in seq_len(nrow(comparisons))) {
    a <- comparisons$a[i]; b <- comparisons$b[i]
    va <- d %>% filter(method == a) %>% arrange(replicate) %>% pull(value)
    vb <- d %>% filter(method == b) %>% arrange(replicate) %>% pull(value)
    p_value <- suppressWarnings(wilcox.test(va, vb, paired = TRUE)$p.value)
    x1 <- match(a, levels(d$method)); x2 <- match(b, levels(d$method)); y <- y_base + (i - 1) * y_step
    p <- p + annotate('segment', x = x1, xend = x1, y = y, yend = y + bracket_height, linewidth = 0.45) +
      annotate('segment', x = x1, xend = x2, y = y + bracket_height, yend = y + bracket_height, linewidth = 0.45) +
      annotate('segment', x = x2, xend = x2, y = y + bracket_height, yend = y, linewidth = 0.45) +
      annotate('text', x = (x1 + x2) / 2, y = y + bracket_height + 0.018, label = format_p(p_value), size = 4.8)
  }
  p
}

plot_recovery <- function(d, title) {
  means <- d %>% group_by(method) %>% summarise(mean_value = mean(value), max_value = max(value), .groups = 'drop') %>% arrange(desc(mean_value))
  d$method <- factor(d$method, levels = means$method)
  means$method <- factor(means$method, levels = means$method)
  y_base <- max(means$max_value) + 0.07
  y_max <- y_base + nrow(comparisons) * 0.12
  p <- ggplot(d, aes(method, value, fill = method)) +
    geom_boxplot(outlier.shape = NA, width = 0.52, colour = 'black', linewidth = 0.7, alpha = 0.82) +
    geom_point(position = position_jitter(width = 0.055, seed = 7), shape = 21, size = 2.25, fill = 'white', colour = '#4d4d4d', stroke = 0.45, alpha = 0.8) +
    geom_point(data = means, aes(x = method, y = mean_value), inherit.aes = FALSE, shape = 23, size = 5, fill = 'white', colour = 'black', stroke = 0.8) +
    geom_text(data = means, aes(x = method, y = max_value + 0.035, label = sprintf('%.3f', mean_value)), inherit.aes = FALSE, fontface = 'bold', size = 6) +
    scale_fill_manual(values = method_colours) +
    scale_y_continuous(breaks = seq(0, y_max, by = 0.25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y_max), clip = 'off') +
    labs(title = title, x = NULL, y = 'Max Pearson r per cell type') +
    theme_classic(base_size = 17) +
    theme(panel.grid.major.y = element_line(colour = '#D9D9D9', linewidth = 0.55), panel.grid.minor.y = element_line(colour = '#EEEEEE', linewidth = 0.35), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = 'none', plot.title = element_text(hjust = 0.5, face = 'plain'), plot.margin = margin(8, 12, 8, 8))
  add_comparisons(p, d, y_base, 0.12, 0.03)
}

k_matched <- bind_rows(
  transmute(clamp_and_null, method = 'CLAMPfull', value = clampfull_k, replicate = paste(dataset, cell_type)),
  transmute(projection, method = 'ARCHS4 projection', value = canon_k, replicate = paste(dataset, cell_type)),
  transmute(clamp_and_null, method = 'ARCHS4 projection null', value = null_k, replicate = paste(dataset, cell_type))
)
full_lvs <- bind_rows(
  transmute(clamp_and_null, method = 'CLAMPfull', value = clampfull_k, replicate = paste(dataset, cell_type)),
  transmute(projection, method = 'ARCHS4 projection', value = canon_1728, replicate = paste(dataset, cell_type)),
  transmute(full_null, method = 'ARCHS4 projection null', value = cor, replicate = paste(dataset, cell_type))
)


for (nm in c("k_matched","full_lvs")) data.table::fwrite(get(nm),file.path(out,paste0("03_pseudobulk_recovery_",nm,".csv")))
