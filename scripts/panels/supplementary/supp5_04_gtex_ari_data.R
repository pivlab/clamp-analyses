library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

result_dir <- here::here('output/03_model_biology/02_archs4/03_projections/05_gtex_ari_archs4')
projection <- read_csv(file.path(result_dir, 'gtex_canon_dim.csv'), show_col_types = FALSE)
null <- read_csv(file.path(result_dir, 'dim_control.csv'), show_col_types = FALSE)
local_ari <- read_csv(here::here('output/03_model_biology/01_gtex/00_kmeans_clustering/ari_data.csv'), show_col_types = FALSE)

tibble(
  method = c('ARCHS4 projection', 'Shared null'),
  mean_ARI_K578 = c(mean(projection$canon_578), mean(null$null_578))
)

method_colours <- c('CLAMPfull' = '#3e348b', 'ARCHS4 projection' = '#923155', 'RNA-Seq' = '#9c9f36', 'ARCHS4 projection null' = '#999999')
comparisons <- tibble(a = c('CLAMPfull', 'CLAMPfull', 'ARCHS4 projection'), b = c('ARCHS4 projection', 'RNA-Seq', 'ARCHS4 projection null'))

format_p <- function(p) if (p < 0.001) format(p, scientific = TRUE, digits = 2) else sprintf('%.3f', p)

add_comparisons <- function(p, d, y_base, y_step, bracket_height) {
  for (i in seq_len(nrow(comparisons))) {
    a <- comparisons$a[i]; b <- comparisons$b[i]
    va <- d %>% filter(method == a) %>% arrange(replicate) %>% pull(ARI)
    vb <- d %>% filter(method == b) %>% arrange(replicate) %>% pull(ARI)
    p_value <- suppressWarnings(wilcox.test(va, vb, paired = TRUE)$p.value)
    x1 <- match(a, levels(d$method)); x2 <- match(b, levels(d$method)); y <- y_base + (i - 1) * y_step
    p <- p + annotate('segment', x = x1, xend = x1, y = y, yend = y + bracket_height, linewidth = 0.45) +
      annotate('segment', x = x1, xend = x2, y = y + bracket_height, yend = y + bracket_height, linewidth = 0.45) +
      annotate('segment', x = x2, xend = x2, y = y + bracket_height, yend = y, linewidth = 0.45) +
      annotate('text', x = (x1 + x2) / 2, y = y + bracket_height + 0.011, label = format_p(p_value), size = 4.8)
  }
  p
}

plot_recovery <- function(d, title) {
  means <- d %>% group_by(method) %>% summarise(mean_ari = mean(ARI), max_ari = max(ARI), .groups = 'drop') %>% arrange(desc(mean_ari))
  d$method <- factor(d$method, levels = means$method)
  means$method <- factor(means$method, levels = means$method)
  y_base <- max(means$max_ari) + 0.025
  y_max <- y_base + nrow(comparisons) * 0.05
  p <- ggplot(d, aes(method, ARI, fill = method)) +
    geom_boxplot(outlier.shape = NA, width = 0.52, colour = 'black', linewidth = 0.7, alpha = 0.82) +
    geom_point(position = position_jitter(width = 0.055, seed = 7), shape = 21, size = 2.25, fill = 'white', colour = '#4d4d4d', stroke = 0.45, alpha = 0.8) +
    geom_point(data = means, aes(x = method, y = mean_ari), inherit.aes = FALSE, shape = 23, size = 5, fill = 'white', colour = 'black', stroke = 0.8) +
    geom_text(data = means, aes(x = method, y = max_ari + 0.010, label = sprintf('%.3f', mean_ari)), inherit.aes = FALSE, fontface = 'bold', size = 6) +
    scale_fill_manual(values = method_colours) +
    scale_y_continuous(breaks = seq(0, y_max, by = 0.2), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y_max), clip = 'off') +
    labs(title = title, x = NULL, y = 'Adjusted Rand index (ARI)') +
    theme_classic(base_size = 17) +
    theme(panel.grid.major.y = element_line(colour = '#D9D9D9', linewidth = 0.55), panel.grid.minor.y = element_line(colour = '#EEEEEE', linewidth = 0.35), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = 'none', plot.title = element_text(hjust = 0.5, face = 'plain'), plot.margin = margin(8, 12, 8, 8))
  add_comparisons(p, d, y_base, 0.05, 0.012)
}

k_matched <- bind_rows(
  transmute(filter(local_ari, method == 'CLAMPfull'), method = 'CLAMPfull', ARI = ari, replicate = row_number()),
  transmute(projection, method = 'ARCHS4 projection', ARI = canon_578, replicate = row_number()),
  transmute(filter(local_ari, method == 'RNA-Seq'), method = 'RNA-Seq', ARI = ari, replicate = row_number()),
  transmute(null, method = 'ARCHS4 projection null', ARI = null_578, replicate = row_number())
)
full_lvs <- bind_rows(
  transmute(filter(local_ari, method == 'CLAMPfull'), method = 'CLAMPfull', ARI = ari, replicate = row_number()),
  transmute(projection, method = 'ARCHS4 projection', ARI = canon_1728, replicate = row_number()),
  transmute(filter(local_ari, method == 'RNA-Seq'), method = 'RNA-Seq', ARI = ari, replicate = row_number()),
  transmute(null, method = 'ARCHS4 projection null', ARI = null_1728, replicate = row_number())
)


for (nm in c("k_matched","full_lvs")) data.table::fwrite(get(nm),file.path(out,paste0("04_gtex_ari_",nm,".csv")))
