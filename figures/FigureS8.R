#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R'); source('R/expression.R')
suppressPackageStartupMessages({library(dplyr);library(cowplot);library(enrichplot)})
source('R/expression_style.R')
cfg <- read_config()
out <- new_stage(cfg,'FigureS8')
objects <- readRDS(file.path(cfg$output_dir,'gene_ontology/enrichment_objects.rds'))
ego_compare <- objects$comparative_mf
go_res <- objects$per_group
compare_go_table <- as.data.frame(ego_compare)

if (nrow(compare_go_table) == 0) {
  stop('No significant terms were identified by compareCluster.')
}

p_compare <- enrichplot::dotplot(
  ego_compare,
  showCategory = 5, includeAll = TRUE,
  color = 'p.adjust', size = 'Count',
  font.size = 10, label_format = 40,
  title = 'GO molecular function enrichment across gene groups'
) +
  scale_color_gradient(
    low = '#DD5129FF',
    high = '#0F7BA2FF',
    name = 'p.adjust'
  ) +
  scale_size_continuous(range = c(2.5, 7), name = 'Gene count') +
  theme_als_v1 +
  theme(
    plot.title.position = 'plot',
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.spacing.y = grid::unit(1, 'lines'),
    legend.position = 'right'
  )

p_compare

go_each_table <- bind_rows(lapply(names(go_res), function(group_name) {
  as.data.frame(go_res[[group_name]]) %>%
    mutate(Group = group_name, .before = 1)
})) %>%
  filter(ONTOLOGY %in% c('BP', 'MF'))

ratio_to_numeric <- function(x) {
  numerator <- as.numeric(sub('/.*', '', x))
  denominator <- as.numeric(sub('.*/', '', x))
  numerator / denominator
}

make_group_go_plot <- function(group_name, n_terms = 5, label_width = 35) {
  plot_df <- go_each_table %>%
    filter(Group == group_name) %>%
    group_by(ONTOLOGY) %>%
    slice_min(order_by = qvalue, n = n_terms, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(GeneRatio_numeric = ratio_to_numeric(GeneRatio)) %>%
    arrange(ONTOLOGY, GeneRatio_numeric, qvalue) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  ggplot(
    plot_df,
    aes(
      x = GeneRatio_numeric, y = Description,
      size = Count, color = p.adjust
    )
  ) +
    geom_point() +
    scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = label_width)
    ) +
    scale_color_gradient(
      low = '#DD5129FF', high = '#0F7BA2FF',
      limits = c(0, 0.05), oob = scales::squish,
      name = 'p.adjust'
    ) +
    scale_size_continuous(range = c(2.5, 7), name = 'Gene count') +
    facet_grid(
      ONTOLOGY ~ .,
      scales = 'free_y', space = 'free_y',
      labeller = as_labeller(c(BP = 'BP', MF = 'MF'))
    ) +
    labs(
      title = paste0('GO enrichment of ', group_name, ' genes'),
      x = 'Gene ratio', y = NULL
    ) +
    theme_als_v1 +
    theme(
      # plot.title.position = 'plot',
      axis.text.y = element_text(size = 9),
      panel.spacing.y = grid::unit(1, 'lines'),
      legend.position = 'right'
    )
}

p_group_a <- make_group_go_plot('Group A')
p_group_a
p_group_b <- make_group_go_plot('Group B')
p_group_c <- make_group_go_plot('Group C')
p_group_d <- make_group_go_plot('Group D')

###### 03. Assemble Figure S8 ######
p_groups_ab <- plot_grid(
  p_group_a, p_group_b, nrow = 1, rel_widths = c(1, 1.2),
  labels = c('B', 'C'), label_size = 16, label_fontfamily = 'Arial', vjust = 1.3
)
p_groups_ab

p_groups_cd <- plot_grid(
  p_group_c, p_group_d, nrow = 1, rel_widths = c(1, 1),
  labels = c('D', 'E'), label_size = 16, label_fontfamily = 'Arial', vjust = 1.3
)
p_groups_cd

p_groups <- plot_grid(
  p_groups_ab, p_groups_cd, ncol = 1, rel_heights = c(1, .9)
)
p_groups

p <- plot_grid(
  p_compare, p_groups, ncol = 2,
  rel_widths = c(1, 1.8),
  labels = c('A', ''),
  label_size = 16, label_fontfamily = 'Arial', vjust = 1.3
)
p

save_figure(p,file.path(out,'FigureS8_before_artwork_edit'),18,7.8,raster=FALSE)
source('R/manuscript_artwork.R')
preserve_artwork(cfg, 'figureS8_artwork', out, 'FigureS8')
record_session(out)
