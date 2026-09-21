#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R'); source('R/expression.R')
suppressPackageStartupMessages({library(dplyr);library(cowplot);library(enrichplot)})
source('R/expression_style.R')
cfg <- read_config()
out <- new_stage(cfg,'FigureS7')
cluster_res <- readRDS(file.path(cfg$output_dir,'expression_clusters/gene_clusters.rds'))
cluster_order <- supercluster_order
fet_results_df <- as.data.frame(fread(file.path(cfg$output_dir,'cell_type_enrichment/enrichment.tsv')))
group_sizes <- sapply(cluster_res, length)
group_labels <- setNames(paste0(names(group_sizes), '\n(n=', group_sizes, ')'), names(group_sizes))

plot_df <- fet_results_df %>%
  mutate(
    log2OR = log2(OddsRatio),
    signif = case_when(
      padj < 0.05 & log2OR > 1 ~ 'padj < 0.05 & log2OR > 1',
      padj < 0.05              ~ 'padj < 0.05',
      P.value < 0.05           ~ 'p < 0.05',
      TRUE                     ~ 'p > 0.05'
    ),
    capped_log2OR = case_when(
      log2OR == Inf  ~ 4,
      log2OR == -Inf ~ -4,
      log2OR > 4     ~ 4,
      log2OR < -4    ~ -4,
      TRUE           ~ log2OR
    ),
    signif = factor(signif, levels = c(
      'padj < 0.05 & log2OR > 1',
      'padj < 0.05',
      'p < 0.05',
      'p > 0.05'
    )),
    Group = factor(Group, levels = names(cluster_res))
  )

p <- ggplot(plot_df, aes(x = factor(Cluster, levels = cluster_order), y = Group)) +

  geom_tile(fill = 'white', linewidth = 0.25, color = '#949494', show.legend = FALSE) +
  coord_fixed(ratio = 1) +

  geom_point(aes(fill = capped_log2OR, size = signif,
                 stroke = ifelse(signif == 'padj < 0.05 & log2OR > 1', 0.9, 0.5)),
             shape = 22, color = 'black') +

  geom_text(data = plot_df %>% filter(padj < 0.05),
            label = '*', size = 5, color = 'black', fontface = 'bold', vjust = .8) +

  scale_fill_gradient2(
    low = '#2166ac', mid = 'white', high = '#b2182b',
    limits = c(-log2(16), log2(16)),
    name = expression(log[2](OR)),
    guide = guide_colorbar(
      theme = theme(
        legend.frame = element_rect(color = 'black', linewidth = .5),
        legend.ticks = element_blank(),
        legend.key.height = unit(0.8, 'lines')
      )
    )
  ) +

  scale_size_manual(
    values = c('padj < 0.05 & log2OR > 1' = 8, 'padj < 0.05' = 6, 'p < 0.05' = 4, 'p > 0.05' = 2),
    labels = c(
      expression('adj.p. < 0.05 & ' * log[2] * '(OR) > 1'),
      'adj.p. < 0.05',
      'p. < 0.05',
      'p. > 0.05'
    ),
    guide = guide_legend(title = NULL)
  ) +

  scale_y_discrete(limits = rev(names(cluster_res)), labels = group_labels) +

  theme_minimal(base_size = 10) +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 9, color = 'black'),
    legend.direction = 'horizontal', legend.position = 'top',
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, colour = 'black'),
    axis.text.y = element_text(size = 10, colour = 'black')
  )
p

save_figure(p,file.path(out,'FigureS7_before_artwork_edit'),10,4.5,raster=FALSE)
source('R/manuscript_artwork.R')
preserve_artwork(cfg, 'figureS7_artwork', out, 'FigureS7')
record_session(out)
