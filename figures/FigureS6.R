#!/usr/bin/env Rscript
# Plot cached silhouette scores without repeating the clustering analysis.
source('R/common.R')
suppressPackageStartupMessages(library(ggplot2))
cfg <- read_config()
path <- cfg$inputs$silhouette
if (is.null(path)) path <- file.path(cfg$output_dir, 'expression_clusters', 'silhouette.tsv')
scores <- fread(path)
require_columns(scores, c('k', 'silhouette'), 'silhouette scores')
stopifnot(!anyDuplicated(scores$k), all(is.finite(scores$silhouette)))
setorder(scores, k)
out <- new_stage(cfg, 'FigureS6')
p <- ggplot(scores, aes(k, silhouette)) + geom_line(linewidth = 1, colour = 'steelblue') +
  geom_point(size = 2, colour = 'steelblue') + scale_x_continuous(breaks = scores$k) +
  labs(x = 'Number of clusters (k)', y = 'Average silhouette score') + theme_classic()
ggsave(file.path(out, 'FigureS6.pdf'), p, width = 3.5, height = 3, device = cairo_pdf)

record_session(out)
