#!/usr/bin/env Rscript
source('R/common.R'); source('R/expression.R'); source('R/plot_style.R')
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize)
  library(RColorBrewer); library(cluster); library(ggplot2)})
cfg <- read_config()
out <- new_stage(cfg, 'expression_clusters')
x <- as.matrix(readRDS(cfg$inputs$expression)); x[is.na(x)] <- 0
map <- expression_gene_map(cfg)
symbols <- map$gene_name[match(colnames(x), map$gene_id)]
symbols[is.na(symbols)] <- colnames(x)[is.na(symbols)]
colnames(x) <- symbols
stopifnot(all(supercluster_order %in% rownames(x)))
gcep <- read_input(cfg, 'clingen', c('Gene', 'Classification'))
gcep[, Gene := trimws(gsub('HGNC:[0-9]+', '', Gene))]
selected <- unique(gcep[Classification %chin% c('Definitive', 'Strong', 'Moderate', 'Limited'), Gene])
genes <- intersect(unique(c(selected, 'CDK15')), colnames(x))
if (length(genes) < 5) stop('Too few selected genes found in the expression matrix')
filtered <- x[supercluster_order, genes, drop = FALSE]
correlation <- cor(filtered, method = 'pearson', use = 'pairwise.complete.obs')
correlation[is.na(correlation)] <- 0; diag(correlation) <- 1
distance <- as.dist(1-correlation)
tree <- hclust(distance, method = 'average')
k <- 2:min(10, ncol(filtered)-1)
scores <- data.table(k = k, silhouette = vapply(k, function(n) {
  mean(cluster::silhouette(cutree(tree, n), distance)[, 3])
}, numeric(1)))
write_tsv(scores, file.path(out, 'silhouette.tsv'))
p <- ggplot(scores, aes(k, silhouette)) + geom_line(linewidth = 1, colour = 'steelblue') +
  geom_point(size = 2, colour = 'steelblue') + scale_x_continuous(breaks = k) +
  labs(x = 'Number of clusters (k)', y = 'Average silhouette score') + theme_classic()
ggsave(file.path(out, 'cluster_silhouette.pdf'), p, width = 3.5, height = 3, device = cairo_pdf)

# Keep the heatmap's displayed order when naming Groups A-D.
# Raw cutree labels are not interchangeable with these displayed group labels.
ht <- Heatmap(filtered, name = 'Mean z-score',
  col = colorRamp2(seq(-2, 2, length.out = 100),
    rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))),
  border_gp = grid::gpar(col = 'black', lwd = 1.5),
  rect_gp = grid::gpar(col = '#c6c6c6', lwd = .5),
  cluster_rows = FALSE, clustering_method_columns = 'average',
  clustering_distance_columns = 'pearson', column_split = 4,
  column_gap = grid::unit(2, 'mm'), column_dend_height = grid::unit(2, 'cm'),
  column_title = 'Cell-type expression patterns of ALS-associated genes',
  column_title_gp = grid::gpar(fontsize = 13, fontfamily = cfg$font_family),
  top_annotation = HeatmapAnnotation(cluster = anno_block(
    gp = grid::gpar(fill = c('#DD5129FF', '#FAB255FF', '#43B284FF', '#0F7BA2FF'), col = NA),
    height = grid::unit(4, 'mm'), labels = LETTERS[1:4],
    labels_gp = grid::gpar(fontsize = 11, fontfamily = cfg$font_family))),
  row_order = supercluster_order, row_split = c(rep('Neuronal', 21), rep('Non neuronal', 10)),
  row_gap = grid::unit(1.5, 'mm'), row_names_side = 'right',
  row_title_gp = grid::gpar(fontsize = 11, fontfamily = cfg$font_family),
  row_names_gp = grid::gpar(fontsize = 10, fontfamily = cfg$font_family),
  column_names_gp = grid::gpar(fontsize = 10, fontfamily = cfg$font_family),
  width = grid::unit(ncol(filtered)*4, 'mm'), height = grid::unit(nrow(filtered)*4, 'mm'),
  heatmap_legend_param = list(direction = 'horizontal', title_position = 'topcenter',
    border = 'black', at = c(-2, 0, 2), legend_width = grid::unit(25, 'mm')))
cairo_pdf(file.path(out, 'expression_heatmap.pdf'), width = 12, height = 7)
drawn <- draw(ht)
dev.off()
groups <- lapply(column_order(drawn), function(i) colnames(filtered)[i])
names(groups) <- paste0('Group ', LETTERS[seq_along(groups)])
saveRDS(groups, file.path(out, 'gene_clusters.rds'))
write_tsv(rbindlist(lapply(names(groups), function(g) data.table(group = g, gene = groups[[g]]))),
          file.path(out, 'gene_clusters.tsv'))
record_session(out)
