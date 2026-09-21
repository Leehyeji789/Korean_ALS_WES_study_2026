#!/usr/bin/env Rscript
source('R/common.R'); source('R/expression.R')
cfg <- read_config()
out <- new_stage(cfg, 'cell_type_enrichment')
deg <- read_input(cfg, 'differential_expression',
  c('names', 'group', 'pvals_adj', 'logfoldchanges', 'pct_nz_group'))
map <- expression_gene_map(cfg)
background <- unique(map[gene_id %chin% unique(deg$names) & gene_type == 'protein_coding', gene_name])
symbols <- map$gene_name[match(deg$names, map$gene_id)]
deg[!is.na(symbols), names := symbols[!is.na(symbols)]]
groups <- readRDS(file.path(cfg$output_dir, 'expression_clusters/gene_clusters.rds'))
missing <- setdiff(unlist(groups, use.names = FALSE), background)
if (length(missing)) stop('Cluster genes absent from DEG background: ', paste(missing, collapse = ', '))
sig <- deg[pvals_adj < .05 & logfoldchanges > .3 & pct_nz_group > .25]
rows <- list()
for (cell in unique(sig$group)) {
  cell_genes <- intersect(sig[group == cell, unique(names)], background)
  for (g in names(groups)) {
    group_genes <- intersect(groups[[g]], background)
    overlap <- intersect(cell_genes, group_genes)
    a <- length(overlap); b <- length(cell_genes)-a; c <- length(group_genes)-a
    d <- length(background)-a-b-c
    stopifnot(min(a,b,c,d) >= 0)
    ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2, byrow = TRUE), alternative = 'greater')
    rows[[paste(cell,g)]] <- data.table(Cluster = cell, Group = g, P.value = ft$p.value,
      OddsRatio = unname(ft$estimate), overlap_risk_deg = a, deg_only = b,
      risk_only = c, bg_only = d, overlapped_genes = paste(overlap, collapse = ','))
  }
}
d <- rbindlist(rows)
d[, padj := p.adjust(P.value, 'BH')]
write_tsv(d, file.path(out, 'enrichment.tsv'))
write_tsv(data.table(gene = background), file.path(out, 'tested_background.tsv'))
record_session(out)
