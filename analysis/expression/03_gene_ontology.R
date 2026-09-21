#!/usr/bin/env Rscript
source('R/common.R')
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
cfg <- read_config()
out <- new_stage(cfg, 'gene_ontology')
groups <- readRDS(file.path(cfg$output_dir, 'expression_clusters/gene_clusters.rds'))
# These are the original GO settings. The background is the package's default
# annotated-gene universe, not the tested-DEG universe used for cell enrichment.
mf <- compareCluster(geneClusters = groups, OrgDb = 'org.Hs.eg.db', fun = 'enrichGO',
  ont = 'MF', pAdjustMethod = 'BH', keyType = 'SYMBOL', qvalueCutoff = .05,
  minGSSize = 10, maxGSSize = 500, readable = TRUE)
each <- lapply(groups, function(genes) enrichGO(gene = genes, OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', ont = 'ALL', pAdjustMethod = 'BH', qvalueCutoff = .05,
  minGSSize = 5, maxGSSize = 500, readable = TRUE))
write_tsv(as.data.table(as.data.frame(mf)), file.path(out, 'molecular_function.tsv'))
rows <- lapply(names(each), function(g) {
  x <- as.data.table(as.data.frame(each[[g]])); x[, Group := g]; x
})
write_tsv(rbindlist(rows, fill = TRUE), file.path(out, 'all_ontologies.tsv'))
saveRDS(list(comparative_mf = mf, per_group = each), file.path(out, 'enrichment_objects.rds'))
record_session(out)
