#!/usr/bin/env Rscript
source('R/common.R')
suppressPackageStartupMessages(library(DiscreteFDR))
cfg <- read_config()
out <- new_stage(cfg, 'burden')
genes <- read_input(cfg, 'genes', c('gene_id', 'gene_symbol', 'chromosome', 'start', 'end'))
stopifnot(!anyDuplicated(genes$gene_id))
results <- list()

for (label in names(cfg$af_thresholds)) {
  a <- fread(file.path(cfg$output_dir, 'counts', paste0('cases_', label, '.tsv')))
  g <- fread(file.path(cfg$output_dir, 'counts', paste0('controls_', label, '.tsv')))
  for (mask in test_masks) {
    d <- merge(a[, .(gene_id, case = get(mask))],
               g[, .(gene_id, control = get(mask))], by = 'gene_id', all = TRUE)
    d[is.na(case), case := 0L]; d[is.na(control), control := 0L]
    d <- d[case >= 1 | control >= 1]
    if (!nrow(d)) stop('No genes to test: ', label, '/', mask)
    setorder(d, gene_id)
    d[, `:=`(nCase = cfg$n_case, nControl = cfg$n_control)]
    counts <- as.matrix(d[, .(case, nCase - case, control, nControl - control)])
    storage.mode(counts) <- 'integer'
    stopifnot(all(counts >= 0))
    support <- fisher.pvalues.support(counts, input = 'noassoc', alternative = 'two.sided')
    raw <- apply(counts, 1, function(x) fisher.test(matrix(x, 2), alternative = 'two.sided')$p.value)
    stopifnot(max(abs(raw - support$raw)) < 1e-12)
    d[, `:=`(af_label = label, mask = mask, raw_p = raw,
      Bonf_mask = p.adjust(raw, 'bonferroni'),
      ADBH_sd = ADBH(support$raw, support$support, direction = 'sd')$Adjusted,
      family_m = .N)]
    d <- cbind(d, odds_ratio_ci(d$case, d$nCase-d$case, d$control, d$nControl-d$control))
    d <- merge(d, genes, by = 'gene_id', all.x = TRUE, sort = FALSE)
    stopifnot(!anyNA(d$start))
    setorder(d, gene_id)
    tag <- paste(label, mask, sep = '_')
    saveRDS(list(gene_id = d$gene_id, support = support$support), file.path(out, paste0('support_', tag, '.rds')))
    results[[tag]] <- d
  }
}
write_tsv(rbindlist(results), file.path(out, 'all_results.tsv'))
record_session(out)
