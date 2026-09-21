#!/usr/bin/env Rscript
source('R/common.R')
cfg <- read_config()
out <- new_stage(cfg, 'counts')

to_wide <- function(x) {
  if (!nrow(x)) return(data.table(gene_id = character(), PTV = integer(), MIS = integer(), DEL = integer()))
  w <- dcast(x, gene_id ~ mask, value.var = 'carriers', fill = 0L)
  for (mask in c('PTV', 'MIS')) if (!mask %in% names(w)) w[, (mask) := 0L]
  # DEL follows the original additive definition in each cohort. A person with
  # both classes contributes to each class; this is not a cross-class union.
  w[, DEL := PTV + MIS]
  w[, .(gene_id, PTV, MIS, DEL)]
}

for (label in names(cfg$af_thresholds)) {
  a <- as.data.table(readRDS(file.path(cfg$output_dir, 'variants', paste0('cases_', label, '.rds'))))
  g <- as.data.table(readRDS(file.path(cfg$output_dir, 'variants', paste0('controls_', label, '.rds'))))
  a <- a[mask %chin% c('PTV', 'MIS')]; g <- g[mask %chin% c('PTV', 'MIS')]
  case_counts <- to_wide(a[, .(carriers = uniqueN(s)), by = .(gene_id, mask)])
  stopifnot(!anyDuplicated(g[, .(variant, gene_id, mask)]))
  g[, p0 := pmin(pmax((AN_gnomad/2 - AC_gnomad + nhomalt_gnomad)/(AN_gnomad/2), 0), 1)]
  control_counts <- to_wide(g[, .(carriers = as.integer(round(cfg$n_control *
    if (any(p0 == 0)) 1 else 1 - exp(sum(log(p0)))))), by = .(gene_id, mask)])
  stopifnot(all(as.matrix(case_counts[, -1]) <= cfg$n_case),
            all(as.matrix(control_counts[, -1]) <= cfg$n_control))
  write_tsv(case_counts, file.path(out, paste0('cases_', label, '.tsv')))
  write_tsv(control_counts, file.path(out, paste0('controls_', label, '.tsv')))
}
record_session(out)
