#!/usr/bin/env Rscript
source('R/common.R')
cfg <- read_config()

# Inputs have passed shared depth-region and upstream sample/variant QC.
# The case key list records the pre-genotype-filter call-rate definition.
annotation_columns <- c('variant', 'gene_id', 'VariantType', 'biotype',
  'kova_AF', 'jpn_AF', 'allofus250k__gvs_max_af', 'cadd_phred', rankscore_columns)
cases <- read_input(cfg, 'cases', c('s', 'internal_AC', 'internal_AN', annotation_columns))
controls <- read_input(cfg, 'controls', c('AC_gnomad', 'AN_gnomad', 'nhomalt_gnomad', annotation_columns))
samples <- read_input(cfg, 'sample_ids', 's')
stopifnot(!anyDuplicated(samples$s), nrow(samples) == cfg$n_case,
          all(cases$s %chin% samples$s))
case_keys <- as.character(read_input(cfg, 'case_pass_keys', 'variant')$variant)
control_exclusions <- as.character(read_input(cfg, 'control_exclusion_keys', 'variant')$variant)
site <- read_input(cfg, 'control_callability', c('variant', 'gnomad_eas_an'))
all_sites <- read_input(cfg, 'all_sites_an', c('locus', 'EAS_AN'))
stopifnot(!anyDuplicated(site$variant), !anyDuplicated(all_sites$locus))
min_an <- ceiling(2 * cfg$n_control * cfg$call_rate)

# A failed allele removes that same allele in both cohorts, not every alternate
# allele at the position. The all-sites filter below is deliberately position-based.
cases <- cases[variant %chin% case_keys & !variant %chin% control_exclusions]
controls <- controls[!variant %chin% control_exclusions &
                       !is.na(AN_gnomad) & AN_gnomad >= min_an]
failed_alleles <- site[is.na(gnomad_eas_an) | gnomad_eas_an < min_an, variant]
cases <- cases[!variant %chin% failed_alleles]
controls <- controls[!variant %chin% failed_alleles]
position_filter <- function(x) {
  locus <- sub('^([^:]+:[^:]+):.*$', '\\1', x$variant)
  x[, EAS_all_sites_AN := all_sites$EAS_AN[match(locus, all_sites$locus)]]
  x[!is.na(EAS_all_sites_AN) & EAS_all_sites_AN >= min_an]
}
cases <- position_filter(cases)
controls <- position_filter(controls)
stopifnot(all(cases$internal_AC >= 0 & cases$internal_AC <= cases$internal_AN),
          all(controls$AC_gnomad >= 0 & controls$AC_gnomad <= controls$AN_gnomad),
          all(2 * controls$nhomalt_gnomad <= controls$AC_gnomad))

out <- new_stage(cfg, 'variants')
summary <- list()
for (label in names(cfg$af_thresholds)) {
  cutoff <- cfg$af_thresholds[[label]]
  external_filter <- function(x) {
    x <- copy(x)[(is.na(kova_AF) | kova_AF < cutoff) &
                  (is.na(jpn_AF) | jpn_AF < cutoff) &
                  (is.na(allofus250k__gvs_max_af) | allofus250k__gvs_max_af < cutoff)]
    x[, mean_rankscore := rowMeans(.SD, na.rm = TRUE), .SDcols = rankscore_columns]
    x[biotype == 'protein_coding' & VariantType %chin% names(variant_masks)]
  }
  a <- external_filter(cases)
  g <- external_filter(controls)[AC_gnomad > 0]
  ac <- unique(a[, .(variant, case_AC = internal_AC, case_AN = internal_AN)])
  gc <- unique(g[, .(variant, control_AC = AC_gnomad, control_AN = AN_gnomad)])
  stopifnot(!anyDuplicated(ac$variant), !anyDuplicated(gc$variant))
  joint <- merge(ac, gc, by = 'variant', all = TRUE)
  for (column in setdiff(names(joint), 'variant')) set(joint, which(is.na(joint[[column]])), column, 0)
  keep <- joint[(case_AC + control_AC)/(case_AN + control_AN) < .01, variant]
  qualifying <- function(x) {
    x <- x[variant %chin% keep & (VariantType != 'missense' |
                                  (!is.na(mean_rankscore) & mean_rankscore >= .7 & cadd_phred >= 20))]
    x[, mask := unname(variant_masks[VariantType])]
    x
  }
  a <- qualifying(a); g <- qualifying(g)
  saveRDS(a, file.path(out, paste0('cases_', label, '.rds')))
  saveRDS(g, file.path(out, paste0('controls_', label, '.rds')))
  summary[[label]] <- data.table(af_label = label, case_rows = nrow(a),
    case_variants = uniqueN(a$variant), control_variants = uniqueN(g$variant))
}
write_tsv(rbindlist(summary), file.path(out, 'filter_summary.tsv'))
record_session(out)
