#!/usr/bin/env Rscript
source('R/common.R')
cfg <- read_config()
out <- new_stage(cfg, 'power')
d <- fread(file.path(cfg$output_dir, 'burden/all_results.tsv'))
set.seed(cfg$seed)

# Genes with identical carrier counts share the same probability model. Simulate
# each distinct table once; these are post-hoc, nominal-significance estimates.
models <- unique(d[, .(case, control, nCase, nControl, OR)])
models[, control_probability := control/nControl]
pc <- pmin(pmax(models$control_probability, 1e-12), 1-1e-12)
models[, case_probability := OR*pc/(1-pc+OR*pc)]
power <- vapply(seq_len(nrow(models)), function(i) {
  m <- models[i]
  a <- rbinom(cfg$power_replicates, m$nCase, m$case_probability)
  c <- rbinom(cfg$power_replicates, m$nControl, m$control_probability)
  key <- paste(a, c, sep = ':'); selected <- match(unique(key), key)
  p <- vapply(selected, function(j) fisher.test(matrix(c(a[j], m$nCase-a[j],
      c[j], m$nControl-c[j]), 2), alternative = 'two.sided')$p.value, numeric(1))
  mean(p[match(key, unique(key))] < .05)
}, numeric(1))
models[, `:=`(power_nominal = power,
  power_MC_SE = sqrt(power*(1-power)/cfg$power_replicates), replicates = cfg$power_replicates)]
d <- merge(d, models, by = c('case', 'control', 'nCase', 'nControl', 'OR'), all.x = TRUE)
write_tsv(d, file.path(out, 'power.tsv'))
record_session(out)
