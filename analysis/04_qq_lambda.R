#!/usr/bin/env Rscript
source('R/common.R')
cfg <- read_config()
out <- new_stage(cfg, 'qq')
results <- fread(file.path(cfg$output_dir, 'burden/all_results.tsv'))

# CoCoRV's empirical lambda is a regression slope, not a chi-square median ratio.
estimate_lambda <- function(p, expected_p) {
  p[abs(p-1) < 1e-12] <- 1
  expected_p[abs(expected_p-1) < 1e-12] <- 1
  observed <- sort(-log10(pmax(p, 1e-30)))
  expected <- -log10(expected_p)
  use <- seq_len(round(length(p) * .95))
  if (nrow(unique(cbind(observed[use], expected[use]))) <= 20)
    return(c(intercept = NA_real_, lambda = NA_real_))
  fit <- coef(lm(observed[use] ~ expected[use]))
  c(intercept = unname(fit[1]), lambda = unname(fit[2]))
}

metadata <- list()
for (label in names(cfg$af_thresholds)) for (mask_name in test_masks) {
  tag <- paste(label, mask_name, sep = '_')
  d <- results[af_label == label & mask == mask_name]; setorder(d, gene_id)
  saved <- readRDS(file.path(cfg$output_dir, 'burden', paste0('support_', tag, '.rds')))
  stopifnot(identical(d$gene_id, saved$gene_id))
  set.seed(cfg$seed)
  null <- matrix(NA_real_, nrow(d), cfg$null_replicates)
  for (i in seq_len(nrow(d))) {
    s <- saved$support[[i]]
    null[i, ] <- sample(s, cfg$null_replicates, replace = TRUE, prob = diff(c(0, s)))
  }
  null <- apply(null, 2, sort, decreasing = TRUE)
  if (is.null(dim(null))) null <- matrix(null, nrow = nrow(d))
  means <- rowMeans(null)
  fit <- estimate_lambda(d$raw_p, means)
  total <- rowSums(null)
  null_lambda <- vapply(seq_len(ncol(null)), function(j)
    estimate_lambda(null[, j], (total-null[, j])/(ncol(null)-1))['lambda'], numeric(1))
  write_tsv(data.table(expected = -log10(means), observed = sort(-log10(pmax(d$raw_p, 1e-30)))),
            file.path(out, paste0('QQ_', tag, '.tsv')))
  write_tsv(data.table(replicate = seq_along(null_lambda), lambda_null = null_lambda),
            file.path(out, paste0('lambda_null_', tag, '.tsv')))
  metadata[[tag]] <- data.table(af_label = label, mask = mask_name, genes = nrow(d),
    intercept = fit['intercept'], lambda_emp = fit['lambda'], replicates = cfg$null_replicates)
  rm(null); gc(verbose = FALSE)
}
write_tsv(rbindlist(metadata), file.path(out, 'metadata.tsv'))
record_session(out)
