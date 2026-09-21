suppressPackageStartupMessages(library(data.table))

read_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop('Supply one configuration file, e.g. config.R')
  path <- normalizePath(args[1], mustWork = TRUE)
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  cfg <- env$config
  resolve <- function(p) {
    if (is.null(p)) return(NULL)
    if (grepl('^(/|[A-Za-z]:)', p)) p else file.path(dirname(path), p)
  }
  cfg$output_dir <- resolve(cfg$output_dir)
  cfg$inputs <- lapply(cfg$inputs, resolve)
  stopifnot(cfg$n_case > 0, cfg$n_control > 0, cfg$call_rate > 0,
            cfg$call_rate <= 1, cfg$null_replicates >= 2, cfg$power_replicates >= 1)
  cfg
}

new_stage <- function(cfg, name) {
  path <- file.path(cfg$output_dir, name)
  if (file.exists(path)) stop('Output already exists: ', path)
  dir.create(path, recursive = TRUE)
  if (!dir.exists(path)) stop('Cannot create output directory: ', path)
  path
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, ' is missing: ', paste(missing, collapse = ', '))
}

read_input <- function(cfg, name, columns = character()) {
  path <- cfg$inputs[[name]]
  if (is.null(path) || !file.exists(path)) stop('Missing input: ', name)
  x <- fread(path)
  require_columns(x, columns, name)
  x
}

write_tsv <- function(x, path) fwrite(x, path, sep = '\t', na = 'NA')
record_session <- function(path) writeLines(capture.output(sessionInfo()), file.path(path, 'sessionInfo.txt'))

rankscore_columns <- c('revel_rankscore', 'metasvm_rankscore', 'gmvp_rankscore',
  'mpc_rankscore', 'primateai_rankscore', 'metarnn_rankscore',
  'fathmm_converted_rankscore', 'provean_converted_rankscore', 'vest4_rankscore',
  'deogen2_rankscore', 'clinpred_rankscore', 'varity_r_rankscore')

variant_masks <- c(protein_truncating = 'PTV', missense = 'MIS', synonymous = 'SYN')
test_masks <- c('PTV', 'MIS', 'DEL')

odds_ratio_ci <- function(a, b, c, d) {
  a <- a + .5; b <- b + .5; c <- c + .5; d <- d + .5
  odds_ratio <- a * d / (b * c)
  se <- sqrt(1/a + 1/b + 1/c + 1/d)
  data.table(OR = odds_ratio, CI_lower = exp(log(odds_ratio) - 1.96 * se),
             CI_upper = exp(log(odds_ratio) + 1.96 * se))
}
