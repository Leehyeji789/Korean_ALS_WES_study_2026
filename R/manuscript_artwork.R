preserve_artwork <- function(cfg, input_name, out, name) {
  input <- cfg$inputs[[input_name]]
  if (is.null(input) || !file.exists(input)) stop('Supply the final artwork: ', input_name)
  output <- file.path(out, paste0(name, '.', tools::file_ext(input)))
  stopifnot(!file.exists(output), file.copy(input, output, overwrite = FALSE))
  stopifnot(identical(unname(tools::md5sum(input)), unname(tools::md5sum(output))))
}
