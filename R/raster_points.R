# Replace point grobs with transparent bitmaps; keep labels and axes as vectors.
rasterize_points <- function(plot, width, height, dpi = 600) {
  if (!requireNamespace('ragg', quietly = TRUE)) stop('Point rasterization requires ragg')
  g <- ggplot2::ggplotGrob(plot)
  points <- list()
  collect <- function(x) {
    if (inherits(x, 'points') && length(x$x)) points[[x$name]] <<- list(grob = x)
    if (!is.null(x$grobs)) invisible(lapply(x$grobs, collect))
    if (length(x$children)) invisible(lapply(x$children, collect))
  }
  collect(g)
  if (!length(points)) return(g)
  grDevices::cairo_pdf(tempfile(fileext = '.pdf'), width = width, height = height)
  grid::grid.newpage(); grid::grid.draw(g); grid::grid.force()
  listing <- grid::grid.ls(print = FALSE, viewports = TRUE, grobs = TRUE)
  for (name in names(points)) {
    i <- which(listing$name == name)
    stopifnot(length(i) == 1)
    path <- strsplit(listing$vpPath[i], '::', fixed = TRUE)[[1]]
    grid::upViewport(0)
    path <- path[path != 'ROOT']
    if (length(path)) grid::downViewport(do.call(grid::vpPath, as.list(path)))
    points[[name]]$width <- grid::convertWidth(grid::unit(1, 'npc'), 'in', valueOnly = TRUE)
    points[[name]]$height <- grid::convertHeight(grid::unit(1, 'npc'), 'in', valueOnly = TRUE)
    points[[name]]$xscale <- grid::current.viewport()$xscale
    points[[name]]$yscale <- grid::current.viewport()$yscale
  }
  grDevices::dev.off()
  for (name in names(points)) {
    spec <- points[[name]]
    if (!is.finite(spec$width) || !is.finite(spec$height) ||
        spec$width*dpi < 1 || spec$height*dpi < 1) next
    capture <- ragg::agg_capture(width = spec$width, height = spec$height,
      units = 'in', res = dpi, background = 'transparent')
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = spec$xscale, yscale = spec$yscale))
    grid::grid.draw(spec$grob)
    bitmap <- capture(); grDevices::dev.off()
    points[[name]]$raster <- grid::rasterGrob(bitmap, width = grid::unit(1, 'npc'),
      height = grid::unit(1, 'npc'), interpolate = FALSE, name = paste0(name, '_raster'))
  }
  replace <- function(x) {
    if (inherits(x, 'points') && !is.null(points[[x$name]]$raster)) return(points[[x$name]]$raster)
    if (!is.null(x$grobs)) x$grobs <- lapply(x$grobs, replace)
    if (length(x$children)) x$children <- do.call(grid::gList, lapply(x$children, replace))
    x
  }
  replace(g)
}
