suppressPackageStartupMessages({library(ggplot2); library(ggpubr); library(ggrepel)})
# Use Cairo for off-screen label measurements as well as final export.
# The standard PDF device cannot measure Arial reliably on all systems.
options(device = function(...) grDevices::cairo_pdf(file = tempfile(fileext = '.pdf'), ...))
if (requireNamespace('extrafont', quietly = TRUE)) extrafont::loadfonts(device = 'pdf', quiet = TRUE)
cowplot::set_null_device(function(...) {
  grDevices::cairo_pdf(file = tempfile(fileext = '.pdf'), ...)
  grDevices::dev.control('enable')
})
source('R/raster_points.R')

af_title <- function(label) if (label == 'AF0_01') 'AF <0.01%' else 'AF <0.001%'
manuscript_theme <- function(font = 'Arial') {
  theme_bw() + theme(text = element_text(family = font, size = 10, colour = 'black'),
    plot.title = element_text(size = 13), axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = 'black'),
    panel.border = element_rect(colour = 'black', linewidth = 1, fill = NA),
    panel.grid.major = element_line(colour = 'grey90'),
    panel.grid.minor = element_blank(), plot.margin = margin(8, 10, 8, 25),
    legend.title = element_text(size = 11), legend.text = element_text(size = 10))
}

panel_grid <- function(plots, labels, ncol = 1, nrow = NULL, font = 'Arial') {
  if (is.null(nrow)) nrow <- ceiling(length(plots)/ncol)
  ggarrange(plotlist = plots, labels = labels, ncol = ncol, nrow = nrow,
    font.label = list(size = 16, face = 'bold', family = font),
    label.x = 0, label.y = 1, hjust = 0, vjust = 1, align = 'hv')
}

save_figure <- function(plot, path, width, height, raster = TRUE) {
  ggsave(paste0(path, '.pdf'), plot, width = width, height = height, device = cairo_pdf)
  if (raster) ggsave(paste0(path, '_raster600.pdf'), rasterize_points(plot,width,height),
                     width=width,height=height,device=cairo_pdf)
  if (nzchar(Sys.which('pdftoppm'))) {
    status <- system2('pdftoppm',c('-singlefile','-r','300','-png',
      shQuote(paste0(path,'.pdf')),shQuote(path)))
    stopifnot(status == 0L)
  } else ggsave(paste0(path,'.png'),plot,width=width,height=height,dpi=300)
}
