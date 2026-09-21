#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R'); source('R/burden_plot_inputs.R')
cfg <- read_config()
out <- new_stage(cfg, 'FigureS5')
inputs <- prepare_burden_plots(cfg, include_power=TRUE)
plot <- with(inputs, {
s5<-arrange(list(powerpanel('AF0_01','PTV'),powerpanel('AF0_01','MIS'),powerpanel('AF0_01','DEL'),
 powerpanel('AF0_001','PTV'),powerpanel('AF0_001','MIS'),powerpanel('AF0_001','DEL')),LETTERS[1:6],ncol=3,nrow=2)

s5
})
save_figure(plot, file.path(out, 'FigureS5'), 13, 8)
record_session(out)
