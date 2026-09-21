#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R'); source('R/burden_plot_inputs.R')
cfg <- read_config()
out <- new_stage(cfg, 'Figure1')
inputs <- prepare_burden_plots(cfg)
plot <- with(inputs, {
main<-arrange(list(mainpanel('AF0_01','MIS'),mainpanel('AF0_01','DEL'),
 mainpanel('AF0_001','MIS'),mainpanel('AF0_001','DEL')),LETTERS[1:4],ncol=2,nrow=2)


main
})
save_figure(plot, file.path(out, 'Figure1'), 13, 6.3)
record_session(out)
