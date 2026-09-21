#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R'); source('R/burden_plot_inputs.R')
cfg <- read_config()
out <- new_stage(cfg, 'FigureS4')
inputs <- prepare_burden_plots(cfg, include_qq=TRUE)
plot <- with(inputs, {
s4<-ggarrange(arrange(list(ptv('AF0_01')),'A'),arrange(list(qq('AF0_01','PTV'),qq('AF0_01','MIS'),qq('AF0_01','DEL')),LETTERS[2:4],ncol=3),
 arrange(list(ptv('AF0_001')),'E'),arrange(list(qq('AF0_001','PTV'),qq('AF0_001','MIS'),qq('AF0_001','DEL')),LETTERS[6:8],ncol=3),
 ncol=1,nrow=4,heights=c(.8,1,.8,1))

s4
})
save_figure(plot, file.path(out, 'FigureS4'), 9.5, 12)
record_session(out)
