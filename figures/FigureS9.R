#!/usr/bin/env Rscript
source('R/common.R'); source('R/plot_style.R')
cfg <- read_config()
cowplot::set_null_device('pdf')
out <- new_stage(cfg, 'FigureS9')
long <- read_input(cfg,'protein_predictions',c('gene','variant','method','ddG','group'))
methods <- c('DUET','DynaMut2','I-Mutant2.0_PDB','I-Mutant2.0_SEQ','PremPS')
stopifnot(all(long$gene %chin% c('SOD1','TBK1')), all(long$method %chin% methods),
  all(long$group %chin% c('Known','New')),is.numeric(long$ddG))
long[,method:=factor(method,levels=methods)]
long[,status:=factor(group,levels=c('Known','New'))]
colors<-c(DUET='#CC3D24',DynaMut2='#F3C558','I-Mutant2.0_PDB'='#6DAE90','I-Mutant2.0_SEQ'='#30B4CC',PremPS='#004F7A')
orders<-list(SOD1=c('p.His44Arg','p.His49Arg','p.His81Arg','p.Ala90Thr','p.Gly142Glu','p.Glu50Val','p.Gly73Arg','p.His121Arg'),TBK1=c('p.Leu94Ser','p.Arg127Gln','p.Arg308Gln','p.Gln52Lys','p.Cys292Tyr'))
panel<-function(gene_name){
 x<-copy(long[gene==gene_name]);x[,variant:=factor(variant,levels=orders[[gene_name]])]
 p<-ggplot(x,aes(variant,ddG,fill=method))+
  geom_col(position=position_dodge(width=.9),width=.9)+
  geom_hline(yintercept=0,linetype='dashed',linewidth=.4)+
  facet_grid(.~status,scales='free_x',space='free_x')+
  scale_fill_manual(values=colors,drop=FALSE)+
  labs(title=paste('Effects on',gene_name,'protein stability'),x='Amino acid change',y=expression(Delta*Delta*G~'(kcal/mol)'),fill=NULL)+
  theme_bw(base_family='Arial',base_size=10)+
  theme(plot.title=element_text(size=12),axis.text=element_text(size=10,color='black'),
   axis.title=element_text(size=11),panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
   panel.grid.major.y=element_line(color='grey90',linewidth=.5),panel.grid.minor.y=element_blank(),panel.border=element_rect(color='black',linewidth=1),
   strip.background=element_rect(fill='grey88',color='black',linewidth=1),strip.text=element_text(size=10),
   legend.background=element_rect(color='black',linewidth=.5),legend.text=element_text(size=10),
   plot.margin=margin(6,6,6,18))
 built<-ggplot_build(p)$data[[1]]
 stopifnot(nrow(built)==nrow(x),isTRUE(all.equal(sort(built$y),sort(x$ddG))))
 p
}
p<-ggarrange(panel('SOD1'),panel('TBK1'),ncol=1,nrow=2,labels=c('A','B'),
 font.label=list(size=14,face='bold',family='Arial'),label.x=0,label.y=1,hjust=0,vjust=1,align='v')
save_figure(p,file.path(out,'FigureS9'),10,5.5,raster=FALSE)
record_session(out)
