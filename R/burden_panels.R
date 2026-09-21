# Panel definitions retained from the final manuscript scripts.
aftext<-function(af)if(af=='AF0_01')'AF <0.01%' else 'AF <0.001%'
label_style<-list(size=16,face='bold',family='Arial',color='black')
arrange<-function(plots,labels,ncol=1,nrow=NULL,heights=NULL){
 ggarrange(plotlist=plots,ncol=ncol,nrow=nrow,labels=labels,font.label=label_style,
 label.x=0,label.y=1,hjust=0,vjust=1,heights=if(is.null(heights))1 else heights,align='hv')
}

paneltheme<-theme(plot.margin=margin(8,10,8,25),plot.caption=element_blank())
theme_als<-theme(
 plot.title=element_text(size=13,family='Arial',color='black'),
 text=element_text(size=10,family='Arial',color='black'),
 axis.title=element_text(size=11,family='Arial',color='black'),
 axis.text=element_text(size=10,family='Arial',color='black'),
 legend.title=element_text(size=11,family='Arial',color='black'),
 legend.text=element_text(size=10,family='Arial',color='black'),
 legend.background=element_rect(color='black',linewidth=.5),
 legend.box='vertical',legend.margin=margin(3,3,3,3),
 panel.border=element_rect(color='black',linewidth=1,fill=NA),
 panel.background=element_blank(),panel.grid.major=element_line(color='grey90'),
 strip.background=element_rect(color='black',linewidth=1,fill='grey88'),
 strip.text=element_text(size=10,family='Arial',color='black'))
manhattan_theme<-theme_bw()+theme(
 title=element_text(size=11,color='black'),plot.title=element_text(vjust=.6),
 axis.title=element_text(size=11,color='black'),axis.text=element_text(size=10,color='black'),
 legend.position='none',panel.border=element_rect(color='black',linewidth=1),
 panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
 panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+paneltheme
mainpanel<-function(af,variant_mask){
 full<-d[af_label==af & mask==variant_mask];shown<-full
 labels<-shown[significance!='Neither']
 ggplot(shown,aes(x,y))+
  geom_point(aes(color=factor(background)),alpha=1,size=1.5)+
  geom_point(data=labels,aes(color=significance),size=2)+
  geom_hline(yintercept=-log10(.05/nrow(full)),linetype=2,color='#b63746',linewidth=.35)+
  geom_text_repel(data=labels,aes(label=gene_symbol),box.padding=.2,
   max.overlaps=Inf,show.legend=FALSE,size=12*.3514321,seed=52)+
  scale_color_manual(values=c('1'='#1f78b4','0'='#a6cee3','ADBH only'='#e57d32','Bonferroni'='#b63746'),guide='none')+
  scale_x_continuous(breaks=lengths$tick,labels=chr)+scale_y_continuous(expand=c(0,.05))+
  coord_cartesian(ylim=c(0,ceiling(max(shown$y)+.3)))+
  labs(title=paste0(variant_mask,', ',aftext(af)),x='Chromosome',y=expression(-log[10](italic(p))))+
  manhattan_theme
}
qq<-function(af,variant_mask){
 q<-fread(file.path(fdr,paste0('QQ_',af,'_',variant_mask,'.tsv')))
 stat<-meta[af_label==af & mask==variant_mask]
 stopifnot(nrow(stat)==1,nrow(q)==nrow(d[af_label==af & mask==variant_mask]))
 lim<-ceiling(max(q$expected,q$observed))
 p<-ggplot(q,aes(expected,observed))+
  geom_point(color='#004F7AFF',size=2)+geom_abline(intercept=0,slope=1)+
  coord_fixed(xlim=c(0,lim),ylim=c(0,lim))+
  labs(title=paste0(variant_mask,' QQ plot, ',aftext(af)),
   x=expression(Expected~-log[10]*italic(P)),y=expression(Observed~-log[10]*italic(P)))+
  theme_bw()+theme_als+theme(plot.title=element_text(vjust=.6),
   panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),
   panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())+paneltheme
 if(is.finite(stat$lambda_emp))p<-p+
  geom_abline(intercept=stat$intercept,slope=stat$lambda_emp,linetype='dashed')+
  annotate('text',x=lim,y=0,label=sprintf('lambda[emp] == %.2f',stat$lambda_emp),
   parse=TRUE,hjust=1,vjust=0,size=3.5,family='Arial')
 p
}
ptv<-function(af){
 z<-d[af_label==af & mask=='PTV']
 threshold<- -log10(.05/nrow(z))
 upper<-ceiling(max(z$y,threshold)+.5)
 message(af, ': Bonferroni line = ', signif(threshold,6), '; y-axis upper limit = ',upper)
 ggplot(z,aes(x,y,color=factor(background)))+geom_point(alpha=1,size=1.5)+
  geom_hline(yintercept=-log10(.05/nrow(z)),linetype='dashed',color='darkred',linewidth=.35)+
  scale_color_manual(values=c('1'='#1f78b4','0'='#a6cee3'),guide='none')+
  scale_x_continuous(breaks=lengths$tick,labels=chr)+scale_y_continuous(expand=c(0,.05))+
  coord_cartesian(ylim=c(0,upper))+
  labs(title=paste0('PTV burden, ',aftext(af)),x='Chromosome',y=expression(-log[10](italic(p))))+
  manhattan_theme
}
powerpanel<-function(af,variant_mask){
 shown<-z[af_label==af & mask==variant_mask];hl<-shown[ADBH_sd<.05 & power_nominal>=.8]
 ggplot(shown,aes(-log10(pmax(raw_p,1e-300)),power_nominal))+
  geom_point(data=shown[!(ADBH_sd<.05 & power_nominal>=.8)],color='black',size=.7,alpha=.7)+geom_jitter(data=hl,color='blue',size=2,position=position_jitter(width=0,height=0,seed=52))+
  geom_text_repel(data=hl,aes(label=gene_symbol),color='blue',size=3,seed=52,max.overlaps=50,segment.size=.3,segment.color='blue',box.padding=.4,point.padding=.2)+
  geom_vline(xintercept=-log10(.05),color='darkred',linetype='dashed',linewidth=.3)+
  geom_hline(yintercept=.8,color='darkred',linetype='dashed',linewidth=.3)+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
  scale_x_continuous(limits=c(0,max(-log10(pmax(z$raw_p,1e-300)))*1.05))+
  labs(title=paste0(variant_mask,', ',aftext(af)),x='-log10(P-value)',y='Post-hoc power')+
  theme_gray()+theme_als+paneltheme
}
