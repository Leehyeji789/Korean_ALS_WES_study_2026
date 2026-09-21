prepare_burden_plots <- function(cfg, include_qq=FALSE, include_power=FALSE) {
d <- fread(file.path(cfg$output_dir, 'burden/all_results.tsv'))
if (include_power) z <- fread(file.path(cfg$output_dir, 'power/power.tsv'))
fdr <- file.path(cfg$output_dir,'qq')
if (include_qq) meta <- fread(file.path(fdr,'metadata.tsv'))
chr<-intersect(c(as.character(1:22),'X','Y','M','MT'),unique(d$chromosome))
lengths<-d[,.(len=max(start,na.rm=TRUE)),by=chromosome][match(chr,chromosome)]
lengths[,offset:=c(0,head(cumsum(len+5e6),-1))];lengths[,tick:=offset+len/2]
d<-merge(d,lengths[,.(chromosome,offset,tick)],by='chromosome',all.x=TRUE)
d[,x:=offset+start];d[,y:=-log10(pmax(raw_p,1e-300))]
d[,significance:=fifelse(Bonf_mask<.05,'Bonferroni',fifelse(ADBH_sd<.05,'ADBH only','Neither'))]
d[,background:=as.integer(factor(chromosome,levels=chr))%%2L]
source('R/burden_panels.R', local=TRUE)
environment()
}
