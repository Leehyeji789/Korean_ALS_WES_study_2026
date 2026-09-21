#!/usr/bin/env Rscript
# Curated post-QC pathogenic-variant carriers, not all burden-test variants.
source('R/common.R')
suppressPackageStartupMessages(library(ggplot2))
cfg <- read_config()
x <- read_input(cfg, 'known_pathogenic_carriers', c('gene_symbol', 's'))
stopifnot(!anyNA(x$s), !anyNA(x$gene_symbol))
out <- new_stage(cfg, 'FigureS3')
output <- file.path(out, 'FigureS3.pdf')
counts<-x[,.(carriers=uniqueN(s)),by=gene_symbol][order(-carriers)]
screened<-uniqueN(x$s); if(sum(counts$carriers)!=screened) stop('Multi-gene carriers require explicit mutually exclusive categories'); stopifnot(screened <= cfg$n_case); counts<-rbind(counts,data.table(gene_symbol="No screened pathogenic variant",carriers=cfg$n_case-screened))
counts[,label:=sprintf("%s: %d",gene_symbol,carriers)]
p<-ggplot(counts,aes(x=2,y=carriers,fill=label))+geom_col(color="white")+coord_polar(theta="y")+xlim(.5,2.5)+
 annotate("text",x=.5,y=0,label=sprintf("%d/%s\n(%.2f%%)",screened,format(cfg$n_case,big.mark=",",scientific=FALSE),100*screened/cfg$n_case),size=5)+
 labs(title="Previously reported pathogenic variants after QC95",fill=NULL)+theme_void(base_size=9)+theme(legend.position="bottom")
ggsave(output,p,width=7,height=5.5,useDingbats=FALSE)

record_session(out)
