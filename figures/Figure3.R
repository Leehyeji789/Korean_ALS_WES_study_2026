#!/usr/bin/env Rscript
# Preserve the supplied final artwork, including coauthor edits.
source('R/common.R'); source('R/manuscript_artwork.R')
cfg <- read_config()
out <- new_stage(cfg, 'Figure3')
preserve_artwork(cfg, 'figure3_artwork', out, 'Figure3')
record_session(out)
