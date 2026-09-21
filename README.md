# Korean ALS WES study

R scripts for rare-variant burden analysis, gene-expression analyses and figures
in the Korean ALS whole-exome sequencing study.

The rare-variant burden analysis workflow was adapted from the
[CoCoRV framework](https://bitbucket.org/Wenan/cocorv/src/master/)
described in *Chen et al.*, 2022.

## Files

- `analysis/`: variant filtering, carrier counts, burden tests, QQ plots and power analysis.
- `analysis/expression/`: gene clustering, cell-type enrichment and GO analysis.
- `figures/`: scripts named by manuscript figure number.
- `R/`: shared functions.
- `config.example.R`: input paths and analysis settings.
