# Korean ALS WES study

R scripts for rare-variant burden analysis, gene-expression analyses and figures
in the Korean ALS whole-exome sequencing study.

## Files

- `analysis/`: variant filtering, carrier counts, burden tests, QQ plots and power analysis.
- `analysis/expression/`: gene clustering, cell-type enrichment and GO analysis.
- `figures/`: scripts named by manuscript figure number.
- `R/`: shared functions.
- `config.example.R`: input paths and analysis settings.

## Running the code

Copy `config.example.R` to `config.R` and set the input paths and output directory.
Run the scripts from the repository root:

```sh
Rscript analysis/01_variant_filtering.R config.R
Rscript analysis/02_sample_counts.R config.R
Rscript analysis/03_rare_variant_burden_test.R config.R
Rscript analysis/04_qq_lambda.R config.R
Rscript analysis/05_power.R config.R
Rscript figures/Figure1.R config.R
Rscript figures/FigureS4.R config.R
Rscript figures/FigureS5.R config.R
```

Run the expression scripts in numerical order. Other figure scripts take the same
`config.R` argument and require their corresponding inputs. Existing output
folders are not overwritten.

## Inputs and packages

The burden analysis starts from annotated variants and QC exports, not raw reads.
Expression analyses use processed expression and differential-expression data.
Participant-level data are not included.

Install the R packages loaded at the top of each script and its shared functions.
The burden analysis uses `data.table` and `DiscreteFDR`; figures use `ggplot2`,
`ggpubr` and related plotting packages. PDF export requires Cairo; manuscript
fonts use Arial.

Figures 2 and 3 use supplied final artwork. Figures S7 and S8 include R plots and
supplied Illustrator-edited artwork; these artwork files are not included.

The control-carrier estimation and discrete-testing methods follow CoCoRV
(Chen et al., 2022), as in the [original study code](https://github.com/Nahyeon203/Korean_ALS_WES).
