# Korean ALS whole-exome sequencing

Code for the rare-variant burden analysis and associated analyses in the Korean
ALS whole-exome sequencing study. This directory is a publication copy of the
analysis code; it does not contain participant data.

## Directory layout

- `analysis/`: burden analysis and expression analysis.
- `figures/`: one script per figure (`Figure1.R`, `Figure2.R`, `Figure3.R`,
  `FigureS2.R`, `FigureS3.R`, `FigureS4.R`, `FigureS5.R`, `FigureS6.R`, `FigureS7.R`, `FigureS8.R`, `FigureS9.R`).
- `R/`: shared calculations, panel definitions and formatting functions.

Run each script from the repository root. Each figure script creates only
its own output directory, such as `results/Figure1/`, and
refuses to overwrite it. Figure 1 does not require QQ or power output; S4 needs
QQ output and S5 needs power output. S2 uses post-QC ancestry coordinates; S3 uses
curated pathogenic-variant carriers; S6 plots cached silhouette scores. See
[inputs and commands for S2, S3 and S6](#figures-s2-s3-and-s6).

## Start here

Run commands from this directory. The numbered scripts retain the order of the
original repository: filter variants, estimate carrier counts, then test genes.
Configuration is kept in one file rather than in the analysis scripts.

```sh
cp config.example.R config.R
# Edit the input paths and output directory in config.R.
Rscript analysis/01_variant_filtering.R config.R
Rscript analysis/02_sample_counts.R config.R
Rscript analysis/03_rare_variant_burden_test.R config.R
Rscript analysis/04_qq_lambda.R config.R
Rscript analysis/05_power.R config.R
Rscript figures/Figure1.R config.R
Rscript figures/FigureS4.R config.R
Rscript figures/FigureS5.R config.R
```

Paths in the configuration are relative to the configuration file, unless they
are absolute. Each script creates its own output subdirectory and refuses to
overwrite an existing stage. Use a new output directory for a new run.

The analysis starts from annotated variants and upstream QC exports, not raw
sequencing reads. Read [input formats](#input-files) before using study data.
Participant-level inputs and intermediate results must remain private.

## What the scripts produce

| Script | Output |
| --- | --- |
| `analysis/01_variant_filtering.R` | Eligible variants at both external AF thresholds |
| `analysis/02_sample_counts.R` | Observed case counts and estimated control counts |
| `analysis/03_rare_variant_burden_test.R` | Two-sided Fisher P-values, Bonferroni and ADBH results, ORs and CIs |
| `analysis/04_qq_lambda.R` | Resampled-null QQ coordinates and inflation estimates |
| `analysis/05_power.R` | Observed-effect power at nominal P <0.05 |
| `figures/Figure1.R` | Figure 1 |
| `figures/FigureS4.R` | Figure S4 |
| `figures/FigureS5.R` | Figure S5 |

Figures are generated in R. Manuscript panel assembly and
externally supplied results are described separately in [scope](#scope-of-this-code-release).
See [statistical details](#counting-and-statistical-testing) for the counting and testing rules.

## Software

The burden workflow uses R, `data.table`, `DiscreteFDR`, `ggplot2`, `ggpubr`,
`ggrepel`. Figure export requires Cairo support and the Arial
font to reproduce manuscript typography. Point rasterization additionally uses
`ragg`. Each stage records `sessionInfo.txt`; this is an environment record,
not a claim that package versions are locked.

Expression analyses additionally use `rtracklayer`, `ComplexHeatmap`, `circlize`,
`RColorBrewer`, `cluster`, `clusterProfiler`, `org.Hs.eg.db`, `dplyr`, `enrichplot`,
`stringr`, and `scales`. `extrafont` supplies the legacy font metrics used for S9.

## Expression and functional analyses

With the corresponding inputs configured, run:

```sh
Rscript analysis/expression/01_gene_clusters.R config.R
Rscript analysis/expression/02_cell_type_enrichment.R config.R
Rscript analysis/expression/03_gene_ontology.R config.R
Rscript figures/FigureS7.R config.R
Rscript figures/FigureS8.R config.R
Rscript figures/FigureS9.R config.R
Rscript figures/Figure2.R config.R
Rscript figures/Figure3.R config.R
```

Figure 1, S4, S5 and S9 use the final manuscript panel definitions. Burden plots
are saved as vector PDFs and as separate `_raster600.pdf` files with only the
points rasterized. When available, `pdftoppm` produces PNGs from the vector PDF.

Figures 2 and 3 contain coauthor-edited artwork; S7 and S8 were also edited in
Illustrator after R export. The individual figure scripts preserve the supplied final
files without changing their contents. The S7/S8 plotting scripts retain
the original R panels, including S8's comparative panel A and group panels B–E,
but do not claim to recreate the subsequent Illustrator edits. These supplied
artwork files are not bundled with the code.

The heatmap uses the 51 ClinGen ALS genes classified as Definitive, Strong,
Moderate or Limited, plus CDK15. It does not perform separate analyses by evidence
level. Groups A–D follow the displayed heatmap order. Cell-type enrichment uses
the protein-coding genes actually tested in the DEG analysis as its background;
GO enrichment retains the original package-default annotated-gene background.

`DiscreteFDR::fisher.pvalues.support()` is retained for compatibility with the
study calculation. Recent package versions mark it deprecated. Raw P-values
are checked against `stats::fisher.test()` before adjustment.

## Provenance and reuse

The original submitted code is available at
[Korean_ALS_WES](https://github.com/Nahyeon203/Korean_ALS_WES).
The control-carrier estimation and discrete-testing workflow follow the
CoCoRV framework (Chen et al., 2022), as acknowledged in the original repository.
Changes to that workflow are listed in [methods](#counting-and-statistical-testing).

No software license has been added here. The repository owner should choose
one and check any third-party code terms before redistribution.

## Input files

All text tables are tab-separated, with a header. Gzip compression is supported.
Missing numeric annotations should be `NA`, not zero. Coordinates are GRCh38
unless stated otherwise. Normalize and split multiallelic variants upstream;
all variant keys must use the same representation, `chr:position:ref:alt`.

### Annotated variants

`cases` contains one row per carrier, variant and gene. `controls` contains one
row per variant and gene, using the gnomAD EAS adjusted counts. Both need:

| Columns | Meaning |
| --- | --- |
| `variant`, `gene_id` | Full allele key and Ensembl gene ID |
| `VariantType` | `protein_truncating`, `missense` or `synonymous` |
| `biotype` | Gene biotype; the analysis retains `protein_coding` |
| `kova_AF`, `jpn_AF`, `allofus250k__gvs_max_af` | External frequencies; All of Us uses the population maximum |
| `cadd_phred` | CADD PHRED score |
| Rank-score columns in `R/common.R` | The 12 scores used for the missense mask |

Cases also need `s`, `internal_AC`, and `internal_AN`. AC and AN are site-level
values repeated across carrier rows, not counts to sum over those rows.
Controls need `AC_gnomad`, `AN_gnomad`, and `nhomalt_gnomad`.

### QC exports

These inputs are required, not optional substitutes for upstream QC.

| Configuration key | Required columns and definition |
| --- | --- |
| `sample_ids` | `s`: all retained cases, including noncarriers |
| `case_pass_keys` | `variant`: alleles passing case variant QC, including the pre-genotype-filter call-rate threshold |
| `control_exclusion_keys` | `variant`: alleles failing reciprocal QC; a header-only file is valid if none failed |
| `control_callability` | `variant`, `gnomad_eas_an`: variant-specific EAS adjusted AN across queried alleles |
| `all_sites_an` | `locus`, `EAS_AN`: position keys `chr:position` and EAS adjusted AN from the all-sites resource |
| `genes` | `gene_id`, `gene_symbol`, `chromosome`, `start`, `end`; one row per gene; chromosome without `chr` |

The shared depth-region restriction and sample/genotype QC must already have
been applied upstream. An absent variant-specific control record is not itself
an AN failure, because a case allele can be absent from gnomAD. Missing all-sites
AN, in contrast, excludes the position in both datasets.

The study uses 1,141 cases and 19,850 controls. EAS adjusted AN must be at least
37,715 (95% of 39,700). Sample-level call rate remains >80%; case variant-level
call rate is at least 95%. Do not recompute the latter from a different
post-genotype-filter checkpoint and assume it is the same filter.

### Expression and external analyses

Expression inputs are a supercluster-by-gene mean z-score RDS, the matching GTF,
the ClinGen ALS GCEP CSV (`Gene`, `Classification`), and supercluster DEG results
(`names`, `group`, `pvals_adj`, `logfoldchanges`, `pct_nz_group`). Gene IDs in the
matrix and DEG file must match the GTF. Use the study annotation releases.

Protein prediction plots use a long-form table with `gene`, `variant`, `method`,
`ddG`, and `group`. Groups are manuscript-defined `Known` and `New`; the plotting
code does not infer clinical classification or novelty from current databases.

Do not commit participant-level inputs, credentials, notebook outputs, or
generated carrier workbooks. The ignore file is a precaution, not a privacy
review.

## Figures S2, S3 and S6

Run each script from the repository root with one configuration argument:

```sh
Rscript figures/FigureS2.R config.R
Rscript figures/FigureS3.R config.R
Rscript figures/FigureS6.R config.R
```

Each writes its own FigureS2, FigureS3 or FigureS6 subdirectory and refuses to overwrite it.

| Figure | Configuration input | Required columns |
| --- | --- | --- |
| S2 | `ancestry` | `#sample_id`, `predicted_ancestry`, `PC1`, `PC2`, `PC3` |
| S2 | `ancestry_labels` | `#sample` (1000 Genomes reference IDs) |
| S3 | `known_pathogenic_carriers` | `gene_symbol`, `s` (case ID) |
| S6 | `silhouette` (optional) | `k`, `silhouette` |

S2 uses post-QC Somalier output. Reference IDs identify 1000 Genomes samples;
other samples are ALS cases. Original panel colors, symbols and theme are retained.
The combined layout is not verified as an exact match to the edited manuscript PDF.

S3 uses the curated post-QC pathogenic-variant carrier table, not the full burden
input. Multi-gene carriers trigger an error because pie categories must be mutually
exclusive. Its plot code retains the final revision donut chart.

S6 defaults to `expression_clusters/silhouette.tsv` under the output directory.
An existing score table can instead be supplied as `inputs$silhouette`. It retains
the original plot and does not repeat clustering. Keep individual-level inputs private.

## Counting and statistical testing

### Filters

The external frequency thresholds are **0.01% and 0.001%**, represented as
`0.0001` and `0.00001` in code. At each threshold, KOVA, JPN and the maximum
All of Us frequency must be below the cutoff or missing. The pooled case/control
frequency must additionally be below 1%, following the original filtering order.
These external thresholds are not filters on the observed case frequency alone.

Missense variants require a mean rank score of at least 0.7 across available
scores and CADD PHRED at least 20. PTV and qualifying missense variants are
retained separately before constructing DEL.

The original shared high-depth regions are retained (at least 10× depth in at
least 90% of samples). The gnomAD depth summary is for the overall exome cohort,
not EAS alone. EAS all-sites adjusted AN at least 95% is an additional
callability criterion, not an ancestry-specific depth measurement. It is applied
to positions in both case and control datasets. Variant-level reciprocal
exclusions use the full allele key, not every allele at a position.

### Carrier counts

Within each gene and each of PTV and MIS, cases are counted as distinct sample
IDs. For controls, the noncarrier probability for an allele is

```
p0 = (AN/2 - AC + nhomalt) / (AN/2)
```

The estimated gene carrier count is `round(n_control * (1 - product(p0)))`.
This aggregate-data approximation does not reconstruct individual genotypes or
linkage between alleles. **DEL is the sum of PTV and MIS counts in each cohort**,
as in the original pipeline; it is not a cross-class deduplicated carrier union.

### Tests and effect sizes

Each AF threshold and each of PTV, MIS and DEL defines a separate test family.
A gene enters a family if it has at least one observed case carrier or at least
one estimated control carrier after rounding. There is no OR-direction filter.

Two-sided Fisher tests use the uncorrected 2×2 counts. Bonferroni correction and
step-down ADBH use all genes in that family. ADBH uses conditional discrete
two-sided Fisher null supports from `DiscreteFDR`. Raw Fisher P-values from that
calculation are checked against `stats::fisher.test()` to tolerance 1e-12.
RBH is not calculated.

For reported ORs and 95% log-Wald CIs, 0.5 is added to each of the four cells.
This correction is not applied to Fisher tests. Noncarrier counts are cohort
size minus carrier count, subtracted once.

### QQ plots and power

QQ expectations come from 10,000 resampled datasets using the discrete null
supports. Inflation is the fitted observed-versus-expected −log10 P slope with
an intercept, excluding the most significant 5% from the fit. All genes remain
in the QQ plot. Resampled inflation estimates use leave-one-out null
expectations. This is the study's regression-based measure, not the median
chi-square definition of genomic inflation.

Observed-effect power uses 10,000 binomial simulations, the observed control
carrier frequency, and the corrected observed OR. A zero control frequency
remains zero for control draws; it is bounded away from zero only for converting
the OR to a case probability. Significance is nominal two-sided Fisher P <0.05,
not Bonferroni significance. Both OR directions are plotted. These are post-hoc
observed-effect estimates, not independent evidence of association.

### Changes from the submitted repository

The original three-step structure and control-carrier approximation are kept.
The revision adds the agreed QC exports, two-sided rather than enrichment-only
Fisher tests, case-or-control rather than case-only testing families, Bonferroni
results, and the separate diagnostic and presentation scripts. No result is
filtered out merely because its estimated OR is below one.

## Scope of this code release

The numbered burden scripts form an executable workflow from annotated,
QC-filtered variants to statistics, diagnostics and figures.
Raw-read alignment, joint calling, VEP annotation and the production Hail QC
checkpoints are upstream inputs; this directory is not yet an end-to-end
raw-read workflow.

Expression and functional-analysis scripts are kept separately from the burden
workflow. They consume the study's processed expression matrix and DEG results;
they do not recreate those inputs from raw single-cell reads.

Externally generated protein predictions, curated clinical classifications,
clinical summaries and coauthor-supplied illustrations require their original
inputs. Their presence in the manuscript does not imply that a plotting script
recomputes the underlying experiments or curation.

Before publication, check the manuscript output inventory against this directory.
In particular, a replacement R plot is not evidence of exact reproduction of a
coauthor's final assembled artwork. The full original manuscript and private
workbooks are intentionally not bundled here.
