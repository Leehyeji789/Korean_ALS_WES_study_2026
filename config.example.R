# Copy this file to config.R and edit the input paths. Paths are relative to
# this file unless absolute. Individual-level data are not distributed here.
config <- list(
  output_dir = 'results',
  n_case = 1141L,
  n_control = 19850L,
  call_rate = 0.95,
  af_thresholds = c(AF0_01 = 0.0001, AF0_001 = 0.00001),
  null_replicates = 10000L,
  power_replicates = 10000L,
  seed = 52L,
  font_family = 'Arial',
  inputs = list(
    cases = 'inputs/ALS_WES.tsv.gz',
    controls = 'inputs/gnomad_exomes.tsv.gz',
    sample_ids = 'inputs/case_samples.tsv',
    case_pass_keys = 'inputs/case_qc95_keys.tsv.gz',
    control_exclusion_keys = 'inputs/control_exclusion_keys.tsv.gz',
    control_callability = 'inputs/control_callability.tsv.gz',
    all_sites_an = 'inputs/eas_all_sites_an.tsv.gz',
    genes = 'inputs/gene_annotations.tsv',
    expression = 'inputs/supercluster_mean_zscore.rds',
    expression_gtf = 'inputs/gb_pri_annot_filtered.gtf.gz',
    clingen = 'inputs/ClinGen_ALS_GCEP.csv',
    differential_expression = 'inputs/supercluster_DEGs.csv',
    ancestry = 'inputs/somalier-ancestry.somalier-ancestry.tsv',
    ancestry_labels = 'inputs/ancestry-labels-1kg.tsv',
    known_pathogenic_carriers = 'inputs/known_pathogenic_carriers.tsv',
    protein_predictions = 'inputs/protein_predictions.tsv',
    figure2_artwork = 'inputs/Figure2.docx',
    figure3_artwork = 'inputs/Figure3.docx',
    figureS7_artwork = 'inputs/FigureS7.pdf',
    figureS8_artwork = 'inputs/FigureS8.pdf'
  )
)
