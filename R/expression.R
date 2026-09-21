supercluster_order <- c(
  'Upper-layer intratelencephalic', 'Deep-layer intratelencephalic',
  'Deep-layer near-projecting', 'Deep-layer corticothalamic and 6b',
  'MGE interneuron', 'CGE interneuron', 'LAMP5-LHX6 and Chandelier', 'Miscellaneous',
  'Hippocampal CA1-3', 'Hippocampal CA4', 'Hippocampal dentate gyrus',
  'Amygdala excitatory', 'Medium spiny neuron', 'Eccentric medium spiny neuron',
  'Splatter', 'Mammillary body', 'Thalamic excitatory', 'Midbrain-derived inhibitory',
  'Upper rhombic lip', 'Cerebellar inhibitory', 'Lower rhombic lip',
  'Oligodendrocyte', 'Committed oligodendrocyte precursor', 'Oligodendrocyte precursor',
  'Astrocyte', 'Ependymal', 'Microglia', 'Vascular', 'Bergmann glia',
  'Fibroblast', 'Choroid plexus')

expression_gene_map <- function(cfg) {
  gtf <- rtracklayer::import(cfg$inputs$expression_gtf)
  x <- data.table(gene_id = gtf$gene_id, gene_name = gtf$gene_name,
                  gene_type = gtf$gene_type)
  unique(x[!is.na(gene_id) & !is.na(gene_name) & !is.na(gene_type)], by = 'gene_id')
}
