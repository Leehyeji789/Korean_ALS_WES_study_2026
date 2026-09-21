#!/usr/bin/env Rscript
# Post-QC ancestry PCA; retain the original Somalier panel definitions.
source('R/common.R'); source('R/plot_style.R')
suppressPackageStartupMessages(library(dplyr))
cfg <- read_config()
ancestry <- as.data.frame(read_input(cfg, 'ancestry', c('#sample_id', 'predicted_ancestry', 'PC1', 'PC2', 'PC3')))
names(ancestry)[names(ancestry) == '#sample_id'] <- 'sample_id'
ancestry_labels <- as.data.frame(read_input(cfg, 'ancestry_labels', '#sample'))
names(ancestry_labels)[names(ancestry_labels) == '#sample'] <- 'sample_id'
ancestry_labels <- unique(ancestry_labels['sample_id'])
theme_als_v1 <- theme(
  plot.title = element_text(size = 13, family = 'Arial', color = 'black'),
  text = element_text(size = 10, family = 'Arial', color = 'black'),
  
  axis.title = element_text(size = 11, family = 'Arial', color = 'black'),
  axis.text = element_text(size = 10, family = 'Arial', color = 'black'),
  # axis.line = element_line(color = 'black', linewidth = .5),
  
  legend.title = element_text(size = 11, family = 'Arial', color = 'black'),
  legend.text = element_text(size = 10, family = 'Arial', color = 'black'),
  legend.background = element_rect(color = 'black', linewidth = .5),
  legend.box = c('vertical'), # arrangement of multiple legends
  legend.margin = margin(3, 3, 3, 3),
  
  panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "grey90"),
  # panel.border = element_blank(),
  
  # Facet
  strip.background = element_rect(color = 'black', linewidth = 1, fill = 'grey88'),
  strip.text = element_text(size = 10, family = 'Arial', color = 'black')
  
)
ancestry_colors <- c(
  'EUR' = '#984EA3',
  'AMR' = '#FF7F00',
  'EAS' = '#4DAF4A',
  'SAS' = '#E41A1C',
  'AFR' = '#377EB8'
)

background_df <- ancestry %>%
  inner_join(ancestry_labels, by = "sample_id") %>%
  mutate(type = '1000 Genome')

query_df <- ancestry %>%
  anti_join(ancestry_labels, by = "sample_id") %>%
  mutate(type = 'ALS')

combined_df <- bind_rows(background_df, query_df) %>%
  mutate(ancestry = factor(predicted_ancestry, levels = c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')),
         across(c(PC1, PC2, PC3), as.numeric))

# PC1 vs PC2
p1 <- ggplot(combined_df, aes(x = PC1, y = PC2, color = ancestry, shape = type)) +
  geom_point(aes(alpha = type), size = 3.5) +
  scale_color_manual(values = ancestry_colors, name = 'Predicted Ancestry', na.value = 'grey50') +
  scale_shape_manual(values = c('1000 Genome' = 16, 'ALS' = 15), name = 'Data Type') +
  scale_alpha_manual(values = c('1000 Genome' = 0.2, 'ALS' = 0.8)) +
  labs(title = 'Ancestry PCA (PC1 vs PC2)', x = 'PC 1', y = 'PC 2') +
  theme_als_v1 +
  guides(alpha = 'none')

# PC1 vs PC3
p2 <- ggplot(combined_df, aes(x = PC1, y = PC3, color = ancestry, shape = type)) +
  geom_point(aes(alpha = type), size = 3.5) +
  scale_color_manual(values = ancestry_colors, name = 'Predicted Ancestry', na.value = 'grey50') +
  scale_shape_manual(values = c('1000 Genome' = 16, 'ALS' = 15), name = 'Data Type') +
  scale_alpha_manual(values = c('1000 Genome' = 0.2, 'ALS' = 0.8)) +
  labs(title = 'Ancestry PCA (PC1 vs PC3)', x = 'PC 1', y = 'PC 3') +
  theme_als_v1 +
  guides(alpha = 'none')


out <- new_stage(cfg, 'FigureS2')
p <- ggarrange(p1, p2, ncol = 2, labels = c('A', 'B'),
  font.label = list(size = 16, face = 'bold', family = 'Arial'),
  common.legend = TRUE, legend = 'bottom')
save_figure(p, file.path(out, 'FigureS2'), 10, 5)
record_session(out)
