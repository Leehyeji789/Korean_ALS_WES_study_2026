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
  strip.text = element_text(size = 10, family = 'Arial', color = 'black'),
  
)
