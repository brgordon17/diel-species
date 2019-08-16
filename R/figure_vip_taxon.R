# Script to construct vip plot of taxon data
# Author: Benjamin R. Gordon
# Date: 2019-04-12library(tidyverse)

library(reshape2)

# Load model
mzrf_taxon <- read_rds("./dev/mzrf_taxon.rds")

# Identify important variables
taxon_impvars <- varImp(mzrf_taxon, scale = FALSE)
taxon_impvars <- taxon_impvars$importance
taxon_impvars <- 
  bind_cols(max_class = factor(names(taxon_impvars)[apply(taxon_impvars, 1, which.max)]),
            importance = apply(taxon_impvars, 1, max),
            mz = str_replace(rownames(taxon_impvars), "mz_", ""),
            taxon_impvars)
taxon_impvars <- taxon_impvars[order(-taxon_impvars$importance), ,drop = FALSE]
taxon_impvars <- taxon_impvars[1:20, ]
taxon_impvars$mz <- round(as.numeric(taxon_impvars$mz), 4)
taxon_impvars$mz <- factor(taxon_impvars$mz, levels = rev(taxon_impvars$mz))
taxon_impvars <- select(taxon_impvars, mz, everything())

# plotting
ggplot(taxon_impvars,
       aes(x = importance, 
           y = mz,
           shape = max_class,
           colour = max_class,
           fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = "Mean Decrease in Accuracy") +
  scale_y_discrete(name = "Important Spectral Features (Da)") +
  scale_color_manual(name = "Species", 
                     values = phdhelpr::qual_colours[c(1:3, 6, 7)]) +
  scale_shape_manual(name = "Species", 
                     values = c(21:25)) +
  scale_fill_manual(name = "Species",
                    values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.3)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank())

# save plot
ggsave("./figs/taxon_vip_plot.pdf",
       width = 10,
       height = 6.5,
       units = "in")
