# Script to construct PCA and RLA plot of sample quenching data
# Author: Benjamin R. Gordon
# Date: 2019-04-12

library(tidyverse)

# load data
load("./data/snapdata.rda")

# PCA
set.seed(1978)
pca <- stats::prcomp(snapdata[-1:-2], 
                     scale = FALSE, 
                     center = TRUE)

# Plotting
exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(snapdata[1:2], pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")

pca <- ggplot(data = filter(scores, class == "methanol" | class == "nitrogen"),
              aes(x = PC1,
                  y = PC2,
                  color = class,
                  fill = class,
                  shape = class)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores$PC1)),
                                        height = 0.01 * diff(range(scores$PC2))
             )) +
  labs(x = x_lab, y = y_lab) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = 
                       grDevices::adjustcolor(phdhelpr::qual_colours[c(6, 2)],
                                              alpha.f = 0.9)) +
  scale_fill_manual(values = 
                      grDevices::adjustcolor(phdhelpr::qual_colours[c(6, 2)],
                                             alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none")

rla <- 
  ggplot(data = reshape2::melt(filter(snapdata, 
                                      class == "methanol" | class == "nitrogen")),
         aes(x = sample_id, 
             y = log(value, 2) - median(log(value, 2)),
             fill = class)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = phdhelpr::qual_colours[c(6, 2)]) +
  scale_x_discrete(name = "sample") +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = "none")

# create grobs and set corner labels
plot_a <-
  gridExtra::arrangeGrob(pca, 
                         top = grid::textGrob("a",
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

plot_b <-
  gridExtra::arrangeGrob(rla,
                         top = grid::textGrob("b",
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

# print plot
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1))

# save plot
grDevices::pdf("./figs/pca_rla_plot.pdf",
               width = 10,
               height = 5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1))
grDevices::dev.off()