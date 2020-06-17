# Script to construct PCA and RLA plot of taxonomic data
# Author: Benjamin R. Gordon
# Date: 2019-04-12

library(tidyverse)

# Load data
load("./data/mzdata.rda")

# PCA plot 
set.seed(1978)
pca <- stats::prcomp(mzdata[-1:-9], 
                     scale = FALSE, 
                     center = TRUE)

exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(mzdata[1:9], pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")

pcaplot <- ggplot(data = filter(scores, taxon != "PBQC"),
                  aes(x = PC1,
                      y = PC2,
                      color = taxon,
                      fill = taxon,
                      shape = taxon)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores$PC1)),
                                        height = 0.01 * diff(range(scores$PC2))
             )) +
  labs(x = x_lab, y = y_lab) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 11),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.position = "top")

# RLA plot
rladata <- reshape2::melt(filter(mzdata[-5:-9], taxon != "PBQC"))
rladata$taxon <- factor(rladata$taxon, 
                        levels = c("M_aequ", "A_aspe", "P_cyli", "P_dami", 
                                   "M_digi"))

rlaplot <- ggplot(data = rladata,
                  aes(x = sample_id, 
                      y = log(value, 2) - median(log(value, 2)),
                      fill = taxon)) +
  geom_boxplot(alpha = 0.7,
               outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = phdhelpr::qual_colours[c(2, 1, 6, 7, 3)]) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 11),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "top")

# create grobs and set corner labels
plot_a <-
  gridExtra::arrangeGrob(pcaplot, 
                         top = grid::textGrob("a",
                                              x = grid::unit(0.017, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 16)
                         ))

plot_b <-
  gridExtra::arrangeGrob(rlaplot,
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
grDevices::pdf("./figs/taxonomy_pca_rla_plot.pdf",
               width = 10.5,
               height = 5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1))
grDevices::dev.off()
