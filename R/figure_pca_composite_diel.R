# Script to construct PCA composite of diel variation
# Author: Benjamin R. Gordon
# Date: 2019-04-12

library(tidyverse)
library(reshape2)

# Load and prep data
load("./data/mzdata.rda")
aspe_data <- filter(mzdata, taxon == "A_aspe")
aequ_data <- filter(mzdata, taxon == "M_aequ")
digi_data <- filter(mzdata, taxon == "M_digi")
cyli_data <- filter(mzdata, taxon == "P_cyli")
dami_data <- filter(mzdata, taxon == "P_dami")

# pca's
set.seed(1978)
aspe_pca <- stats::prcomp(aspe_data[-1:-9], 
                          scale = FALSE, 
                          center = TRUE)
aspe_scores <- as_tibble(data.frame(aspe_data[1:9], aspe_pca$x/10^6))

set.seed(1978)
aequ_pca <- stats::prcomp(aequ_data[-1:-9], 
                          scale = FALSE, 
                          center = TRUE)
aequ_scores <- as_tibble(data.frame(aequ_data[1:9], aequ_pca$x/10^6))

set.seed(1978)
digi_pca <- stats::prcomp(digi_data[-1:-9], 
                          scale = FALSE, 
                          center = TRUE)
digi_scores <- as_tibble(data.frame(digi_data[1:9], digi_pca$x/10^6))

set.seed(1978)
cyli_pca <- stats::prcomp(cyli_data[-1:-9], 
                          scale = FALSE, 
                          center = TRUE)
cyli_scores <- as_tibble(data.frame(cyli_data[1:9], cyli_pca$x/10^6))

set.seed(1978)
dami_pca <- stats::prcomp(dami_data[-1:-9], 
                          scale = FALSE, 
                          center = TRUE)
dami_scores <- as_tibble(data.frame(dami_data[1:9], dami_pca$x/10^6))

# plots
aspe_plot <- ggplot(data = aspe_scores,
                    aes(x = PC1,
                        y = PC2,
                        color = time_fac,
                        fill = time_fac,
                        shape = time_fac)) +
  geom_point(size = 2.5, stroke = 0.7) +
  scale_shape_manual(name = "Time (hh:mm)",
                     values = c(21:25)) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(axis.ticks = element_blank(),panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.position = "right")

aequ_plot <- ggplot(data = aequ_scores,
                    aes(x = PC1,
                        y = PC2,
                        color = time_fac,
                        fill = time_fac,
                        shape = time_fac)) +
  geom_point(size = 2.5, stroke = 0.7) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")

digi_plot <- ggplot(data = digi_scores,
                    aes(x = PC1,
                        y = PC2,
                        color = time_fac,
                        fill = time_fac,
                        shape = time_fac)) +
  geom_point(size = 2.5, stroke = 0.7) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")

cyli_plot <- ggplot(data = cyli_scores,
                    aes(x = PC1,
                        y = PC2,
                        color = time_fac,
                        fill = time_fac,
                        shape = time_fac)) +
  geom_point(size = 2.5, stroke = 0.7) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")

dami_plot <- ggplot(data = dami_scores,
                    aes(x = PC1,
                        y = PC2,
                        color = time_fac,
                        fill = time_fac,
                        shape = time_fac)) +
  geom_point(size = 2.5, stroke = 0.7) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")

# Construct composite plot
# Extract Legend 
g_tab <- ggplot_gtable(ggplot_build(aspe_plot))
leg <- which(sapply(g_tab$grobs, function(x) x$name) == "guide-box")
legend <- g_tab$grobs[[leg]]
rm(g_tab, leg)

# Set corner labels
plot_aspe <-
  gridExtra::arrangeGrob(aspe_plot +
                           theme(legend.position = "none"),
                         top = grid::textGrob("A. aspera",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)
                         ))
plot_aequ <-
  gridExtra::arrangeGrob(aequ_plot,
                         top = grid::textGrob("M. aequituberculata",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)
                         ))
plot_digi <-
  gridExtra::arrangeGrob(digi_plot,
                         top = grid::textGrob("M. digitata",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)
                         ))

plot_cyli <-
  gridExtra::arrangeGrob(cyli_plot,
                         top = grid::textGrob("P. cylindrica",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)
                         ))

plot_dami <-
  gridExtra::arrangeGrob(dami_plot,
                         top = grid::textGrob("P. damicornis",
                                              x = grid::unit(0.5, "npc"),
                                              y = grid::unit(0.5, "npc"),
                                              gp = grid::gpar(fontsize = 12)
                         ))

# create x- and y-axis text grobs
xgrob <- grid::textGrob(expression(PC1%*%10^4),
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        gp = grid::gpar(fontsize = 12)
)

ygrob <- grid::textGrob(expression(PC2%*%10^4),
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        rot = 90,
                        gp = grid::gpar(fontsize = 12)
)

gridExtra::grid.arrange(plot_aspe,
                        plot_aequ,
                        plot_cyli,
                        plot_digi,
                        plot_dami,
                        legend,
                        nrow = 2,
                        left = ygrob,
                        bottom = xgrob)

# save
grDevices::pdf("./figs/pca_plot_diurnal.pdf",
               width = 10,
               height = 8,
               useDingbats = FALSE)
gridExtra::grid.arrange(plot_aspe,
                        plot_aequ,
                        plot_cyli,
                        plot_digi,
                        plot_dami,
                        legend,
                        nrow = 2,
                        left = ygrob,
                        bottom = xgrob)
grDevices::dev.off()
