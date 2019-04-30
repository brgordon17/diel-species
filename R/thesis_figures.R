# Script to construct figures in chapter 6 (heron 24 h)
# Author: Benjamin R. Gordon
# Date: 2019-04-12

# PCA and RLA plot of quenching variation --------------------------------------
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
                       grDevices::adjustcolor(gordon01::qual_colours[c(6, 2)],
                                              alpha.f = 0.9)) +
  scale_fill_manual(values = 
                      grDevices::adjustcolor(gordon01::qual_colours[c(6, 2)],
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
  scale_fill_manual(values = gordon01::qual_colours[c(6, 2)]) +
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

# plot of quneched sample feature intensities ----------------------------------
library(tidyverse)

# Load data
load("./data/snapdata.rda")
mzrf <- read_rds("./dev/quenching_model.rds")

# identify the 9 most important variables
mzrf_impvars <- varImp(mzrf, scale = FALSE)
mzrf_impvars <- mzrf_impvars$importance
mzrf_impvars <- cbind(importance = apply(mzrf_impvars, 1, max), mzrf_impvars)
mzrf_impvars <- mzrf_impvars[order(-mzrf_impvars$importance), ,drop = FALSE]
mzrf_impvars <- as_tibble(mzrf_impvars[1:12, ], rownames = "mz")

# get abundances
sum_snapdata <- 
  snapdata %>% 
  filter(class == "methanol" | class == "nitrogen") %>%
  select(1:2, c(mzrf_impvars$mz)) %>%
  reshape2::melt(.) %>%
  as_tibble()
sum_snapdata$variable <- str_replace_all(sum_snapdata$variable, "mz_", "")

# plot
ggplot(sum_snapdata, aes(x = class,
                         y = value/10^4,
                         colour = class,
                         fill = class)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0.1, linetype = 3) +
  facet_wrap(vars(variable), scales = "free") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = expression(Intensity%*%10^4)) +
  scale_colour_manual(values = gordon01::qual_colours[c(6, 2)]) +
  theme(axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
  )

# save plot
ggsave("./figs/quenching_intensities.pdf",
       width = 10,
       height = 6.5,
       units = "in")

# taxonomy PCA and RLA plot ----------------------------------------------------
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 11),
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
  scale_fill_manual(values = gordon01::qual_colours[c(2, 1, 6, 7, 3)]) +
  scale_x_discrete(name = "Sample") +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 11),
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
               width = 10,
               height = 5,
               useDingbats = FALSE)
gridExtra::grid.arrange(gridExtra::arrangeGrob(plot_a,
                                               plot_b,
                                               nrow = 1))
grDevices::dev.off()

# PCA composite plot of diurnal variation --------------------------------------
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
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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

# PCA composite plot of diurnal variation (LOG data) ---------------------------
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
aspe_pca <- stats::prcomp(log(aspe_data[-1:-9], 2), 
                          scale = FALSE, 
                          center = TRUE)
aspe_scores <- as_tibble(data.frame(aspe_data[1:9], aspe_pca$x/10^6))

set.seed(1978)
aequ_pca <- stats::prcomp(log(aequ_data[-1:-9], 2),
                          scale = FALSE, 
                          center = TRUE)
aequ_scores <- as_tibble(data.frame(aequ_data[1:9], aequ_pca$x/10^6))

set.seed(1978)
digi_pca <- stats::prcomp(log(digi_data[-1:-9], 2),
                          scale = FALSE, 
                          center = TRUE)
digi_scores <- as_tibble(data.frame(digi_data[1:9], digi_pca$x/10^6))

set.seed(1978)
cyli_pca <- stats::prcomp(log(cyli_data[-1:-9], 2),
                          scale = FALSE, 
                          center = TRUE)
cyli_scores <- as_tibble(data.frame(cyli_data[1:9], cyli_pca$x/10^6))

set.seed(1978)
dami_pca <- stats::prcomp(log(dami_data[-1:-9], 2),
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
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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
grDevices::pdf("./figs/pca_plot_diurnal_logged.pdf",
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

# VIP plot of taxonomic differences --------------------------------------------
library(tidyverse)
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
                     values = gordon01::qual_colours[c(1:3, 6, 7)]) +
  scale_shape_manual(name = "Species", 
                     values = c(21:25)) +
  scale_fill_manual(name = "Species",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
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

# plot of taxonomic feature intensities ----------------------------------------
library(tidyverse)
library(caret)

# Load data
load("./data/mzdata.rda")
mzrf <- read_rds("./dev/mzrf_taxon.rds")

# identify the 12 most important variables
mzrf_impvars <- varImp(mzrf, scale = FALSE)
mzrf_impvars <- mzrf_impvars$importance
mzrf_impvars <- cbind(importance = apply(mzrf_impvars, 1, max), mzrf_impvars)
mzrf_impvars <- mzrf_impvars[order(-mzrf_impvars$importance), ,drop = FALSE]
mzrf_impvars <- as_tibble(mzrf_impvars[1:12, ], rownames = "mz")

# get abundances
sum_mzdata <- 
  mzdata %>% 
  filter(taxon != "PBQC") %>%
  select(1, 3, c(mzrf_impvars$mz)) %>%
  reshape2::melt(.) %>%
  as_tibble()
sum_mzdata$variable <- str_replace_all(sum_mzdata$variable, "mz_", "")
sum_mzdata <- droplevels(sum_mzdata)

# plot
ggplot(sum_mzdata, 
       aes(x = taxon,
           y = value/10^4,
           colour = taxon,
           fill = taxon)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0.1, linetype = 3) +
  facet_wrap(vars(variable), scales = "free") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = expression(Intensity%*%10^4)) +
  scale_colour_manual(values = gordon01::qual_colours[c(1:3, 6, 7)]) +
  scale_fill_manual(values = gordon01::qual_colours[c(1:3, 6, 7)]) +
  theme(axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
  )

# save plot
ggsave("./figs/taxon_intensities.pdf",
       width = 10,
       height = 6.5,
       units = "in")

# VIP plot of diurnal differences ----------------------------------------------
library(tidyverse)
library(reshape2)
library(caret)

# Load models
mzrf_aspe <- read_rds("./dev/mzrf_aspe.rds")
mzrf_aequ <- read_rds("./dev/mzrf_aequ.rds")
mzrf_digi <- read_rds("./dev/mzrf_digi.rds")
mzrf_cyli <- read_rds("./dev/mzrf_cyli.rds")
mzrf_dami <- read_rds("./dev/mzrf_dami.rds")

# Identify important variables
# for aspe
aspe_impvars <- varImp(mzrf_aspe, scale = FALSE)
aspe_impvars <- aspe_impvars$importance
aspe_impvars <- 
  bind_cols(max_class = factor(names(aspe_impvars)[apply(aspe_impvars, 1, which.max)]),
            importance = apply(aspe_impvars, 1, max),
            mz = str_replace(rownames(aspe_impvars), "mz_", ""),
            aspe_impvars)
aspe_impvars <- aspe_impvars[order(-aspe_impvars$importance), ,drop = FALSE]
aspe_impvars <- aspe_impvars[1:20, ]
aspe_impvars$mz <- round(as.numeric(aspe_impvars$mz), 4)
aspe_impvars$mz <- factor(aspe_impvars$mz, levels = rev(aspe_impvars$mz))
aspe_impvars <- select(aspe_impvars, mz, everything())

# for aequ
aequ_impvars <- varImp(mzrf_aequ, scale = FALSE)
aequ_impvars <- aequ_impvars$importance
aequ_impvars <- 
  bind_cols(max_class = factor(names(aequ_impvars)[apply(aequ_impvars, 1, which.max)]),
            importance = apply(aequ_impvars, 1, max),
            mz = str_replace(rownames(aequ_impvars), "mz_", ""),
            aequ_impvars)
aequ_impvars <- aequ_impvars[order(-aequ_impvars$importance), ,drop = FALSE]
aequ_impvars <- aequ_impvars[1:20, ]
aequ_impvars$mz <- round(as.numeric(aequ_impvars$mz), 4)
aequ_impvars$mz <- factor(aequ_impvars$mz, levels = rev(aequ_impvars$mz))
aequ_impvars <- select(aequ_impvars, mz, everything())

# for digi
digi_impvars <- varImp(mzrf_digi, scale = FALSE)
digi_impvars <- digi_impvars$importance
digi_impvars <- 
  bind_cols(max_class = factor(names(digi_impvars)[apply(digi_impvars, 1, which.max)]),
            importance = apply(digi_impvars, 1, max),
            mz = str_replace(rownames(digi_impvars), "mz_", ""),
            digi_impvars)
digi_impvars <- digi_impvars[order(-digi_impvars$importance), ,drop = FALSE]
digi_impvars <- digi_impvars[1:20, ]
digi_impvars$mz <- round(as.numeric(digi_impvars$mz), 4)
digi_impvars$mz <- factor(digi_impvars$mz, levels = rev(digi_impvars$mz))
digi_impvars <- select(digi_impvars, mz, everything())

# for cyli
cyli_impvars <- varImp(mzrf_cyli, scale = FALSE)
cyli_impvars <- cyli_impvars$importance
cyli_impvars <- 
  bind_cols(max_class = factor(names(cyli_impvars)[apply(cyli_impvars, 1, which.max)]),
            importance = apply(cyli_impvars, 1, max),
            mz = str_replace(rownames(cyli_impvars), "mz_", ""),
            cyli_impvars)
cyli_impvars <- cyli_impvars[order(-cyli_impvars$importance), ,drop = FALSE]
cyli_impvars <- cyli_impvars[1:20, ]
cyli_impvars$mz <- round(as.numeric(cyli_impvars$mz), 4)
cyli_impvars$mz <- factor(cyli_impvars$mz, levels = rev(cyli_impvars$mz))
cyli_impvars <- select(cyli_impvars, mz, everything())

# for dami
dami_impvars <- varImp(mzrf_dami, scale = FALSE)
dami_impvars <- dami_impvars$importance
dami_impvars <- 
  bind_cols(max_class = factor(names(dami_impvars)[apply(dami_impvars, 1, which.max)]),
            importance = apply(dami_impvars, 1, max),
            mz = str_replace(rownames(dami_impvars), "mz_", ""),
            dami_impvars)
dami_impvars <- dami_impvars[order(-dami_impvars$importance), ,drop = FALSE]
dami_impvars <- dami_impvars[1:20, ]
dami_impvars$mz <- round(as.numeric(dami_impvars$mz), 4)
dami_impvars$mz <- factor(dami_impvars$mz, levels = rev(dami_impvars$mz))
dami_impvars <- select(dami_impvars, mz, everything())

# plotting
aspe_plot <- ggplot(aspe_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class,
                        fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank())
aspe_plot

aequ_plot <- ggplot(aequ_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class,
                        fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
aequ_plot

digi_plot <- ggplot(digi_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class,
                        fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
digi_plot

cyli_plot <- ggplot(cyli_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class,
                        fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
cyli_plot

dami_plot <- ggplot(dami_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class,
                        fill = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)",
                     values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(name = "Time (hh:mm)",
                    values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
dami_plot

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
xgrob <- grid::textGrob("Mean Decrease in Accuracy",
                        x = grid::unit(0.5, "npc"),
                        y = grid::unit(0.5, "npc"),
                        just = c("centre"),
                        gp = grid::gpar(fontsize = 12)
                        )

ygrob <- grid::textGrob("Important Spectral Features (Da)",
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
grDevices::pdf("./figs/vip_plot_diurnal.pdf",
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

# plot of important diurnal feature intensities --------------------------------
library(tidyverse)
library(caret)

# Load data
load("./data/mzdata.rda")

# feature abundances
sum_mzdata <- 
  mzdata %>% 
  filter(taxon != "PBQC") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778) %>%
  reshape2::melt(.) %>%
  group_by(variable, taxon, time_fac) %>%
  summarise(mean = mean(value/10^4), 
            sd = sd(value/10^4), 
            se = sd(value/10^4)/sqrt(n())
            ) %>%
  ungroup()
sum_mzdata$variable <- as.factor(str_replace_all(sum_mzdata$variable, "mz_", ""))

# named vector and scaling values
names <- c(
  "784.53491" = "m/z 784.53491",
  "277.18008" = "m/z 277.18008",
  "1019.77533" = "m/z 1019.77533",
  "347.25863" = "m/z 347.25863",
  "815.55253" = "m/z 815.55253",
  "135.04778" = "m/z 135.04778")

# create plot
ggplot(sum_mzdata, aes(x = time_fac,
                       y = mean,
                       fill = taxon,
                       colour = taxon,
                       shape = taxon,
                       group = taxon)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = .1,
                show.legend = FALSE) +
  geom_path(size = 0.8, show.legend = FALSE) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.1, linetype = 3) +
  facet_wrap(vars(variable), labeller = as_labeller(names), scales = "free") +
  scale_x_discrete(name = "Time (hh:mm)") +
  scale_y_continuous(name = Mean~Intensity%*%10^4) +
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(values = c(21:25)) +
  theme(axis.text = element_text(size = 12, colour = "grey65"),
        axis.title = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        strip.text = element_text(size = 16),
        strip.background = element_blank()
  )

ggsave("./figs/diurnal_intensities.pdf")
