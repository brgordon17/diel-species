# Script to construct figures in chapter 6 (heron 24 h)
# Author: Benjamin R. Gordon
# Date: 2019-04-12

# VIP plot of diurnal differences ----------------------------------------------
library(tidyverse)
library(reshape2)

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
                        colour = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(name = "Time (hh:mm)", 
                     values = gordon01::qual_colours[c(1:3, 6)]) +
  scale_shape_manual(name = "Time (hh:mm)", 
                     values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank())
aspe_plot

aequ_plot <- ggplot(aequ_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(values = gordon01::qual_colours[c(1:3, 6)]) +
  scale_shape_manual(values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")
aequ_plot

digi_plot <- ggplot(digi_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(values = gordon01::qual_colours[c(1:3, 6)]) +
  scale_shape_manual(values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")
digi_plot

cyli_plot <- ggplot(cyli_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(values = gordon01::qual_colours[c(1:3, 6)]) +
  scale_shape_manual(values = c(21:24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
        legend.position = "none")
cyli_plot

dami_plot <- ggplot(dami_impvars, 
                    aes(x = importance, 
                        y = mz,
                        shape = max_class,
                        colour = max_class)) +
  geom_point(size = 3, stroke = 1) +
  scale_x_continuous(name = NULL) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(values = gordon01::qual_colours[c(1, 3, 6)]) +
  scale_shape_manual(values = c(21, 23, 24)) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90",
                                          size = 0.4),
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
