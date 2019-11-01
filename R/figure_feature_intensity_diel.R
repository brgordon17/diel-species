# Script to construct feature intensity plot of diel data
# Author: Benjamin R. Gordon
# Date: 2019-04-12library(tidyverse)library(tidyverse)library(tidyverse)

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
  "784.53491" = "arachidonolthio-PC",
  "277.18008" = "methyl montiporate B",
  "1019.77533" = "lyso-PAF C18",
  "347.25863" = "10-hydroxydocosa-\npentaenoic acid",
  "815.55253" = "pyrophaeophytin-a",
  "135.04778" = "dimethylsulfonio-\npropionate")

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
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(values = c(21:25)) +
  theme(axis.text = element_text(size = 12, colour = "grey65"),
        axis.title = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_blank()
  )

ggsave("./figs/diurnal_intensities.pdf")

# plot of mdigitata m/z 135. Not printed
ggplot(filter(sum_mzdata, taxon == "M_digi" & variable == "135.04778"),
       aes(x = time_fac,
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
  scale_x_discrete(name = "Time (hh:mm)") +
  scale_y_continuous(name = Mean~Intensity%*%10^4) +
  scale_color_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(phdhelpr::qual_colours[c(1:3, 6, 7)],
                                                    alpha.f = 0.5)) +
  scale_shape_manual(values = c(21:25)) +
  theme(axis.text = element_text(size = 12, colour = "grey65"),
        axis.title = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_blank()
  )
