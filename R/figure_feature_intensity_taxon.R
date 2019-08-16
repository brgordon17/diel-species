# Script to construct feature intensities of taxon data
# Author: Benjamin R. Gordon
# Date: 2019-04-12library(tidyverse)

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
  scale_colour_manual(values = phdhelpr::qual_colours[c(1:3, 6, 7)]) +
  scale_fill_manual(values = phdhelpr::qual_colours[c(1:3, 6, 7)]) +
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
