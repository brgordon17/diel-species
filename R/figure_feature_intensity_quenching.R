# Script to construct feature intensities of sample quenching data
# Author: Benjamin R. Gordon
# Date: 2019-04-12

library(tidyverse)
library(caret)

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
  scale_colour_manual(values = phdhelpr::qual_colours[c(6, 2)]) +
  theme(axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
  )

# save plot
ggsave("./figs/quenching_intensities.pdf",
       width = 10,
       height = 6.5,
       units = "in")
