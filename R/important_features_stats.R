# Script to conduct stats anlaysis of VIP features.
# Author: Benjamin R Gordon.
# Date: 2019-04-27

library(tidyverse)

# Load and prep data
load("./data/mzdata.rda")
matches <- read_rds("./dev/impvar_matches.rds")

features <- 
  mzdata %>%
  filter(taxon != "PBQC") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

# test for normality (normal if p > 0.05)
apply(features[5:ncol(features)], 2, shapiro.test)

# all of the variables are non-normal so perform a Kruskal-Wallis test
# on all variables within each species
apply(features[5:ncol(features)], 2, kruskal.test, g = features$time_fac)
