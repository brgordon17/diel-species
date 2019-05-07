# Script to conduct stats anlaysis of VIP features.
# Author: Benjamin R Gordon.
# Date: 2019-04-27

library(tidyverse)

# Load and prep data
load("./data/mzdata.rda")
matches <- read_rds("./dev/impvar_matches.rds")

# one way anovas (kruskal-wallis test) for each species ------------------------
# for A. aspera
features_aspe <- 
  mzdata %>%
  filter(taxon != "PBQC" & taxon == "A_aspe") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

# test for normality (normal if p > 0.05)
apply(features_aspe[5:ncol(features_aspe)], 2, shapiro.test)

# some of the variables are non-normal so perform a Kruskal-Wallis test
apply(features_aspe[5:ncol(features_aspe)], 2, kruskal.test, 
      g = features_aspe$time_fac)

# For M. digitata
features_digi <- 
  mzdata %>%
  filter(taxon != "PBQC" & taxon == "M_digi") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

# test for normality (normal if p > 0.05)
apply(features_digi[5:ncol(features_digi)], 2, shapiro.test)

# all of the variables are non-normal so perform a Kruskal-Wallis test
# on all variables within each species
apply(features_digi[5:ncol(features_digi)], 2, kruskal.test, 
      g = features_digi$time_fac)
