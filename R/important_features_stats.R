# Script to conduct stats anlaysis of VIP features.
# Author: Benjamin R Gordon.
# Date: 2019-04-27

library(tidyverse)

# Load and prep data
load("./data/mzdata.rda")
matches <- read_rds("./dev/impvar_matches.rds")

# one way anovas for A. aspera -------------------------------------------------
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

formula <- 
  as.formula(paste0("cbind(", paste(names(features_aspe)[5:ncol(features_aspe)], 
                                    collapse = ","), ") ~ time_fac"))
summary(aov(formula, data = features_aspe))

# one way anovas for M. aequ ---------------------------------------------------
features_aequ <- 
  mzdata %>%
  filter(taxon != "PBQC" & taxon == "M_aequ") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

formula <- 
  as.formula(paste0("cbind(", paste(names(features_aequ)[5:ncol(features_aequ)], 
                                    collapse = ","), ") ~ time_fac"))
summary(aov(formula, data = features_aequ))

# one way anovas for M. digitata -----------------------------------------------
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

formula <- 
  as.formula(paste0("cbind(", paste(names(features_digi)[5:ncol(features_digi)], 
                                    collapse = ","), ") ~ time_fac"))
summary(aov(formula, data = features_digi))

# one way anovas for P. cylindrica ---------------------------------------------
features_cyli <- 
  mzdata %>%
  filter(taxon != "PBQC" & taxon == "P_cyli") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

formula <- 
  as.formula(paste0("cbind(", paste(names(features_cyli)[5:ncol(features_cyli)], 
                                    collapse = ","), ") ~ time_fac"))
summary(aov(formula, data = features_cyli))

# one way anovas for P. damicornis ---------------------------------------------
features_dami <- 
  mzdata %>%
  filter(taxon != "PBQC" & taxon == "P_dami") %>%
  droplevels() %>%
  select(1:4,
         mz_784.53491,
         mz_277.18008,
         mz_1019.77533,
         mz_347.25863,
         mz_815.55253,
         mz_135.04778)

formula <- 
  as.formula(paste0("cbind(", paste(names(features_dami)[5:ncol(features_dami)], 
                                    collapse = ","), ") ~ time_fac"))
summary(aov(formula, data = features_dami))

# Kruskal wallis tests ---------------------------------------------------------
apply(features_cyli[5:ncol(features_cyli)], 2, kruskal.test, 
      g = features_cyli$time_fac)
