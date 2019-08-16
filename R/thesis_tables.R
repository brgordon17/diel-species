# ---------------- DEPRECATED --------------------------------------------------
# Use individual table_*.R scripts

# Script for producing tables used in Chapter 6 (Heron 24 h)
# Author: Benjamin R. Gordon
# Date: 2019-04-12

# diurnal model results --------------------------------------------------------
# Table built in excel from results below
library(tidyverse)
library(caret)

# load models
mzrf_aspe <- read_rds("./dev/mzrf_aspe.rds")
mzrf_aequ <- read_rds("./dev/mzrf_aequ.rds")
mzrf_digi <- read_rds("./dev/mzrf_digi.rds")
mzrf_cyli <- read_rds("./dev/mzrf_cyli.rds")
mzrf_dami <- read_rds("./dev/mzrf_dami.rds")

# conmats and stats
conmat_aspe <- confusionMatrix(mzrf_aspe$pred$pred[mzrf_aspe$pred$mtry == 
                                                     mzrf_aspe$bestTune$mtry],
                               mzrf_aspe$pred$obs[mzrf_aspe$pred$mtry ==
                                                    mzrf_aspe$bestTune$mtry])

conmat_aequ <- confusionMatrix(mzrf_aequ$pred$pred[mzrf_aequ$pred$mtry == 
                                                     mzrf_aequ$bestTune$mtry],
                               mzrf_aequ$pred$obs[mzrf_aequ$pred$mtry ==
                                                    mzrf_aequ$bestTune$mtry])

conmat_digi <- confusionMatrix(mzrf_digi$pred$pred[mzrf_digi$pred$mtry == 
                                                     mzrf_digi$bestTune$mtry],
                               mzrf_digi$pred$obs[mzrf_digi$pred$mtry ==
                                                    mzrf_digi$bestTune$mtry])

conmat_cyli <- confusionMatrix(mzrf_cyli$pred$pred[mzrf_cyli$pred$mtry == 
                                                     mzrf_cyli$bestTune$mtry],
                               mzrf_cyli$pred$obs[mzrf_cyli$pred$mtry ==
                                                    mzrf_cyli$bestTune$mtry])

conmat_dami <- confusionMatrix(mzrf_dami$pred$pred[mzrf_dami$pred$mtry == 
                                                     mzrf_dami$bestTune$mtry],
                               mzrf_dami$pred$obs[mzrf_dami$pred$mtry ==
                                                    mzrf_dami$bestTune$mtry])

# annotated diurnal features table ---------------------------------------------
library(tidyverse)

# Load and prep data
matches <- read_rds("./dev/impvar_matches.rds")

# remove generic referencing and replace commas
matches <-
  matches %>%
  mutate(endnote_ref = stringr::str_replace(endnote_ref, ",", ";")) %>%
  select(-ref)


# Save csv 
readr::write_csv(matches, "./tables/annotated_diurnal_features.txt")

# NOTE: Fragment ions must be removed after pasting into word.

