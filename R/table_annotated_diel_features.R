# Script for producing table data of annotated diel features
# Author: Benjamin R. Gordon
# Date: 2019-04-12
# Table built in excel from results below

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
