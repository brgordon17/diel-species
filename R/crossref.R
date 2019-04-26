# Script to cross reference important diurnal features with the coral
# research literature.
# Author: Benjamin R. Gordon
# Date: 2019-04-26

library(tidyverse)
library(caret)

# Load data --------------------------------------------------------------------
mzrf_aspe <- read_rds("./dev/mzrf_aspe.rds")
mzrf_aequ <- read_rds("./dev/mzrf_aequ.rds")
mzrf_digi <- read_rds("./dev/mzrf_digi.rds")
mzrf_cyli <- read_rds("./dev/mzrf_cyli.rds")
mzrf_dami <- read_rds("./dev/mzrf_dami.rds")

load("./data/litmz.rda")

mzdata_raw  <-  readr::read_csv("./data-raw/mzdata-raw.csv", na = "0")
colnames(mzdata_raw)[1] <- "mz_raw"

# Crossref A. aspera -----------------------------------------------------------
# retrieve impvars
impvars <- varImp(mzrf_aspe, scale = FALSE)
impvars <- impvars$importance
impvars <- 
  bind_cols(mz = as.numeric(str_replace(rownames(impvars), "mz_", "")),
            importance = apply(impvars, 1, max))
impvars <- impvars[order(-impvars$importance), ,drop = FALSE]
impvars <- impvars[1:20, ]

# Add variables for 50 ppm error ranges 
ppm <- 50
matches <-
  impvars %>%
  rowwise() %>%
  mutate(adduct = "none",
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz
matches <-
  matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "A_aspe") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# identify adducts from impvars
adduct_matches <- 
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(mz,
         importance,
         adduct) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference adducts with litmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "A_aspe") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# Join
aspe_matches <- bind_rows(matches, adduct_matches)

# Crossref M_aequ --------------------------------------------------------------
# retrieve impvars
impvars <- varImp(mzrf_aequ, scale = FALSE)
impvars <- impvars$importance
impvars <- 
  bind_cols(mz = as.numeric(str_replace(rownames(impvars), "mz_", "")),
            importance = apply(impvars, 1, max))
impvars <- impvars[order(-impvars$importance), ,drop = FALSE]
impvars <- impvars[1:20, ]

# Add variables for 50 ppm error ranges 
ppm <- 50
matches <-
  impvars %>%
  rowwise() %>%
  mutate(adduct = "none",
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz
matches <-
  matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "M_aequ") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# identify adducts from impvars
adduct_matches <- 
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(mz,
         importance,
         adduct) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference adducts with litmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "M_aequ") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# Join
aequ_matches <- bind_rows(matches, adduct_matches)

# Crossref M_digi --------------------------------------------------------------
# retrieve impvars
impvars <- varImp(mzrf_digi, scale = FALSE)
impvars <- impvars$importance
impvars <- 
  bind_cols(mz = as.numeric(str_replace(rownames(impvars), "mz_", "")),
            importance = apply(impvars, 1, max))
impvars <- impvars[order(-impvars$importance), ,drop = FALSE]
impvars <- impvars[1:20, ]

# Add variables for 50 ppm error ranges 
ppm <- 50
matches <-
  impvars %>%
  rowwise() %>%
  mutate(adduct = "none",
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz
matches <-
  matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "M_digi") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# identify adducts from impvars
adduct_matches <- 
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(mz,
         importance,
         adduct) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference adducts with litmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "M_digi") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# Join
digi_matches <- bind_rows(matches, adduct_matches)

# Crossref M_cyli --------------------------------------------------------------
# retrieve impvars
impvars <- varImp(mzrf_cyli, scale = FALSE)
impvars <- impvars$importance
impvars <- 
  bind_cols(mz = as.numeric(str_replace(rownames(impvars), "mz_", "")),
            importance = apply(impvars, 1, max))
impvars <- impvars[order(-impvars$importance), ,drop = FALSE]
impvars <- impvars[1:20, ]

# Add variables for 50 ppm error ranges 
ppm <- 50
matches <-
  impvars %>%
  rowwise() %>%
  mutate(adduct = "none",
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz
matches <-
  matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "P_cyli") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# identify adducts from impvars
adduct_matches <- 
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(mz,
         importance,
         adduct) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference adducts with litmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "P_cyli") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# Join
cyli_matches <- bind_rows(matches, adduct_matches)

# Crossref P_dami --------------------------------------------------------------
# retrieve impvars
impvars <- varImp(mzrf_dami, scale = FALSE)
impvars <- impvars$importance
impvars <- 
  bind_cols(mz = as.numeric(str_replace(rownames(impvars), "mz_", "")),
            importance = apply(impvars, 1, max))
impvars <- impvars[order(-impvars$importance), ,drop = FALSE]
impvars <- impvars[1:20, ]

# Add variables for 50 ppm error ranges 
ppm <- 50
matches <-
  impvars %>%
  rowwise() %>%
  mutate(adduct = "none",
         mz_neutral = mz - 1.007276,
         mz_low = mz_neutral - (mz * ppm/10^6),
         mz_high = mz_neutral + (mz * ppm/10^6)) %>%
  ungroup()

# Cross reference with litmz
matches <-
  matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "P_dami") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# identify adducts from impvars
adduct_matches <- 
  impvars %>%
  mutate(dummy = TRUE) %>%
  left_join(mzdata_raw %>% mutate(dummy = TRUE))  %>%
  filter(near(mz, mz_raw, tol = .0001)) %>%
  select(mz,
         importance,
         adduct) %>%
  mutate(mz_neutral = as.numeric(str_sub(adduct, str_length(adduct)-6, -1))) %>%
  mutate(adduct = str_sub(adduct, 1, str_length(adduct)-8)) %>%
  filter(!is.na(mz_neutral))

# add variables for 50ppm error
adduct_matches <-
  adduct_matches %>%
  rowwise() %>%
  mutate(mz_low = mz_neutral - (mz_neutral * ppm/10^6),
         mz_high = mz_neutral + (mz_neutral * ppm/10^6)) %>%
  ungroup()

# Cross reference adducts with litmz
adduct_matches <-
  adduct_matches %>%
  mutate(dummy = TRUE) %>%
  left_join(litmz %>% mutate(dummy = TRUE))  %>%
  filter(monoiso_mass <= mz_high, monoiso_mass >= mz_low) %>%
  mutate(sample = "P_dami") %>%
  select(-dummy,
         -mz_neutral,
         -mz_low,
         -mz_high,
         -importance) %>%
  select(sample, everything())

# Join
dami_matches <- bind_rows(matches, adduct_matches)

# save -------------------------------------------------------------------------
all_matches <- bind_rows(aspe_matches, 
                         aequ_matches, 
                         cyli_matches, 
                         digi_matches, 
                         dami_matches)

write_rds(all_matches, "./dev/impvar_matches.rds")
