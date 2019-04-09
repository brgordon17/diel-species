# Script to create mzdata.rda for the heron 24 hr experiment
# Author: Benjamin R. Gordon
# Date: 2019-04-09

library(tidyverse)
library(metabolomics)
library(missForest)
library(reshape2)

# Load data --------------------------------------------------------------------
mzdata  <-  read_csv("./data-raw/mzdata-raw.csv", na = "0")
load("./data/metadata.rda")

# remove isotopes --------------------------------------------------------------
mzdata <- mzdata[-grep("[M+1]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+2]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+3]", mzdata$isotopes, fixed = TRUE),]
mzdata <- mzdata[-grep("[M+4]", mzdata$isotopes, fixed = TRUE),]

# Clean up ---------------------------------------------------------------------
mzdata <- dplyr::select(mzdata, -isotopes, -adduct, -pcgroup)
mzdata <- mzdata[-2:-32]

mz_names <- round(mzdata[, 1], 5)
mz_names <- paste("mz", mz_names$mz, sep = "_")
mz_names <- make.names(mz_names, unique = TRUE)

mzdata <- tibble::as_tibble(t(mzdata))
colnames(mzdata) <- mz_names
mzdata <- mzdata[-1, ]
mzdata <- type_convert(mzdata) # convert to dbl
mzdata <- bind_cols(metadata, mzdata)

# Impute noise and remove unreliable features ----------------------------------
round(mean(is.na(mzdata))*100, 2)
mzdata_filt <- MissingValues(mzdata[c(-1, -3:-9)],
                             column.cutoff = 0.6,
                             group.cutoff = 0.8,
                             complete.matrix = FALSE,
                             seed = 1978)
mzdata <- bind_cols(metadata,
                    mzdata_filt$output[-1])
round(mean(is.na(mzdata))*100, 2)

rm(mzdata_filt, mz_names)

# UP TO HERE
# Impute remaining missing values ----------------------------------------------
# doMC::registerDoMC()
# set.seed(1978)
# mzdata.imp <- missForest(as.data.frame(mzdata[-1:-6]),
#                          ntree = 100,
#                          verbose = TRUE,
#                          maxiter = 8,
#                          parallelize = "variables")
# 
# mzdata <- bind_cols(metadata, mzdata.imp$ximp)
# round(mean(is.na(mzdata[-1:-6]))*100, 2)
