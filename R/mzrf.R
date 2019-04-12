# Script to perform random forests predicition of time and taxon for Heron 24 h
# experiment
# Author: Benjamin R. Gordon 
# Date: 2019-04-10

library(tidyverse)
library(caret)

# Load and prep data -----------------------------------------------------------
load("./data/mzdata.rda")
mzdata <- filter(mzdata, taxon != "PBQC")
mzdata <- data.frame(droplevels(mzdata))

# Partition data into training, test and validation sets
set.seed(16)
train_index <- createDataPartition(mzdata$class,
                                   p = 0.8,
                                   list = FALSE)
train_data <- mzdata[train_index, ]
test_data  <- mzdata[-train_index, ]

# Model of taxonomy ------------------------------------------------------------
tunegrid <- expand.grid(.mtry = seq(from = 25, to = 200, by = 25))

# set seeds
set.seed(1978)
seeds <- vector(mode = "list", length = 31)
for(i in 1:30) seeds[[i]] <- sample.int(1000, length(tunegrid[,1]))
seeds[[31]] <- sample.int(1000, 1)

# train control
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3,
                            summaryFunction = defaultSummary,
                            seeds = seeds,
                            savePredictions = "all",
                            selectionFunction = "oneSE")

# model
doMC::registerDoMC()
set.seed(1978)
mzrf_taxon <- caret::train(x = train_data[-1:-9],
                           y = train_data$taxon,
                           method = "rf",
                           trControl = ctrl,
                           #preProc = c("center", "scale"),
                           allowParallel = TRUE,
                           importance = TRUE,
                           tuneGrid = tunegrid
                           )

preds <- confusionMatrix(predict(mzrf_taxon, newdata = test_data),
                         test_data$taxon)

# Model of diurnal variation ---------------------------------------------------
# Load and prep data
load("./data/mzdata.rda")

mzdata_aspe <- filter(mzdata, taxon == "A_aspe")
mzdata_aspe <- data.frame(droplevels(mzdata_aspe))

mzdata_aequ <- filter(mzdata, taxon == "M_aequ")
mzdata_aequ <- data.frame(droplevels(mzdata_aequ))

mzdata_digi <- filter(mzdata, taxon == "M_digi")
mzdata_digi <- data.frame(droplevels(mzdata_digi))

mzdata_cyli <- filter(mzdata, taxon == "P_cyli")
mzdata_cyli <- data.frame(droplevels(mzdata_cyli))

mzdata_dami <- filter(mzdata, taxon == "P_dami")
mzdata_dami <- data.frame(droplevels(mzdata_dami))

# tuning grid
tunegrid <- expand.grid(.mtry = seq(from = 25, to = 200, by = 25))

# set seeds
set.seed(1978)
seeds <- vector(mode = "list", length = 31)
for(i in 1:30) seeds[[i]] <- sample.int(1000, length(tunegrid[,1]))
seeds[[31]] <- sample.int(1000, 1)

# train control
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 3,
                            summaryFunction = defaultSummary,
                            seeds = seeds,
                            savePredictions = "all",
                            selectionFunction = "oneSE")

# A_aspe model
doMC::registerDoMC()
set.seed(1978)
mzrf_aspe <- caret::train(x = mzdata_aspe[-1:-9],
                           y = mzdata_aspe$time_fac,
                           method = "rf",
                           trControl = ctrl,
                           #preProc = c("center", "scale"),
                           allowParallel = TRUE,
                           importance = TRUE,
                           tuneGrid = tunegrid
                           )

confusionMatrix(mzrf_aspe$pred$pred[mzrf_aspe$pred$mtry == 
                                      mzrf_aspe$bestTune$mtry],
                mzrf_aspe$pred$obs[mzrf_aspe$pred$mtry ==
                                     mzrf_aspe$bestTune$mtry])

saveRDS(mzrf_aspe, "./dev/mzrf_aspe.rds")

# M_aequ model
doMC::registerDoMC()
set.seed(1978)
mzrf_aequ <- caret::train(x = mzdata_aequ[-1:-9],
                          y = mzdata_aequ$time_fac,
                          method = "rf",
                          trControl = ctrl,
                          #preProc = c("center", "scale"),
                          allowParallel = TRUE,
                          importance = TRUE,
                          tuneGrid = tunegrid
                          )

confusionMatrix(mzrf_aequ$pred$pred[mzrf_aequ$pred$mtry == 
                                      mzrf_aequ$bestTune$mtry],
                mzrf_aequ$pred$obs[mzrf_aequ$pred$mtry ==
                                     mzrf_aequ$bestTune$mtry])

saveRDS(mzrf_aequ, "./dev/mzrf_aequ.rds")

# M_digi model
doMC::registerDoMC()
set.seed(1978)
mzrf_digi <- caret::train(x = mzdata_digi[-1:-9],
                          y = mzdata_digi$time_fac,
                          method = "rf",
                          trControl = ctrl,
                          #preProc = c("center", "scale"),
                          allowParallel = TRUE,
                          importance = TRUE,
                          tuneGrid = tunegrid
                          )

confusionMatrix(mzrf_digi$pred$pred[mzrf_digi$pred$mtry == 
                                      mzrf_digi$bestTune$mtry],
                mzrf_digi$pred$obs[mzrf_digi$pred$mtry ==
                                     mzrf_digi$bestTune$mtry])

saveRDS(mzrf_digi, "./dev/mzrf_digi.rds")

# P_cyli model
doMC::registerDoMC()
set.seed(1978)
mzrf_cyli <- caret::train(x = mzdata_cyli[-1:-9],
                          y = mzdata_cyli$time_fac,
                          method = "rf",
                          trControl = ctrl,
                          #preProc = c("center", "scale"),
                          allowParallel = TRUE,
                          importance = TRUE,
                          tuneGrid = tunegrid
                          )

confusionMatrix(mzrf_cyli$pred$pred[mzrf_cyli$pred$mtry == 
                                      mzrf_cyli$bestTune$mtry],
                mzrf_cyli$pred$obs[mzrf_cyli$pred$mtry ==
                                     mzrf_cyli$bestTune$mtry])

saveRDS(mzrf_cyli, "./dev/mzrf_cyli.rds")

# P_dami model
doMC::registerDoMC()
set.seed(1978)
mzrf_dami <- caret::train(x = mzdata_dami[-1:-9],
                          y = mzdata_dami$time_fac,
                          method = "rf",
                          trControl = ctrl,
                          #preProc = c("center", "scale"),
                          allowParallel = TRUE,
                          importance = TRUE,
                          tuneGrid = tunegrid
                          )

confusionMatrix(mzrf_dami$pred$pred[mzrf_dami$pred$mtry == 
                                      mzrf_dami$bestTune$mtry],
                mzrf_dami$pred$obs[mzrf_dami$pred$mtry ==
                                     mzrf_dami$bestTune$mtry])

saveRDS(mzrf_dami, "./dev/mzrf_dami.rds")
