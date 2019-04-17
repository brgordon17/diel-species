# Script to analyse the quenching data.
# Author: Benjamin R. Gordon 
# Date: 2019-04-16

# PCA --------------------------------------------------------------------------
library(tidyverse)

# load data
load("./data/snapdata.rda")

# PCA
set.seed(1978)
pca <- stats::prcomp(snapdata[-1:-2], 
                     scale = FALSE, 
                     center = TRUE)

# Plotting
exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(snapdata[1:2], pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")

ggplot(data = filter(scores, class == "methanol" | class == "nitrogen"),
       aes(x = PC1,
           y = PC2,
           color = class,
           fill = class,
           shape = class)) +
  geom_point(size = 2.5,
             stroke = 0.7,
             position = position_jitter(width = 0.01 * diff(range(scores$PC1)),
                                        height = 0.01 * diff(range(scores$PC2))
             )) +
  labs(x = x_lab, y = y_lab) +
  scale_shape_manual(values = c(21:25)) +
  scale_color_manual(values = 
                       grDevices::adjustcolor(gordon01::qual_colours[c(6, 2)],
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = 
                      grDevices::adjustcolor(gordon01::qual_colours[c(6, 2)],
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey65"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 10))

# RLA --------------------------------------------------------------------------
load("./data/snapdata.rda")

ggplot(data = reshape2::melt(filter(snapdata, 
                                    class == "methanol" | class == "nitrogen")),
       aes(x = sample_id, 
           y = log(value, 2) - median(log(value, 2)),
           fill = class)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = gordon01::qual_colours[c(6, 2)]) +
  scale_x_discrete(name = "Replicate") +
  scale_y_continuous(name = "Relative log Abundance") +
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10, colour = "grey65"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank())

# RF model ---------------------------------------------------------------------
library(caret)
library(tidyverse)

# Load and prep data 
load("./data/snapdata.rda")
snapdata <- filter(snapdata, class == "methanol" | class == "nitrogen")
snapdata <- data.frame(droplevels(snapdata))

# Model of two classes ---------------------------------------------------------
tunegrid <- expand.grid(.mtry = seq(from = 25, to = 100, by = 25))

# set seeds
set.seed(1978)
seeds <- vector(mode = "list", length = 31)
for(i in 1:30) seeds[[i]] <- sample.int(1000, length(tunegrid[,1]))
seeds[[31]] <- sample.int(1000, 1)

# check for best k. Ideally, 4/20 samples (1 of each class) to be held out.
# createFolds(y = snapdata$class, k = 5)

# train control
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 5,
                            repeats = 6,
                            summaryFunction = defaultSummary,
                            seeds = seeds,
                            savePredictions = "all",
                            selectionFunction = "oneSE")

# model
doMC::registerDoMC()
set.seed(1978)
mzrf <- caret::train(x = snapdata[-1:-2],
                     y = snapdata$class,
                     method = "rf",
                     trControl = ctrl,
                     #preProc = c("center", "scale"),
                     allowParallel = TRUE,
                     importance = TRUE,
                     tuneGrid = tunegrid)

confusionMatrix(mzrf$pred$pred[mzrf$pred$mtry == mzrf$bestTune$mtry],
                mzrf$pred$obs[mzrf$pred$mtry == mzrf$bestTune$mtry])

saveRDS(mzrf, "./dev/quenching_model.rds")

# Pull out important variables
mzrf_impvars <- varImp(mzrf, scale = FALSE)
mzrf_impvars <- mzrf_impvars$importance
mzrf_impvars <- cbind(importance = apply(mzrf_impvars, 1, max), mzrf_impvars)
mzrf_impvars <- mzrf_impvars[order(-mzrf_impvars$importance), ,drop = FALSE]
mzrf_impvars <- as_tibble(mzrf_impvars[1:9, ], rownames = "mz")

# plot top 6 variable intensities
sum_snapdata <- 
  snapdata %>% 
  filter(class == "methanol" | class == "nitrogen") %>%
  select(1:2, c(mzrf_impvars$mz)) %>%
  reshape2::melt(.) %>%
  as_tibble()
  
sum_snapdata$variable <- str_replace_all(sum_snapdata$variable, "mz_", "")

ggplot(sum_snapdata, aes(x = class,
                         y = value/10^4,
                         colour = class,
                         fill = class)) +
  geom_boxplot(alpha = 0.1) +
  facet_wrap(vars(variable)) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = expression(Intensity%*%10^4)) +
  scale_colour_manual(values = gordon01::qual_colours[c(6, 2)]) +
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
