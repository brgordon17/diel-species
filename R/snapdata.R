# Script to create snapdata.rda for the heron 24 hr experiment
# Author: Benjamin R. Gordon
# Date: 2019-04-15

library(tidyverse)
library(metabolomics)
library(missForest)
library(reshape2)

# Load and prep data --------------------------------------------------------------------
snapdata <- read_csv("./data-raw/snapdata-raw.csv", na = "0")
phenodata <- read_csv("./data-raw/snap-phenodata-raw.csv", na = "0")
colnames(phenodata)[1] <- "sample_id"

# remove isotopes --------------------------------------------------------------
snapdata <- snapdata[-grep("[M+1]", snapdata$isotopes, fixed = TRUE),]
snapdata <- snapdata[-grep("[M+2]", snapdata$isotopes, fixed = TRUE),]
snapdata <- snapdata[-grep("[M+3]", snapdata$isotopes, fixed = TRUE),]
snapdata <- snapdata[-grep("[M+4]", snapdata$isotopes, fixed = TRUE),]

# Clean up ---------------------------------------------------------------------
snapdata <- select(snapdata, -isotopes:-pcgroup)
snapdata <- select(snapdata, -mzmin:-npeaks)
snapdata <- select(snapdata, -methanol:-snap)

mz_names <- round(snapdata[, 1], 5)
mz_names <- paste("mz", mz_names$mz, sep = "_")
mz_names <- make.names(mz_names, unique = TRUE)

snapdata <- tibble::as_tibble(t(snapdata))
colnames(snapdata) <- mz_names
snapdata <- snapdata[-1, ]
snapdata <- type_convert(snapdata) # convert to dbl

snapdata <- bind_cols(phenodata, snapdata)


# Impute noise and remove unreliable features ----------------------------------
round(mean(is.na(snapdata))*100, 2)
snapdata_filt <- MissingValues(snapdata[-1],
                             column.cutoff = 0.6,
                             group.cutoff = 0.7,
                             complete.matrix = FALSE,
                             seed = 1978)
snapdata <- bind_cols(phenodata, snapdata_filt$output[-1])
round(mean(is.na(snapdata))*100, 2)

rm(snapdata_filt, mz_names)

# Impute remaining missing values ----------------------------------------------
doMC::registerDoMC()
set.seed(1978)
snapdata.imp <- missForest(as.data.frame(snapdata[-1:-2]),
                         ntree = 100,
                         verbose = TRUE,
                         maxiter = 8,
                         parallelize = "variables")

snapdata <- bind_cols(phenodata, snapdata.imp$ximp)
round(mean(is.na(snapdata[-1:-9]))*100, 2)
saveRDS(snapdata.imp, "./dev/snapdata_imp.rds")

# Check data structure ---------------------------------------------------------
# RLA plot
ggplot(data = reshape2::melt(filter(snapdata, 
                                    class == "methanol" | class == "snap")),
       aes(x = sample_id, 
           y = log(value, 2) - median(log(value, 2)),
           fill = class)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = gordon01::qual_colours) +
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

# PCA plot 
set.seed(1978)
pca <- stats::prcomp(snapdata[-1:-2], 
                     scale = FALSE, 
                     center = TRUE)

exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(phenodata, pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")

ggplot(data = filter(scores, class == "methanol" | class == "snap"),
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
  scale_color_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
                                                     alpha.f = 0.9)) +
  scale_fill_manual(values = grDevices::adjustcolor(gordon01::qual_colours,
                                                    alpha.f = 0.5)) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text = element_text(size = 10, colour = "grey50"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 10))

# write data -------------------------------------------------------------------
save(snapdata, file = "./data/snapdata.rda", compress = "bzip2")
