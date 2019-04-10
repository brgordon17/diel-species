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

# Impute remaining missing values ----------------------------------------------
doMC::registerDoMC()
set.seed(1978)
mzdata.imp <- missForest(as.data.frame(mzdata[-1:-9]),
                         ntree = 100,
                         verbose = TRUE,
                         maxiter = 8,
                         parallelize = "variables")

mzdata <- bind_cols(metadata, mzdata.imp$ximp)
round(mean(is.na(mzdata[-1:-9]))*100, 2)
saveRDS(mzdata.imp, "./dev/mzdata_imp.rds")

# Check data structure ---------------------------------------------------------
# RLA plot
ggplot(data = reshape2::melt(filter(mzdata[-5:-9], taxon != "PBQC")),
       aes(x = sample_id, 
           y = log(value, 2) - median(log(value, 2)),
           fill = taxon)) +
  geom_boxplot(outlier.alpha = 0.4,
               outlier.size = 1) +
  scale_fill_manual(values = gordon01::qual_colours) +
  scale_x_discrete(name = NULL) +
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
pca <- stats::prcomp(mzdata[-1:-9], 
                     scale = FALSE, 
                     center = TRUE)

exp_var <- summary(pca)$importance[2 ,]
scores <- data.frame(metadata, pca$x)
x_lab <- paste("PC1", " (", round(exp_var[1] * 100, 2), "%)", sep =  "")
y_lab <- paste("PC2", " (", round(exp_var[2] * 100, 2), "%)", sep =  "")

ggplot(data = filter(scores, taxon != "PBQC"),
         aes(x = PC1,
             y = PC2,
             color = taxon,
             fill = taxon,
             shape = taxon)) +
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
save(mzdata, file = "./data/mzdata.rda", compress = "bzip2")
