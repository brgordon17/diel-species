# Feature detection and annotation using xcms and CAMERA
# Benjamin Gordon 8-JAN-2019
# See ~/gordon-C5/dev/mzdata-raw-dev.R for the developmental version of this 
# code.

library(xcms)
library(magrittr)
library(CAMERA)

# Register parallel processing -------------------------------------------------
register(bpstart(MulticoreParam(workers = 4)))

# Data import ------------------------------------------------------------------

# path to files
mzfiles <- list.files(path = "~/Big Data/LCMS Data/Heron 24 hour/mzxml_data/",
                      full.names = TRUE,
                      recursive = TRUE)

# create metadata or phenodata
phenod <- data.frame(sample_name = sub(basename(mzfiles), pattern = ".mzXML",
                                       replacement = "", fixed = TRUE),
                     class = c(rep("aspe_05", 5),
                               rep("aspe_10", 5),
                               rep("aspe_16", 5),
                               rep("aspe_21", 5),
                               rep("aequ_05", 5),
                               rep("aequ_10", 5),
                               rep("aequ_16", 5),
                               rep("aequ_21", 5),
                               rep("digi_05", 5),
                               rep("digi_10", 5),
                               rep("digi_16", 5),
                               rep("digi_21", 5),
                               rep("cyli_05", 5),
                               rep("cyli_10", 5),
                               rep("cyli_16", 5),
                               rep("cyli_21", 5),
                               rep("dami_05", 5),
                               rep("dami_16", 5),
                               rep("dami_21", 5),
                               rep("PBQC_aequ", 2),
                               rep("PBQC_all", 10),
                               rep("PBQC_aspe", 2),
                               rep("PBQC_cyli", 2),
                               rep("PBQC_dami", 2),
                               rep("PBQC_digi", 2)
                     ),
                     stringsAsFactors = FALSE)

# load raw data and save
raw_data <- readMSData(files = mzfiles, pdata = new("NAnnotatedDataFrame", 
                                                    phenod),
                       mode = "onDisk")
saveRDS(raw_data, "./data-raw/temp-saves/OnDiskMSnExp.rds")

# data inspection --------------------------------------------------------------
# ID problematic LCMS runs using boxplots of the TIC
tics <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tics, ylab = "intensity", main = "TIC")
rm(tics)

# Peak detection ---------------------------------------------------------------
# define the parameters 
cwp <- CentWaveParam(ppm = 30, 
                     peakwidth = c(10, 60),
                     snthresh = 10,
                     prefilter = c(3, 1100),
                     mzCenterFun = "wMean",
                     integrate = 1,
                     mzdiff = -0.001,
                     fitgauss = FALSE)
# Do peak detection and save
mzdata <- findChromPeaks(raw_data, param = cwp)
rm(cwp)
saveRDS(mzdata, "./data-raw/temp-saves/XCMSnExp.rds")

# Alignment --------------------------------------------------------------------
mzdata <- adjustRtime(mzdata, param = ObiwarpParam(binSize = 0.5))
# plot retention time adjustment
plotAdjustedRtime(mzdata)

# Correspondence or Grouping ---------------------------------------------------
pdp <- PeakDensityParam(mzdata$class,
                        bw = 5,
                        minFraction = 0.5,
                        minSamples = 1,
                        binSize = 0.025, #same as mzwid
                        maxFeatures = 50)
mzdata <- groupChromPeaks(mzdata, param = pdp)

# fill missing peaks -----------------------------------------------------------
mzdata <- fillChromPeaks(mzdata)
saveRDS(mzdata, "./data-raw/temp-saves/XCMSnExp.rds")
processHistory(mzdata)

# Isotope and adduct annotation ------------------------------------------------
# convert to xcmsset object
mzdata <- as(mzdata, "xcmsSet")
# Create report with with isotope and adduct info
mzdata <- annotate(mzdata, 
                   perfwhm = 0.7, 
                   cor_eic_th = 0.75, 
                   ppm = 10, 
                   polarity = "positive")
# save
saveRDS(mzdata, "./data-raw/temp-saves/xsAnnotate.rds")

# Tidy and save ----------------------------------------------------------------
# create tibble
mzdata_raw <- tibble::as_tibble(CAMERA::getPeaklist(mzdata))
# reattach sample names
names(mzdata_raw)[33:147] <- phenod$sample_name
# save mzdata-raw as csv
readr::write_csv(mzdata_raw, "./data-raw/mzdata-raw.csv")
# save phenodata as csv
readr::write_csv(phenod, "./data-raw/phenodata-raw.csv")
# remove temporary saves
file.remove(list.files("./data-raw/temp-saves", full.names = TRUE))
