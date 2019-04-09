# Script to create metadata.rda for the Heron 24hr experiment
# Authour: Benjamin R. Gordon
# Date: 2019-04-09

library(tidyverse)

# Load and prep data -----------------------------------------------------------
temp_data <- read_csv("./data-raw/water_temp_0.3m.csv")
salinity_data <- read_csv("./data-raw/salinity_5.4m.csv")
par_data <- read_csv("./data-raw/PAR.csv")
phenodata <- read_csv("./data-raw/phenodata-raw.csv")

# pull out readings for each collection time
temp <- c(rep(25.183, 5),
          rep(25.178, 5),
          rep(26.377, 5),
          rep(24.728, 5),
          rep(25.183, 5),
          rep(25.178, 5),
          rep(26.377, 5),
          rep(24.728, 5),
          rep(25.183, 5),
          rep(25.178, 5),
          rep(26.377, 5),
          rep(24.728, 5),
          rep(25.183, 5),
          rep(25.178, 5),
          rep(26.377, 5),
          rep(24.728, 5),
          rep(25.183, 5),
          rep(26.377, 5),
          rep(24.728, 5),
          rep(NA, 20)
          )

time_fac <- factor(c(rep("05:00", 5),
                     rep("10:00", 5),
                     rep("16:00", 5),
                     rep("21:00", 5),
                     rep("05:00", 5),
                     rep("10:00", 5),
                     rep("16:00", 5),
                     rep("21:00", 5),
                     rep("05:00", 5),
                     rep("10:00", 5),
                     rep("16:00", 5),
                     rep("21:00", 5),
                     rep("05:00", 5),
                     rep("10:00", 5),
                     rep("16:00", 5),
                     rep("21:00", 5),
                     rep("05:00", 5),
                     rep("16:00", 5),
                     rep("21:00", 5),
                     rep("PBQC", 20)
                     ))

taxon <- factor(c(rep("A_aspe", 20),
                  rep("M_aequ", 20),
                  rep("M_digi", 20),
                  rep("P_cyli", 20),
                  rep("P_dami", 15),
                  rep("PBQC", 20)
                  ))

rep <- c(rep(seq(1:5), 19),
         1, 2, 1, 10,
         seq(from = 2, to = 9, by = 1),
         rep(seq(1:2), 4)
         )

time_date <- as.POSIXct(c(rep("05:00 2014-04-21", 5),
                          rep("10:00 2014-04-21", 5),
                          rep("16:00 2014-04-21", 5),
                          rep("21:00 2014-04-21", 5),
                          rep("05:00 2014-04-21", 5),
                          rep("10:00 2014-04-21", 5),
                          rep("16:00 2014-04-21", 5),
                          rep("21:00 2014-04-21", 5),
                          rep("05:00 2014-04-21", 5),
                          rep("10:00 2014-04-21", 5),
                          rep("16:00 2014-04-21", 5),
                          rep("21:00 2014-04-21", 5),
                          rep("05:00 2014-04-21", 5),
                          rep("10:00 2014-04-21", 5),
                          rep("16:00 2014-04-21", 5),
                          rep("21:00 2014-04-21", 5),
                          rep("05:00 2014-04-21", 5),
                          rep("16:00 2014-04-21", 5),
                          rep("21:00 2014-04-21", 5),
                          rep("12:00 2014-04-21", 20)
                          ), format = "%H:%M %F")

par <- c(rep(0.115941, 5),
         rep(1495.510000, 5),
         rep(649.798000, 5),
         rep(0.011594, 5),
         rep(0.115941, 5),
         rep(1495.510000, 5),
         rep(649.798000, 5),
         rep(0.011594, 5),
         rep(0.115941, 5),
         rep(1495.510000, 5),
         rep(649.798000, 5),
         rep(0.011594, 5),
         rep(0.115941, 5),
         rep(1495.510000, 5),
         rep(649.798000, 5),
         rep(0.011594, 5),
         rep(0.115941, 5),
         rep(649.798000, 5),
         rep(0.011594, 5),
         rep(NA, 20)
         )

salinity <- c(rep(35.24663, 5),
              rep(35.23519, 5),
              rep(35.24205, 5),
              rep(35.27422, 5),
              rep(35.24663, 5),
              rep(35.23519, 5),
              rep(35.24205, 5),
              rep(35.27422, 5),
              rep(35.24663, 5),
              rep(35.23519, 5),
              rep(35.24205, 5),
              rep(35.27422, 5),
              rep(35.24663, 5),
              rep(35.23519, 5),
              rep(35.24205, 5),
              rep(35.27422, 5),
              rep(35.24663, 5),
              rep(35.24205, 5),
              rep(35.27422, 5),
              rep(NA, 20)
              )

metadata <- tibble(sample_id = phenodata$sample_name,
                   class = phenodata$class,
                   taxon, 
                   time_fac, 
                   rep,
                   temp, 
                   par, 
                   salinity,
                   time_date)

save(metadata, file = "./data/metadata.rda", compress = "bzip2")
