# clear space
# rm(list = ls())

# Download and visualize ir/ biological (surface) temperature readings for Niwot
library(neonUtilities)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats) #fct_reorder
library(arm) #display

# Load temperature data
load(file = "data_derived/ir_bio_temp_NIWO.Rdata")
list2env(ir_bio_temp, .GlobalEnv)  
rm(ir_bio_temp)

#View(IRBT_30_minute)

# visualize sensor positions (three sensors at two geographic locations at Niwot)
ir_bio_temp$sensor_positions_00005 %>%
  ggplot(aes(referenceLongitude, referenceLatitude)) +
  geom_point()

# visualize 30 minute IR temperature readings 
ir_bio_temp$IRBT_30_minute %>%
  ggplot(aes(startDateTime, bioTempMean)) +
  geom_point(size = .1) +
  facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))

