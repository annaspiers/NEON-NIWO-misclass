# Download and visualize soil water content and ion content for Niwot
# this includes 30 minute and 1 minute data at 5 locations
# Max Joseph (reorganized by Anna)

library(neonUtilities)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats) #fct_reorder
library(arm) #display

# Load carabid data
# AIS load soil water content data
#load(file="data_raw/beetles/carabids_NIWO.Rdata")
#list2env(carabids_NIWO, .GlobalEnv)  


# visualize sensor positions (five locations across Niwot)
soil_wc$sensor_positions_00094 %>%
  ggplot(aes(referenceLongitude, referenceLatitude)) + 
  geom_point()

# visualize 30 minute volumetric soil water content (VSWC)
soil_wc$SWS_30_minute %>%
  ggplot(aes(startDateTime, VSWCMean)) + 
  geom_point(size = .1) + 
  facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))

# same for volumetric soil ion content (VSIC)
soil_wc$SWS_30_minute %>%
  ggplot(aes(startDateTime, VSICMean)) + 
  geom_point(size = .1) + 
  facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))