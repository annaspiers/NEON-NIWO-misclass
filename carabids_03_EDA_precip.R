# Download and visualize precip for Niwot

library(neonUtilities)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats) #fct_reorder
library(arm) #display


# Load soil water content data
load(file="data_derived/precip_NIWO.Rdata")
list2env(precip, .GlobalEnv)  
rm(precip)

# Explore data


# # visualize sensor positions (five locations across Niwot)
# soil_wc$sensor_positions_00094 %>%
#   ggplot(aes(referenceLongitude, referenceLatitude)) + 
#   geom_point()
# 
# # visualize 30 minute volumetric soil water content (VSWC)
# soil_wc$SWS_30_minute %>%
#   ggplot(aes(startDateTime, VSWCMean)) + 
#   geom_point(size = .1) + 
#   facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))
# 
# # same for volumetric soil ion content (VSIC)
# soil_wc$SWS_30_minute %>%
#   ggplot(aes(startDateTime, VSICMean)) + 
#   geom_point(size = .1) + 
#   facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))