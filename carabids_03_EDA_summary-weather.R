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
load(file = "data_derived/summ_weath_NIWO.Rdata")
list2env(summ_weath, .GlobalEnv)  
rm(summ_weath)

View(wss_daily_temp)

load(file = "data_raw/air_temp_NIWO.Rdata")
list2env(air_temp, .GlobalEnv)  
rm(air_temp)

# visualize sensor positions (two geographic locations at Niwot)
summ_weath$sensor_positions_00001 %>%
    ggplot(aes(referenceLongitude, referenceLatitude)) +
    geom_point()

# visualize daily air temp readings 
wss_daily_temp %>%
    ggplot(aes(date, wssTempTripleMean)) +
    geom_point(size = .1) +
    facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))


#OR
TAAT_30min %>%
    ggplot(aes(startDateTime, tempTripleMean)) +
    geom_point(size = .1) +
    facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))

# visualize daily precip readings 
wss_daily_precip %>%
    ggplot(aes(date, wssPrecipTotal)) +
    geom_point(size = .1) +
    facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))

