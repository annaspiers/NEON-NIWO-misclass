# Download and visualize precip for Niwot

library(neonUtilities)
library(dplyr)
library(ggplot2)
library(lubridate)


# Load soil water content data
load(file="data_derived/precip_NIWO.Rdata")
list2env(precip, .GlobalEnv)  
rm(precip)


# Explore data

# look at precip over time
PRIPRE_30min %>%
    ggplot() +
    geom_line(aes(x = startDateTime, y = priPrecipBulk))

# daily sums
PRIPRE_30min %>%
    group_by(day = floor_date(startDateTime, unit = 'day')) %>%
    summarise(daily_precip = sum(priPrecipBulk, na.rm = T)) %>%
    ggplot() +
    geom_line(aes(x = day, y = daily_precip))

# yearly total precipitation
PRIPRE_30min %>%
    group_by(lubridate::year(startDateTime)) %>%
    summarise(year_precip = sum(priPrecipBulk, na.rm = T))


# Seems like the precipitation gauge has too many missing time points to be
# super useful to us. Perhaps we could take an average across the data that are
# available and use it to look at seasonality, but it seems like a lot of work
# for something that wouldn't be super reliable...
