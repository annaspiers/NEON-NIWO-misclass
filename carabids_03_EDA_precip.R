# Download and visualize precip for Niwot

library(neonUtilities)
library(dplyr)
library(ggplot2)
library(lubridate)

# Load data
# Load soil water content data
load(file="data_derived/precip_NIWO.Rdata") # NEON precip 
list2env(precip, .GlobalEnv)  
rm(precip)
load(file="data_derived/merged_C1-saddle_precip.Rdata") # Niwot LTER C1 and saddle precip

# NEON Precipitation gauge ------------------------------------------------

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


# Seems like the NEON precipitation gauge has too many missing time points to be
# super useful to us. Perhaps we could take an average across the data that are
# available and use it to look at seasonality, but it seems like a lot of work
# for something that wouldn't be super reliable...



# Niwot Ridge LTER precipitation ------------------------------------------

# Merge df's

# Look at precip over time (by day)
ggplot(merged_precip, aes(x=date, y=ppt_tot)) +
    geom_smooth(aes(colour=local_site, alpha=.6))
# We see higher precip at the saddle
    
# Monthly sums
merged_precip %>%
    select(date, ppt_tot, local_site) %>%
    group_by(local_site, month = floor_date(date, unit = 'month')) %>%
    summarise(monthly_precip = sum(ppt_tot, na.rm = T)) %>%
    ggplot() +
    geom_line(aes(x = month, y = monthly_precip, colour=local_site))
# Wow, soo much more precip at the saddle than at C1

