# Load required packages
library(tidyverse)
library(lubridate)

# load Niwot airtemp data from saddle and C1 (and change date format)
airtemp_saddle <- read.csv(file = "data_derived/Niwot_airtemp_saddle.csv") %>%
    mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
    filter(date >= as.Date('2015-01-01'))

airtemp_C1 <- read.csv(file = "data_derived/Niwot_airtemp_C1.csv") %>%
    mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
    filter(date >= as.Date('2015-01-01'))

# combine and clip the datasets
airtemp_combined <- rbind(airtemp_saddle, airtemp_C1)

# compare temperature at the two sites 
# Look at precip over time (by day)
ggplot(airtemp_combined, aes(x=date, y=airtemp_avg)) +
    geom_line(aes(colour=local_site, alpha=.6))

t.test(airtemp_avg ~ local_site, data = airtemp_combined)
# We see higher (and significantly different) temps at C1
# Also - sadly - it looks like the "ongoing" saddle data ends at the end of 2017
# I guess we just use C1 temps to generate GDD for the model?  Not 100% accurate, but... 
