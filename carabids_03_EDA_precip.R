# Download and visualize precip for Niwot

library(neonUtilities)
library(dplyr)
library(ggplot2)
library(lubridate)

# Load data
load(file="data_derived/precip_NIWO.Rdata") # NEON precip 
list2env(precip, .GlobalEnv)  
rm(precip)
load(file="data_derived/merged_C1-saddle_precip.Rdata") # Niwot LTER C1 and saddle precip
merged_precip <- merged_precip %>% 
    mutate(collectDate = as.Date(NA,format="%Y-%m-%d")) 
model_df <- read.csv("data_derived/model_df_by_species_in_sample.csv") %>%
    mutate(collectDate = as.Date(collectDate,format="%Y-%m-%d")) #load model df to help summarize precip data

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
    geom_line(aes(colour=local_site, alpha=.6))
# We see higher precip at the saddle
    
# Monthly sums
merged_precip %>%
    select(date, ppt_tot, local_site) %>%
    group_by(local_site, month = floor_date(date, unit = 'month')) %>%
    summarise(monthly_precip = sum(ppt_tot, na.rm = T)) %>%
    ggplot() +
    geom_line(aes(x = month, y = monthly_precip, colour=local_site))
# Wow, soo much more precip at the saddle than at C1


# These data seem usable, so let's summarize them to be able to plug right into the model df in the 03_zcompiled script
# Assign which collection window a precipitation date falls under, if any
uni_collDates <- data.frame("collectDate"=unique(model_df$collectDate)) %>%
    mutate(two_wk_int = interval(ymd(collectDate-14), ymd(collectDate)))
for (i in 1:nrow(merged_precip)) {
    for (j in 1:nrow(uni_collDates)) {
        if (merged_precip$date[i] %within% uni_collDates$two_wk_int[j]) {
            merged_precip$collectDate[i] <- uni_collDates$collectDate[j]
        }
    }
}

# Summarize the accumulated precip for 2 weeks prior to each collection date
summ_precip <- merged_precip %>%
    select(local_site, ppt_tot, collectDate) %>%
    group_by(local_site, collectDate) %>%
    summarise(precip_2weeks = sum(ppt_tot, na.rm = T)) %>%
    filter(is.na(collectDate)==FALSE) %>%
    mutate(nlcdClass = ifelse(local_site=="c1","evergreenForest","grasslandHerbaceous")) %>%
    ungroup() %>%
    mutate(col_year = lubridate::year(collectDate), 
           col_month = lubridate::month(collectDate),
           col_day = lubridate::day(collectDate))

# Visualize accummulated precip for each collection date
summ_precip %>%
    ggplot() +
    geom_point(aes(x = collectDate, y = precip_2weeks, colour=local_site))

save(summ_precip, file="data_derived/summarized_precip.Rdata") 

# In carabids_03_zcompiled script, we will assign tundra plots saddle precip data and forest plots C1 precip data

