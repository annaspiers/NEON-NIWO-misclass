# Script to calculate Growing Degree Days [GDD] for Niwot based on the NEON air temp
# dataset. 

# GDD calculation formula developed and shared by Lauren Buckley and Cesar Nufio (citation("TrenchR"))

###########################################################

# Load required packages
library(tidyverse)
library(lubridate)

# Could not install the GDD function from TrenchR (where it is housed); per developer 
# Lauren, they are working on the package, which is breaking the dl. Instead, copied 
# and pasted into a separate script, called below

# install.packages("devtools")   
# library("devtools")   
# devtools::install_github(build_vignettes = TRUE,repo = "trenchproject/TrenchR")

# Call GDD function
source(file = "carabids_03_function_GDD.R")

# Load data: Niwot LTER airtemp data
# Niwot air temp has a longer record than NEON: unfortunately, the saddle data ends in 2017.
# C1 data is posted through the end of 2019

airtemp_C1 <- read.csv(file = "data_derived/Niwot_airtemp_C1.csv") %>%
    mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
    filter(date >= as.Date('2015-01-01'))

# visualize
airtemp_C1 %>%
    ggplot() +
    geom_line(aes( x = date, y = airtemp_max, col = 'max')) +
    geom_line(aes(x = date, y = airtemp_min, col = 'min')) +
    theme_bw()
ggsave("output/temp_time.png", width = 5, height = 3, dpi = 'retina')


# calculate GDD 
# Set parameters:
LDT <- 12 # used by Nufio et al. for Niwot Ridge grasshoppers - likely not accurate for our beetles
UDT <- 33 # also used by Nufio et al. for Niwot grasshoppers

# Call function to calculate GDD for each day, then sum cumulative GDD per year

airtemp_C1 <- airtemp_C1 %>%
    group_by(date) %>%
    dplyr::select(local_site, date, airtemp_max, airtemp_min, airtemp_avg) %>%
    na.omit() %>%
    mutate(daily_GDD = degree_days(T_min = airtemp_min, T_max = airtemp_max, 
                                   LDT = LDT, UDT = UDT, method = "single.sine")) %>%
    group_by(year(date))%>%
    mutate(GDD_cum = cumsum(daily_GDD)) 

# Visualize to make sure it's working 
airtemp_C1 %>% 
    ggplot() +
    geom_line(aes ( x = date, y = GDD_cum)) +
    theme_bw()
ggsave("output/degree_days_cumulative.png", width = 5, height = 3, dpi = 'retina')

# it's not perfect (would have to clip after accumulation stopped) but it is functional/ 
# good enough for our purposes. SO:
save(airtemp_C1, file="data_derived/summarized_C1temp_GDD.Rdata") 

# Adding the data to the main df here, because it wasn't working for me in the 
# 'carabids_03_zcompiled_df.R' script. Will clean/ amend script later

# load model df so can add gdd data
model_df <- read.csv("data_derived/model_df_by_species_in_sample.csv") %>%
    mutate(collectDate = as.Date(collectDate,format="%Y-%m-%d")) 

model_df1 <- left_join(x = model_df, y = airtemp_C1, by = c('collectDate' = 'date')) %>%
    dplyr::select(-c('local_site', 'year(date)', 'airtemp_max', 'airtemp_min'))
    
# Save change
write.csv(model_df1, file="data_derived/model_df_by_species_in_sample.csv") #df
