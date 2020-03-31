library(geoNEON)
library(httr)
library(jsonlite)
library(ggplot2)
library(raster)

source('query_api_for_coords.R')
source("extr_car_plot_coords.R")

load("data_raw/carabids_NIWO.Rdata")

car_plot_coords = extr_car_plot_coords() #gets the plot level coords using the geoNEON package

#make list of strings to paste on the url to query NEON's API
trap_strings = carabid_abund$bet_fielddata %>%
    dplyr::select(namedLocation, trapID) %>%
    mutate(string_to_query = paste(namedLocation, trapID, sep = ".")) %>%
    dplyr::select(string_to_query) %>%
    unique() %>%
    pull(string_to_query)

car_trap_coords = query_api_for_coords(strings_to_query = trap_strings) #queries the NEON api


car_coords = car_trap_coords %>% separate(locationName, into = c('plotID', 'junk', 'junk2', 'trapID'), sep = "\\.", remove = FALSE) %>%
    select(-junk, -junk2) %>%
    left_join(car_plot_coords, by = 'plotID')


#plot over the dem to check
dem = raster("data_derived/neonDEM_50x50.grd")
dem_df = dem %>% as.data.frame(xy = TRUE)

car_coords %>%
    ggplot() +
    #geom_raster(data = dem_df, aes(x = x, y = y, fill = layer)) +
    geom_point(aes(x = plot.Easting, y = plot.Northing), shape = 1, size = 3, color = 'red') +
    geom_point(aes(x = trap.Easting, y = trap.Northing), shape = 1) 
# hard to see without zooming, but looks good without the dem

#save the coords data for later use
save(car_coords, file = "data_derived/car_coords.RData")






