# Here, we explore the NEON Slope and Aspect data products at the Niwot Ridge site and summarize slope and aspect at the trap- and plot-level.

library(dplyr)
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)


# Load data ---------------------------------------------------------------

slope_raster_2017 = raster("data_derived/neonSlope_1x1_2017.grd")
aspect_raster_2017 = raster("data_derived/neonAspect_1x1_2017.grd")
slope_raster_2018 = raster("data_derived/neonSlope_1x1_2018.grd")
aspect_raster_2018 = raster("data_derived/neonAspect_1x1_2018.grd")
slope_raster_2019 = raster("data_derived/neonSlope_1x1_2019.grd")
aspect_raster_2019 = raster("data_derived/neonAspect_1x1_2019.grd")

load("data_derived/car_coords.RData")
SA_proj = projection(slope_raster_2017) #the same across years
trap_spdf <- SpatialPointsDataFrame(coords = car_coords %>%
                                        dplyr::select(trap.Easting, trap.Northing), 
                                    data=car_coords %>%
                                        dplyr::select(trap.Easting, trap.Northing), 
                                    proj4string = CRS(SA_proj))
plot_spdf <- SpatialPointsDataFrame(coords = car_coords %>%
                                        dplyr::select(plot.Easting, plot.Northing), 
                                    data=car_coords %>%
                                        dplyr::select(plot.Easting, plot.Northing), 
                                    proj4string = CRS(SA_proj))


# Crop rasters ---------------------------------------------------------

buff <- 300 #set arbitrary buffer - I think this is in m per the earthdatascience resource
carabid_extent <- extent(min(car_coords$trap.Easting) - buff, #xmin
                         max(car_coords$trap.Easting) + buff, #xmax
                         min(car_coords$trap.Northing) - buff, #ymin
                         max(car_coords$trap.Northing) + buff) #ymax
slope_2017_crop <- crop(slope_raster_2017, carabid_extent)
aspect_2017_crop <- crop(aspect_raster_2017, carabid_extent)
slope_2018_crop <- crop(slope_raster_2018, carabid_extent)
aspect_2018_crop <- crop(aspect_raster_2018, carabid_extent)
slope_2019_crop <- crop(slope_raster_2019, carabid_extent)
aspect_2019_crop <- crop(aspect_raster_2019, carabid_extent)


# What's the difference in LAI between years? -----------------------------

par(mfrow=c(1, 3))
plot(slope_2017_crop, main="2017 Slope")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)
plot(slope_2018_crop, main="2018 Slope")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)
plot(slope_2019_crop, main="2019 Slope")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)

hist(slope_2017_crop, main="2017 Slope")
hist(slope_2018_crop, main="2018 Slope")
hist(slope_2019_crop, main="2019 Slope")
# Slope of all 3 years looks pretty much the same. For this project, let's just choose 2017 since it's within the carabid data span of 2015-18

par(mfrow=c(1, 3))
plot(aspect_2017_crop, main="2017 Aspect")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)
plot(aspect_2018_crop, main="2018 Aspect")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)
plot(aspect_2019_crop, main="2019 Aspect")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)

hist(aspect_2017_crop, main="2017 Aspect")
hist(aspect_2018_crop, main="2018 Aspect")
hist(aspect_2019_crop, main="2019 Aspect")
# Aspect, too, looks the same through the years. Choose 2017


# Summarize slope/aspect at trap- and plot-level ------------------------------

# Plots are 40m x 40m. Traps lay within the plot at the middle of each edge.  If they laid exactly on plot edge, they would be about 28m from the traps on the neighboring sides. How large should we summarize LAI for each trap? As to not overlap with neighboring traps' buffers, the buffer should be at most 14m. Max, min, avg LAI? for now choose avg

#AIS experiment with different buffer sizes. 
# Mean slope summarized to trap-level
trap_slope_2017 <- raster::extract(x = slope_2017_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
trap_slope_2018 <- raster::extract(x = slope_2018_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
trap_slope_2019 <- raster::extract(x = slope_2019_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)

# Mean aspect summarized to trap-level
trap_aspect_2017 <- raster::extract(x = aspect_2017_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
trap_aspect_2018 <- raster::extract(x = aspect_2018_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
trap_aspect_2019 <- raster::extract(x = aspect_2019_crop, #raster
                               y = trap_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)

# Mean slope summarized to plot-level
plot_slope_2017 <- raster::extract(x = slope_2017_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
plot_slope_2018 <- raster::extract(x = slope_2018_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
plot_slope_2019 <- raster::extract(x = slope_2019_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)

# Mean aspect summarized to trap-level
plot_aspect_2017 <- raster::extract(x = aspect_2017_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
plot_aspect_2018 <- raster::extract(x = aspect_2018_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)
plot_aspect_2019 <- raster::extract(x = aspect_2019_crop, #raster
                               y = plot_spdf, buffer = 10, # specify a 10 m radius
                               fun = mean, sp = TRUE, stringsAsFactors = FALSE)

# Merge slope/aspect into df
slope_aspect_17 <- car_coords %>% 
    left_join(trap_aspect_2017@data) %>% 
    rename(trap17aspect=layer) %>% 
    left_join(plot_aspect_2017@data) %>% 
    rename(plot17aspect=layer) %>% 
    left_join(trap_slope_2017@data) %>% 
    rename(trap17slope=layer) %>% 
    left_join(plot_slope_2017@data) %>% 
    rename(plot17slope=layer) %>% 
    distinct() %>%
    arrange(plotID)

# Save LAI at trap-level to incorporate
save(slope_aspect_17, file="data_derived/slope_aspect_17.Rdata")
