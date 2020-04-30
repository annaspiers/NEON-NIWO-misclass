# Here, we explore the NEON LAI product at the Niwot Ridge site and summarize LAI at the trap-level.
# resources: 
    # https://www.earthdatascience.org/courses/earth-analytics/remote-sensing-uncertainty/extract-data-from-raster/
    # https://www.neonscience.org/dc-plot-raster-data-r 
    # https://www.earthdatascience.org/courses/earth-analytics/remote-sensing-uncertainty/extract-data-from-raster/

library(dplyr)
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)


# Load data ---------------------------------------------------------------

LAI_raster_2017 = raster("data_derived/neonLAI_1x1_2017.grd")
LAI_raster_2018 = raster("data_derived/neonLAI_1x1_2018.grd")
LAI_raster_2019 = raster("data_derived/neonLAI_1x1_2019.grd")

load("data_derived/car_coords.RData")
LAI_proj = projection(LAI_raster_2019) #the same across years
trap_spdf <- SpatialPointsDataFrame(coords = car_coords %>%
                                        dplyr::select(trap.Easting, trap.Northing), 
                                    data=car_coords %>%
                                        dplyr::select(trap.Easting, trap.Northing), 
                                    proj4string = CRS(LAI_proj))
plot_spdf <- SpatialPointsDataFrame(coords = car_coords %>%
                                        dplyr::select(plot.Easting, plot.Northing), 
                                    data=car_coords %>%
                                        dplyr::select(plot.Easting, plot.Northing), 
                                    proj4string = CRS(LAI_proj))


# Crop LAI raster ---------------------------------------------------------

buff <- 500 #set arbitrary buffer - I think this is in m per the earthdatascience resource
car_extent <- extent(min(car_coords$trap.Easting) - buff, #xmin
                     max(car_coords$trap.Easting) + buff, #xmax
                     min(car_coords$trap.Northing) - buff, #ymin
                     max(car_coords$trap.Northing) + buff) #ymax
LAI_2017_crop <- crop(LAI_raster_2017, car_extent)
LAI_2018_crop <- crop(LAI_raster_2018, car_extent)
LAI_2019_crop <- crop(LAI_raster_2019, car_extent)
# AIS need to save these to data_derived

plot(LAI_2017_crop, main="2017 LAI")
points(car_coords$trap.Easting,car_coords$trap.Northing, pch=20, cex=0.5)
plot(LAI_2018_crop, main="2018 LAI")
plot(LAI_2019_crop, main="2019 LAI")


# What's the difference in LAI between years? -----------------------------

#save a lower resolution raster for plotting
LAI_agg_2017 <- LAI_2017_crop %>% aggregate(fact = 10) #make it 25x25m
LAI_agg_2018 <- LAI_2018_crop %>% aggregate(fact = 10)
LAI_agg_2019 <- LAI_2019_crop %>% aggregate(fact = 10)

hist(LAI_agg_2017, main="2017 LAI")
hist(LAI_agg_2018, main="2018 LAI")
hist(LAI_agg_2019, main="2019 LAI")
# 2018 and 2019 don't look right

plot(LAI_2017_crop, main="2017 LAI",
     breaks=c(0,.5,1,1.5,2,2.5,3),
     col=terrain.colors(3))
plot(LAI_2018_crop, main="2018 LAI",
     breaks=c(0,.5,1,1.5,2,2.5,3),
     col=terrain.colors(3))
plot(LAI_2019_crop, main="2019 LAI",
     breaks=c(0,.5,1,1.5,2,2.5,3),
     col=terrain.colors(3))
# AIS add plot coords too and label so we can visualize the names of the plots

# AIS contacted NEON to ask about the difference in LAI values between years


# Create LAI value around traps -------------------------------------------

# Plots are 40m x 40m. Traps lay within the plot at the middle of each edge.  If they laid exactly on plot edge, they would be about 28m from the traps on the neighboring sides. How large should we summarize LAI for each trap? As to not overlap with neighboring traps' buffers, the buffer should be at most 14m. Max, min, avg LAI? for now choose avg

#AIS experiment with different buffer sizes. 
# Mean ummarized to trap-level
trap_LAI_avg_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_avg_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_avg_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)

# Max summarized to trap-level
trap_LAI_max_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_max_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_max_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)

# Stand dev summarized to trap-level
trap_LAI_sd_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_sd_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
trap_LAI_sd_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = trap_spdf,
                               buffer = 10, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)

# Mean summarized to plot-level
plot_LAI_avg_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_avg_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_avg_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
# AIS something seems wonky with 2019 LAI, values are too low

# Max summarized to plot-level
plot_LAI_max_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_max_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_max_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = max, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)

# Max summarized to plot-level
plot_LAI_sd_2017 <- raster::extract(x = LAI_2017_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_sd_2018 <- raster::extract(x = LAI_2018_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_LAI_sd_2019 <- raster::extract(x = LAI_2019_crop, #raster
                               y = plot_spdf,
                               buffer = 30, # specify a 10 m radius
                               fun = sd, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)

# AIS while figuring out what's up with 2019...
# Just take the 2017 average values for trap and plot
plot_trap_LAI_2017

# Save LAI at trap-level to incorporate
save(plot_trap_LAI_2017, file="data_derived/plot_trap_LAI_2017.Rdata")
