library(dplyr)
library(ggplot2)
library(raster)
library(ggthemes)

library(rgeos)

load("data_derived/car_coords.RData")
#canopy = raster("data_derived/neonCHM_25x25.grd")
canopy = raster("data_derived/neonCHM_1x1.grd") #full resolution
chm_proj = projection(canopy)

bord = 500 #number of pixels (or meters?) for border
xmin = min(car_coords$trap.Easting) - bord
xmax = max(car_coords$trap.Easting) + bord
ymin = min(car_coords$trap.Northing) - bord
ymax = max(car_coords$trap.Northing) + bord

ext = extent(xmin, xmax, ymin, ymax)
chm_crop = crop(canopy, ext)

plot(chm_crop)
points(car_coords$trap.Easting,car_coords$trap.Northing)


#plot level buffer
#get coords and calc min distance between plots
plot_coords = car_coords %>% dplyr::select(plot.Easting, plot.Northing, plotID) %>% unique() #unique plots
min_dist = min(dist(cbind(plot_coords$plot.Easting, plot_coords$plot.Northing), method = "euclidean")) #minimum distance between plots

#make spatial points dataframe for use with gBuffer
plot_spat = SpatialPointsDataFrame(cbind(plot_coords$plot.Easting, plot_coords$plot.Northing),
                                   data = plot_coords, proj4string = CRS(chm_proj))


plt_buffer_width = 100
plot_CHM_sp <- raster::extract(chm_crop,
                               plot_spat,
                               buffer = plt_buffer_width, # specify a 10 m radius
                               fun = mean, # extract the MEAN value from each trap
                               sp = TRUE, # create spatial object
                               stringsAsFactors = FALSE)
plot_CHM_sp@data = plot_CHM_sp@data %>%
    rename(plot_CH=layer)

save(plot_CHM_sp, file="data_derived/plot_CHM_sp.Rdata")



#trap level buffer
trap_coords = car_coords %>% dplyr::select(trap.Easting, trap.Northing, plotID, trapID) %>% unique() #unique plots
min_dist = min(dist(cbind(trap_coords$trap.Easting, trap_coords$trap.Northing), method = "euclidean")) #minimum distance between plots

#make spatial points dataframe 
trap_spat = SpatialPointsDataFrame(cbind(trap_coords$trap.Easting, trap_coords$trap.Northing),
                                   data = trap_coords, proj4string = CRS(chm_proj))


trap_buffer_width = 10
trap_CHM_sp <- raster::extract(chm_crop,
                                trap_spat,
                                buffer = trap_buffer_width, # specify a 10 m radius
                                fun = mean, # extract the MEAN value from each trap
                                sp = TRUE, # create spatial object
                                stringsAsFactors = FALSE)

trap_CHM_sp@data = trap_CHM_sp@data %>%
    rename(trap_CH=layer)

save(trap_CHM_sp, file="data_derived/trap_CHM_sp.Rdata")


