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

plot(chm_crop, legend.args = list(text = 'Canopy height (m)'))
points(car_coords$trap.Easting,car_coords$trap.Northing)
#saved via zoom window


slope = raster("data_derived/neonSlope_1x1_2017.grd") %>% crop(ext)

plot(slope, legend.args = list(text = 'Slope'))
points(car_coords$trap.Easting,car_coords$trap.Northing)


aspect = raster("data_derived/neonAspect_1x1_2017.grd") %>% crop(ext)

plot(aspect, legend.args = list(text = 'Aspect (degrees from N)'))
points(car_coords$trap.Easting,car_coords$trap.Northing)




