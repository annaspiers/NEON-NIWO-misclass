####### Plotting data in space #########
library(neonUtilities)
library(geoNEON)
library(dplyr)
library(ggplot2)
library(raster)
library(ggthemes)



# load neon DEM
dem = raster("data_derived/neonDEM_50x50.grd")
dem_crs = crs(dem)


# gather point data for carabids and predictors, comes in single data frame
#source("point_data.R") #script for downloading the spatial data
all_pts_spat = read.csv("data_derived/all_pts_spat.csv")


# project data frame to DEM projection
spat = SpatialPoints(all_pts_spat[,2:1], proj4string = CRS("+proj=longlat")) %>%
  spTransform(dem_crs) 
all_pts_df = data.frame('Easting' = spat@coords[,1], 'Northing' = spat@coords[,2]) %>%
  cbind('data_type' = all_pts_spat$data_type)



# plot full extent
dem_df = dem %>% as.data.frame(xy = TRUE)
all_pts_df = all_pts_df %>%
  mutate(carabid = case_when(data_type == 'carabid' ~ 1,
                             data_type != 'carabid' ~ 0))# make a binary column so that carabids can be plotted bigger?

ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = layer)) + #raster
  geom_point(data = all_pts_df %>% filter(!data_type %in% c('ir_temp','soil_temp')), 
             aes(x = Easting, y = Northing, 
                 colour = data_type, 
                 shape = data_type, 
                 size = as.factor(carabid))) +
  labs(fill = "Elevation (m)") +
  scale_fill_gradientn(colours = c('#222222', 'darkgreen', 'white')) +
  theme_map() +
  scale_color_manual(values = c('yellow', 'blue', rep('black',6))) + #line depends on number of data_types
  scale_shape_discrete(solid = FALSE) +
  guides(size = FALSE, fill = FALSE) 

#ggsave(file = "../map.png", width = 5, height = 5, dpi = 'retina')

# crop DEM to carabid traps
bord = 300 #number of pixels (or meters?) for border
xmin = min(all_pts_df$Easting) - bord
xmax = max(all_pts_df$Easting) + bord
ymin = min(all_pts_df$Northing) - bord
ymax = max(all_pts_df$Northing) + bord

ext = extent(xmin, xmax, ymin, ymax)
dem_crop = crop(dem, ext)
dem_crop_df = dem_crop %>% as.data.frame(xy=TRUE)

ggplot() +
  geom_raster(data = dem_crop_df, aes(x = x, y = y, fill = layer)) + #raster
  geom_point(data = all_pts_df %>% filter(!data_type %in% c('ir_bio_temp','precip')), aes(x = Easting, y = Northing, colour = data_type, shape = data_type, size = as.factor(carabid))) +
  labs(fill = "Elevation (m)") +
  scale_fill_gradientn(colours = c('#222222', 'darkgreen', 'white')) +
  theme_map() +
  scale_color_manual(values = c('yellow', 'blue', rep('black',6))) + #line depends on number of data_types
  scale_shape_discrete(solid = FALSE) +
  guides(size = FALSE, fill = FALSE) 

#ggsave(file = "../crop_map.png", width = 5, height = 5, dpi = 'retina')

# Plot carabid traps in proximity to Niwot Ridge C1 and saddle weather stations. 
# Resources for C1 and saddle coordinates:
  # https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.416.10
  # https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.401.3
  # convert lat/long to easting/northing: https://www.latlong.net/lat-long-utm.html
weather_st_sp <- data.frame(rbind(c("c1",449737.62,4434020.06),
                                  c("saddle",453573.88, 4431887.38)))
colnames(weather_st_sp) <- c("weather_st","Easting", "Northing") 

ggplot() + 
  geom_raster(data = dem_crop_df, aes(x = x, y = y, fill = layer)) + 
  geom_point(data = all_pts_df %>% filter(data_type %in% c('carabid')), aes(x = Easting, y = Northing)) +
  geom_point(data=weather_st_sp, aes(x = as.numeric(as.character(Easting)), y = as.numeric(as.character(Northing)), label=weather_st), color="red") +
  labs(fill = "Elevation (m)") +
  scale_fill_gradientn(colours = c('#222222', 'darkgreen', 'white')) +
  theme_map()
