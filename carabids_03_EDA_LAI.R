# Here, we download the NEON LAI product at the Niwot Ridge site, explore the data, and try to summarize LAI at the trap-level.

library(dplyr)
library(neonUtilities)
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)

# Download LAI  -----------------------------------------------------------

byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2017, check.size=T, savepath="data_raw/") # 708 MB
byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2018, check.size=T, savepath="data_raw/") # 795 MB
byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2019, check.size=T, savepath="data_raw/") # 1800 MB
# From NEON LAI data product details: AOP legacy data (those collected in 2013 through 2016) currently has partial availability, and will be completely available by April of 2019.


# Merge tiles -------------------------------------------------------------

# get vector of all .tif file names
path_to_tifs = "data_raw/DP3.30012.001/2017/FullSite/D13/2017_NIWO_1/L3/Spectrometer/LAI/"
file_list = list.files(path_to_tifs)
file_paths = paste0(path_to_tifs, file_list)

# make list of rasters from each individual .tif
rast_list = list()
for (i in 1:length(file_paths)){
    rast_list[[i]] = raster(file_paths[i])
}

# merge all rasters
full_lai_1x1 = do.call(merge, rast_list)
writeRaster(full_lai_1x1, filename = "data_derived/neonLAI_1x1.grd")


# Crop LAI raster ---------------------------------------------------------

# Convert trap coordinates (car_coords) to shapefile
# AIS only 3 traps from plot 13, get the last one (N)
# Change column titles so online converter to shp can interpret csv: https://mygeodata.cloud/converter/csv-to-shp - limit of 3 free datasets per month
# Another source, but didn't use: https://www.gisconvert.com/
# AIS use this more reproducible and automated method:  https://www.neonscience.org/dc-csv-to-shapefile-r
load("data_derived/car_coords.Rdata")
# Make each plot and plot_trap a different row with their Northings/Eastings
#AIS probably a more automated way to do this, but good for now
# AIS incorporate plot coord uncertainty into this. for now, drop trap.coord_uncertainty
trap_coords <- car_coords %>%
    dplyr::select(-c(plot.Easting, plot.Northing, trap.coord_uncertainty), 
                  Easting=trap.Easting, Northing=trap.Northing)
plot_coords <- car_coords %>%
    dplyr::select(-c(trap.Easting, trap.Northing, trap.coord_uncertainty), 
                  Easting=plot.Easting, Northing=plot.Northing) %>%
    mutate(trapID = NA) %>%
    distinct()
plot_trap_coords <- rbind(trap_coords, plot_coords)
write.csv(plot_trap_coords, file = "data_derived/plot_trap_coords.csv")
# assuming CRS NAD83 / UTM zone 13N (EPSG:26913)
# Visualize new shp
niwo_trap_coords <- readOGR("data_derived/trap_coord_geo/car_coords_xy-point.shp")

# Crop to extent of carabid traps
buff <- 300 #set arbitrary buffer
mod_site_extent <- extent(extent(niwo_trap_coords)[1]-buff, #xmin
                     extent(niwo_trap_coords)[2]+buff, #xmax
                     extent(niwo_trap_coords)[3]-buff, #ymin
                     extent(niwo_trap_coords)[4]+buff) #ymax
full_1x1_crop <- crop(full_1x1, mod_site_extent)
plot(full_1x1_crop, 
     main = "LAI cropped to traps with buffer",
     col = gray.colors(80, start = .9, end = .3))
plot(niwo_trap_coords, 
     pch = ".",
     col=2,
     add = TRUE)
#AIS add plot coords too and label so we can visualize the names of the plots
#AIS change LAI color spectrum no monochrome gradient


#save a lower resolution raster for plotting
full_25x25 = full_1x1 %>%
    aggregate(fact = 25) #make it 25x25m
writeRaster(full_25, filename = "data_derived/neonLAI_25x25.grd")

# Create buffer around traps and grab LAI value
# resource: https://www.earthdatascience.org/courses/earth-analytics/remote-sensing-uncertainty/extract-data-from-raster/
# Plots are 40m x 40m. Traps lay within the plot at the middle of each edge.  If they laid exactly on plot edge, they would be about 28m from the traps on the neighboring sides. How large should we summarize LAI for each trap? As to not overlap with enighboring traps' buffers, the buffer should be at most 14m. For now, max, min, avg LAI? for now choose avg
#AIS experiment with different buffer sizes. 
#AIS summarize to trap and plot-levels
trap_LAI_sp <- raster::extract(full_1x1_crop,
                    niwo_trap_coords,
                    buffer = 10, # specify a 10 m radius
                    fun = mean, # extract the MEAN value from each trap
                    sp = TRUE, # create spatial object
                    stringsAsFactors = FALSE)
# Rename layer column
trap_LAI_sp@data <- trap_LAI_sp@data %>%
    dplyr::select(-X, trap_LAI=layer)
trap_LAI_sp@data
min(trap_LAI_sp@data$trap_LAI)
max(trap_LAI_sp@data$trap_LAI)

# Save LAI at trap-level to incorporate
save(trap_LAI_sp, file="data_derived/trap_LAI_sp.Rdata")
