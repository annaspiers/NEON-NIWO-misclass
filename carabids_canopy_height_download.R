library(dplyr)
library(neonUtilities)
library(raster)


###### Canopy height explo ######
byFileAOP(dpID = "DP3.30015.001", site = "NIWO", check.size = F) # 830 MB #called ecosystem structure in data portal

#only 2017 downloaded for some reason - will need to figure this out at some point
path_to_tifs = "data_raw/DP3.30015.001/2017/FullSite/D13/2017_NIWO_1/L3/DiscreteLidar/CanopyHeightModelGtif/"
file_list = list.files(path_to_tifs)
files = paste0(path_to_tifs, file_list)

# make list of rasters from each individual .tif
rast_list = list()
for (i in 1:length(files)){
    rast_list[[i]] = raster(files[i])
}

# merge all rasters
full = do.call(merge, rast_list)
writeRaster(full, filename = "data_derived/neonCHM_1x1.grd")

full_25 = full %>%
    aggregate(fact = 25) #make it 25x25m


#saving a lower resolution raster for plotting
writeRaster(full_25, filename = "data_derived/neonCHM_25x25.grd")


