### Download DEM data and save the full raster ###
# downloads the full DEM for NIWO, saves a 50mx50m lower resolution raster
library(neonUtilities)
library(raster)
library(dplyr)


byFileAOP(dpID = "DP3.30024.001", site = "NIWO") #dem data is 1.5GB, takes a while (~40-60min for me)
# not sure if this all downloaded for me, as I only got 2017 data, but works for now

# get vector of all .tif file names
path_to_tifs = "../DP3.30024.001/2017/FullSite/D13/2017_NIWO_1/L3/DiscreteLidar/DTMGtif/"
file_list = list.files(path_to_tifs)
files = paste0(path_to_tifs, file_list)

# make list of rasters from each individual .tif
rast_list = list()
for (i in 1:length(files)){
  rast_list[[i]] = raster(files[i])
}

# merge all rasters
full = do.call(merge, rast_list)

full_50 = full %>%
  aggregate(fact = 50) #make it 50x50m


#saving a lower resolution raster for plotting
writeRaster(full_50, filename = "neonDEM_50x50.grd")




