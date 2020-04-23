library(dplyr)
library(neonUtilities)
library(raster)


# Canopy height explo -----------------------------------------------------

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



# LAI  ------------------------------------------------------------

# Dwonload data
byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2017, check.size=T, savepath="data_raw/") # 708 MB
byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2018, check.size=T, savepath="data_raw/") # 795 MB
byFileAOP(dpID="DP3.30012.001", site="NIWO", year=2019, check.size=T, savepath="data_raw/") # 1800 MB
# AIS From NEON LAI data product details: AOP legacy data (those collected in 2013 through 2016) currently has partial availability, and will be completely available by April of 2019.


# Merge tiles

years <- c("2017","2018","2019")

for (i in years) {
    
    # get vector of all .tif file names
    path_to_tifs = paste0("data_raw/DP3.30012.001/",i,"/FullSite/D13/",i,"_NIWO_",which(years==i),"/L3/Spectrometer/LAI/")
    
    file_list <- list.files(path_to_tifs)
    file_paths <- paste0(path_to_tifs, file_list)
    
    # make list of rasters from each individual .tif
    rast_list = list()
    for (j in 1:length(file_paths)){
        rast_list[[j]] <- raster(file_paths[j])
    }
    
    # merge all rasters
    full_lai_1x1 <- do.call("merge", rast_list)
    writeRaster(full_lai_1x1, filename = paste0("data_derived/neonLAI_1x1_",i,".grd"))
}



# Slope/Aspect ------------------------------------------------------------
# Download data
years <- c("2017","2018","2019")
byFileAOP(dpID="DP3.30025.001",site="NIWO",year=years[1],check.size=T,savepath="data_raw/") # 1.4 GB
byFileAOP(dpID="DP3.30025.001",site="NIWO",year=years[2],check.size=T,savepath="data_raw/") # 1.5 GB
byFileAOP(dpID="DP3.30025.001",site="NIWO",year=years[3],check.size=T,savepath="data_raw/") # 1.6 GB

# Merge tiles
product <- c("Aspect", "Slope")
for (i in years) {
    for (j in product) {
    
        # get vector of all .tif file names
        path_to_tifs = paste0("data_raw/DP3.30025.001/",i,"/FullSite/D13/",i,"_NIWO_",which(years==i),"/L3/DiscreteLidar/",j,"Gtif/")

        file_list <- list.files(path_to_tifs)
        file_paths <- paste0(path_to_tifs, file_list)

        # make list of rasters from each individual .tif
        rast_list = list()
        for (k in 1:length(file_paths)){
            rast_list[[k]] <- raster(file_paths[k])
        }

        # merge all rasters
        full_1x1 <- do.call("merge", rast_list)
        writeRaster(full_1x1, filename = paste0("data_derived/neon",j,"_1x1_",i,".grd"))
    }
    print(which(years==i))
    if (which(years==i) == length(years)) cat("DONE!\n")
} #started 1:30

