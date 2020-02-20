####### Gather point data for plotting ######
library(neonUtilities)
library(geoNEON)
library(dplyr)
library(raster)


##### Data Downloads #####
# these data products vary in size. Those with check.size=FALSE are small enough
# to be downloaded in less than 1 minute. Others take quite a lot of time.
# At the end, a csv is written with all of the spatial locations of each data
# set in this document.

### carabid data ###
carabids_all <- loadByProduct(dpID = "DP1.10022.001", 
                              site = "NIWO", 
                              check.size = FALSE) 
carabid_spat <- def.extr.geo.os(data = carabids_all$bet_fielddata, 
                                   'namedLocation', 
                                   locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude) %>%
  cbind('data_type' = rep('carabid'))


### soil water ###
soil_water_all = loadByProduct(dpID = "DP1.10098.001", 
                           site = "NIWO", 
                           check.size = FALSE)
soil_water_spat = def.extr.geo.os(data = soil_water_all$vst_perplotperyear, 
                                  'namedLocation', 
                                  locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude) %>%
  cbind('data_type' = rep('soil_water'))


### woody veg ###
woody_all = loadByProduct(dpID = "DP1.10098.001", 
                                 site = "NIWO", 
                                 check.size = FALSE)
list2env(woody_all, .GlobalEnv)
woody_spat = def.extr.geo.os(data = woody_all$vst_perplotperyear, 
                             'namedLocation', 
                             locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude) %>%
  cbind('data_type' = rep('woody_veg'))


### litter and woodfall ###
litter_all = loadByProduct(dpID = 'DP1.10033.001', 
                       site = "NIWO", 
                       check.size = FALSE)
litter_spat = def.extr.geo.os(data = litter_all$ltr_pertrap, 
                              'namedLocation', 
                              locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude) %>%
  cbind('data_type' = rep('litter'))



### niwot long term forest plots ###
forest_niwot_all = read.csv("data_raw/PP_plot_data_1982-2016.tv.data.csv") #just the plot info
forest_niwot_spat = data.frame('api.decimalLatitude' = forest_niwot_all$lat, 
                               'api.decimalLongitude' = forest_niwot_all$long, 
                               'data_type' = rep('forest_niwot')) #uses same naming as NEON data





ir_temp_all = loadByProduct(dpID = 'DP1.00005.001', site = 'NIWO', check.size = FALSE) # takes about 1-2 minute
ir_temp_spat = ir_temp_all$sensor_positions_00005 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('data_type' = rep('ir_temp')) %>%
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)


soil_temp_all = loadByProduct(dpID = 'DP1.00041.001', site = 'NIWO') # 860 MB
#save(soil_temp_all, file = '../soil_temp_all.RData')
soil_temp_spat = soil_temp_all$sensor_positions_00041 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('data_type' = rep('soil_temp')) %>%
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)


### Make single data frame of long lat points ###

all_pts = rbind(carabid_spat, litter_spat, woody_spat, soil_water_spat, forest_niwot_spat, ir_temp_spat, soil_temp_spat)
write.csv(x = all_pts, file = "data_raw/all_pts_spat.csv", row.names = FALSE)
