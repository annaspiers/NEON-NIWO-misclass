# This script downloads mostly NEON and one Niwot LTER data product to be used as variables in carabid species distribution modeling.

library(dplyr)
library(neonUtilities)
library(geoNEON)
library(raster)

### Functions ###

# Removes columns with NA
remove_NA_cols <- function(list_of_lists){
    for (i in 1:length(list_of_lists))
        list_of_lists[[i]] <- list_of_lists[[i]][, colSums(is.na(list_of_lists[[i]])) != 
                                                     nrow(list_of_lists[[i]])] 
    return(list_of_lists)
}


### Load raw data ###

load("data_raw/carabids_NIWO.Rdata")
# load("data_raw/soil_wc_NIWO.Rdata")
load("data_raw/woody_veg_NIWO.Rdata")
load("data_raw/litter_woodfall_NIWO.Rdata")
load("data_raw/ir_bio_temp_NIWO.Rdata")
#load("data_raw/soil_temp_NIWO.Rdata")
load("data_raw/precip_NIWO.Rdata")
load("data_raw/precip_NIWO_LTER_saddle.Rdata") #precip_LTER_saddle
load("data_raw/precip_NIWO_LTER_C1.Rdata") #precip_LTER_C1
load("data_raw/rad_net_NIWO.Rdata")
load("data_raw/rad_short_dir_diff_NIWO.Rdata")
forest_niwot_all = read.csv("data_raw/PP_plot_data_1982-2016.tv.data.csv") #just the plot info


### Clean data ###

# Remove columns with NAs from every data package
all_data_packages <- c(carabid_abund, carabid_barcode, woody_veg, litter_woodfall,
                       ir_bio_temp, precip, rad_net, rad_short_dir_diff) #soil_temp, soil_wc
for (i in 1:length(all_data_packages)) {
    all_data_packages[i] <- remove_NA_cols(all_data_packages[i])
}

### Extract spatial data ###

# carabid abundance 
carabid_spat <- def.extr.geo.os(data = carabid_abund$bet_fielddata, 'namedLocation', locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude, Value.for.Plot.ID) %>%
  cbind('data_type' = rep('carabid'))

# # soil water content
# soil_water_spat = def.extr.geo.os(data = soil_wc$vst_perplotperyear,'namedLocation', locOnly=T) %>%
#   dplyr::select(api.decimalLatitude, api.decimalLongitude) %>%
#   cbind('data_type' = rep('soil_water'))

# woody veg
woody_spat = def.extr.geo.os(data = woody_veg$vst_perplotperyear, 'namedLocation', locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude, Value.for.Plot.ID) %>%
  cbind('data_type' = rep('woody_veg'))

# litter and woodfall
litter_spat = def.extr.geo.os(data = litter_woodfall$ltr_pertrap, 'namedLocation', locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude, Value.for.Plot.ID) %>%
  cbind('data_type' = rep('litter'))

# infrared biological temp
ir_bio_temp_spat = ir_bio_temp$sensor_positions_00005 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('Value.for.Plot.ID' = NA, 'data_type' = rep('ir_bio_temp')) %>% #site-level, not plot-level
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)

# # soil temperature
# soil_temp_spat = soil_temp$sensor_positions_00041 %>%
#   dplyr::select(referenceLatitude, referenceLongitude) %>%
#   cbind('data_type' = rep('soil_temp')) %>%
#   rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)

# NEON precipitation
precip_spat <- precip$sensor_positions_00006 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('Value.for.Plot.ID' = NA, 'data_type' = rep('precip')) %>% #site-level, not plot-level
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)

# NIWO LTER precipitation - see download script for data/code source
for (i in 1:2) {
  if (i == 1) {precip_obj <- precip_LTER_C1} else {precip_obj <- precip_LTER_saddle}

  if (class(precip_obj$LTER_site)!="factor") precip_obj$LTER_site<- as.factor(precip_obj$LTER_site)
  if (class(precip_obj$local_site)!="factor") precip_obj$local_site<- as.factor(precip_obj$local_site)
  tmpDateFormat<-"%Y-%m-%d"
  tmp1date<-as.Date(precip_obj$date,format=tmpDateFormat)
  if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){precip_obj$date <- tmp1date } else {print("Date conversion failed for precip_obj$date. Please inspect the data and do the date conversion yourself.")}
  rm(tmpDateFormat,tmp1date) 
  if (class(precip_obj$ppt_tot)=="factor") precip_obj$ppt_tot <-as.numeric(levels(precip_obj$ppt_tot))[as.integer(precip_obj$ppt_tot) ]
  if (class(precip_obj$qdays)=="factor") precip_obj$qdays <-as.numeric(levels(precip_obj$qdays))[as.integer(precip_obj$qdays) ]
  if (class(precip_obj$flag_ppt_tot)!="factor") precip_obj$flag_ppt_tot<- as.factor(precip_obj$flag_ppt_tot)
  precip_obj$LTER_site <- as.factor(ifelse((trimws(as.character(precip_obj$LTER_site))==trimws("NaN")),NA,as.character(precip_obj$LTER_site)))
  precip_obj$local_site <- as.factor(ifelse((trimws(as.character(precip_obj$local_site))==trimws("NaN")),NA,as.character(precip_obj$local_site)))
  precip_obj$ppt_tot <- ifelse((trimws(as.character(precip_obj$ppt_tot))==trimws("NaN")),NA,precip_obj$ppt_tot)
  precip_obj$qdays <- ifelse((trimws(as.character(precip_obj$qdays))==trimws("NaN")),NA,precip_obj$qdays)
  precip_obj$flag_ppt_tot <- as.factor(ifelse((trimws(as.character(precip_obj$flag_ppt_tot))==trimws("NaN")),NA,as.character(precip_obj$flag_ppt_tot)))
  precip_obj <- subset(precip_obj, date >= "2015-01-01") # Keep data since 2015
  
  if (i == 1) {precip_LTER_C1 <- precip_obj} else {precip_LTER_saddle <- precip_obj}
}

# net radiation (shortwave and longwave)
rad_net_spat <- rad_net$sensor_positions_00006 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('Value.for.Plot.ID' = NA, 'data_type' = rep('rad_net')) %>% #site-level, not plot-level
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)

# shortwave radiation - diffuse and direct
rad_short_dir_diff_spat <- rad_short_dir_diff$sensor_positions_00023 %>%
  dplyr::select(referenceLatitude, referenceLongitude) %>%
  cbind('Value.for.Plot.ID' = NA, 'data_type' = rep('rad_short_dir_diff')) %>% #site-level, not plot-level
  rename('api.decimalLatitude' = referenceLatitude,'api.decimalLongitude' = referenceLongitude)

# niwot lter forest plots
forest_niwot_spat = data.frame('api.decimalLatitude' = forest_niwot_all$lat, 
                               'api.decimalLongitude' = forest_niwot_all$long, 
                               'Value.for.Plot.ID' = NA, #not a NEON data product
                               'data_type' = rep('forest_niwot')) #uses same naming as 





# Join spatial data (lat/long)
all_spat <- rbind(carabid_spat, woody_spat, litter_spat, ir_bio_temp_spat, 
                  precip_spat, rad_net_spat, rad_short_dir_diff_spat, forest_niwot_spat) #soil_water_spat, soil_temp_spat
write.csv(x = all_spat, file = "data_derived/all_pts_spat.csv", row.names = FALSE)


### Save cleaned data do data_derived ###
save(carabid_abund, carabid_barcode, file="data_derived/carabids_NIWO.Rdata")
#save(soil_wc, file="data_derived/soil_wc_NIWO.Rdata")
save(woody_veg, file="data_derived/woody_veg_NIWO.Rdata")
save(litter_woodfall, file="data_derived/litter_woodfall_NIWO.Rdata")
save(ir_bio_temp, file="data_derived/ir_bio_temp_NIWO.Rdata")
#save(soil_temp, file="data_derived/soil_temp_NIWO.Rdata")
save(precip, file="data_derived/precip_NIWO.Rdata")
save(precip_LTER_saddle, file="data_derived/precip_NIWO_LTER_saddle.Rdata") 
save(precip_LTER_C1, file="data_derived/precip_NIWO_LTER_C1.Rdata") 
save(rad_net, file="data_derived/rad_net_NIWO.Rdata")
save(rad_short_dir_diff, file="data_derived/rad_short_dir_diff_NIWO.Rdata")
