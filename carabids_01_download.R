# This script downloads mostly NEON and one Niwot LTER data product to be used as variables in carabid species distribution modeling.

library(neonUtilities)

### Download data ###

# carabid abundance 
carabid_abund <- loadByProduct(dpID = "DP1.10022.001", site = "NIWO", check.size = FALSE)

# carabid barcode
carabid_barcode <- loadByProduct(dpID = "DP1.10020.001", site = "NIWO", check.size = FALSE) 

# soil water content
soil_wc <- loadByProduct("DP1.00094.001",  site = "NIWO") #476 MB, takes ?? to download and stack
#AIS stalls after stacking 1min swc

# woody veg
woody_veg <- loadByProduct(dpID = "DP1.10098.001", site = "NIWO", check.size = FALSE)

# litter and woodfall
litter_woodfall <- loadByProduct(dpID = 'DP1.10033.001', site = "NIWO", check.size = FALSE)

# infrared biological temp
ir_bio_temp <- loadByProduct(dpID = 'DP1.00005.001', site = 'NIWO', check.size = FALSE) # takes about 1-2 minute

# soil temperature
soil_temp <- loadByProduct(dpID = 'DP1.00041.001', site = 'NIWO') # 860 MB, took ?? min

# Precipitation
precip <- loadByProduct(dpID = 'DP1.00006.001', site = 'NIWO', check.size = FALSE) 

# Rad net
rad_net <- loadByProduct(dpID = 'DP1.00006.001', site = 'NIWO', check.size = FALSE) 

# Rad short direct diffuse
rad_short_dir_diff <- loadByProduct(dpID = 'DP1.00023.001', site = 'NIWO', check.size = FALSE) #78 MB, took a minute 

### Save to data_raw ###
save(carabid_abund, carabid_barcode, file="data_raw/carabids_NIWO.Rdata")
#save(soil_wc, file="data_raw/soil_wc_NIWO.Rdata")
save(woody_veg, file="data_raw/woody_veg_NIWO.Rdata")
save(litter_woodfall, file="data_raw/litter_woodfall_NIWO.Rdata")
save(ir_bio_temp, file="data_raw/ir_bio_temp_NIWO.Rdata")
#save(soil_temp, file="data_raw/soil_temp_NIWO.Rdata")
save(precip, file="data_raw/precip_NIWO.Rdata")
save(rad_net, file="data_raw/rad_net_NIWO.Rdata")
save(rad_short_dir_diff, file="data_raw/rad_short_dir_diff_NIWO.Rdata")

