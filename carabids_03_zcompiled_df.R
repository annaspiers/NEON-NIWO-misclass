# This script compiles carabid abundance and richness data with variables
# describing layers of data aggregation as well as environmental variables. 

library(dplyr)
library(stringr) #word()
library(geoNEON)
library(ggplot2)
library(forcats) #fct_reorder
library(tidyr) #separate()

### Load data ####
load("data_derived/carabids_NIWO.Rdata")
list2env(carabid_abund, .GlobalEnv)
#load("data_derived/soil_wc_NIWO.Rdata")
#load("data_derived/woody_veg_NIWO.Rdata")
#load("data_derived/litter_woodfall_NIWO.Rdata")
#load("data_derived/ir_bio_temp_NIWO.Rdata")
#load("data_derived/soil_temp_NIWO.Rdata")
load("data_derived/summarized_precip.Rdata") 
#load("data_derived/rad_net_NIWO.Rdata")
#load("data_derived/rad_short_dir_diff_NIWO.Rdata")
load("data_derived/trap_LAI_1718avg.Rdata")
load("data_derived/trap_CHM_sp.Rdata")
load("data_derived/slope_aspect_17.Rdata")
load("data_derived/summarized_C1temp_GDD.Rdata")

# Select 7 most abundant spp.
select_spp <- c("Amara alpina", "Amara quenseli", "Calathus advena", "Carabus taedatus", "Cymindis unicolor", "Harpalus nigritarsis", "Pterostichus restrictus")


# I. Create df where each row is an individual beetle ------------------------

model_df_by_ind <- bet_parataxonomistID %>%
    filter(scientificName %in% select_spp) %>% #filter to selected 7 species
    dplyr::select(individualID, siteID, plotID, trapID, collectDate, taxonRank, morphospeciesID, "para_sciname" = scientificName) %>% #select desired columns
    mutate(col_year = lubridate::year(collectDate), 
           col_month = lubridate::month(collectDate), 
           col_day = lubridate::day(collectDate),
           dayofyear = as.numeric(strftime(collectDate, format = "%j"))) %>%
    left_join(distinct(bet_expertTaxonomistIDProcessed %>% 
            mutate(expert_sciname = paste(genus, specificEpithet)) %>%
            dplyr::select(individualID, expert_sciname, taxonID, sex))) %>%
    left_join(distinct(bet_fielddata %>% #join carabid field data columns: nlcdClass, lat, long, elev
            dplyr::select(sampleID, plotID, trapID, collectDate, trappingDays, decimalLatitude, decimalLongitude, elevation, nlcdClass))) %>%
    left_join(def.extr.geo.os(data = carabid_abund$bet_fielddata, 'namedLocation', locOnly=T) %>% #join plot soil type
            dplyr::select("plotID" = Value.for.Plot.ID, "soilOrder" = api.soilTypeOrder))



# II. Create df where each row is a species in a trap from a sing --------

# Initialize df organized by each species by collectionDate by trap by plot
model_df_by_sample <- model_df_by_ind %>%
  distinct(plotID, trapID, collectDate, para_sciname) %>%
  arrange(plotID, trapID, collectDate, para_sciname) %>%
  complete(plotID, trapID, collectDate, para_sciname) %>%
  left_join(model_df_by_ind %>% 
              group_by(plotID,trapID,collectDate,para_sciname) %>%
              summarize(sp_abund = n()) ) %>% # add in non-zero species abundance
  mutate(sp_abund = ifelse(is.na(sp_abund),0,sp_abund)) %>% # turn NA abundance into 0
  left_join(model_df_by_ind %>% 
              dplyr::select(siteID,plotID,plot_lat=decimalLatitude,plot_long=decimalLongitude,elev=elevation,nlcdClass,soilOrder) %>%
              distinct()) %>% # attach useful spatial variables from model_df_by_ind
  left_join(model_df_by_ind %>% 
              dplyr::select(collectDate,col_year,col_month,col_day,DOY=dayofyear) %>%
              distinct()) %>% # attach useful temporal variables from model_df_by_ind
  left_join(model_df_by_ind %>%
           distinct(col_year, collectDate) %>%
           group_by(col_year) %>%
           mutate(col_index = 1:n())) %>% #index the collection date within each year
  filter(!(plotID == "NIWO_004" & col_year == 2018), # Remove rows for plot 4 in 2018
         !(plotID == "NIWO_013" & (col_year==2015 | col_year==2016 | col_year==2017 ))) %>% # Remove rows for plot 13 in 2015-2017
  mutate(plot_trap = as.factor(paste(plotID, trapID, sep="")),
           occ = ifelse(sp_abund>0,1,0)) 

  
# III. Add in more predictor variables -----------------------------------------

# Add predictors to model df
model_df_by_sample <- model_df_by_sample %>%
  left_join(trap_CHM_sp@data %>%
              rename(trap_CHM = trap_CH)) %>% #add CHM
  left_join(trap_LAI_1718avg %>%
              dplyr::select(trap.Easting, trap.Northing, LAI_1718avg = avg1718)) %>%
  left_join(summ_precip %>%
              dplyr::select(-c(local_site,collectDate))) %>%
  left_join(slope_aspect_17) 

# Not quite the right syntax, will fix (in progress)
# %>%
#     left_join(summarized_C1temp_GDD %>% 
#                   dplyr::select(c(airtemp_avg, GDD_cum)))
#     

# Save dataframes to csv --------------------------------------------------
  
write.csv(model_df_by_ind, file="data_derived/model_df_by_individual_beetle.csv") #df by individuals
write.csv(model_df_by_sample, file="data_derived/model_df_by_species_in_sample.csv") #df by para_sciname in sample


