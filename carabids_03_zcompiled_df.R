# This script compiles carabid abundance and richness data with variables
# describing layers of data aggregation as well as environmental variables. 

library(dplyr)
library(stringr) #word()
library(geoNEON)
library(ggplot2)
library(forcats) #fct_reorder

### Load data ####
load("data_derived/carabids_NIWO.Rdata")
list2env(carabid_abund, .GlobalEnv)
#load("data_derived/soil_wc_NIWO.Rdata")
#load("data_derived/woody_veg_NIWO.Rdata")
#load("data_derived/litter_woodfall_NIWO.Rdata")
#load("data_derived/ir_bio_temp_NIWO.Rdata")
#load("data_derived/soil_temp_NIWO.Rdata")
#load("data_derived/precip_NIWO.Rdata")
#load("data_derived/rad_net_NIWO.Rdata")
#load("data_derived/rad_short_dir_diff_NIWO.Rdata")

### Combine into one df ###
select_spp <- c("Amara alpina", "Amara quenseli", "Calathus advena", "Carabus taedatus", "Cymindis unicolor", "Harpalus nigritarsis", "Pterostichus restrictus")

# each row is individual beetle
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

# each row is a para_sciname in a sample with 0 abundance rows included
all_samples_df <- model_df_by_ind %>%
    group_by(siteID, plotID, trapID, sampleID, collectDate, col_year, col_month, col_day, dayofyear, trappingDays, decimalLatitude, decimalLongitude, elevation, nlcdClass, soilOrder, para_sciname) %>%
  summarize(sp_abund=n())

# add absences
all_samps_na_sp_df <- all_samples_df %>%
  mutate(para_sciname = NA) %>%
  dplyr::select(-sp_abund)
for(i in 1:length(select_spp)) {
    df_i <- all_samples_df %>%
      mutate(para_sciname = NA) %>%
      dplyr::select(-sp_abund)
    df_i$para_sciname <- select_spp[i]
    all_samps_na_sp_df <- rbind(all_samps_na_sp_df, df_i)
    rm(df_i)
}
# remove initial df chunk with para_sciname values as NA
all_samps_all_sp_df <- all_samps_na_sp_df %>% filter(is.na(para_sciname)==FALSE) 

# Join with df organized by collectDate  in trap in plot and summarize for sp_abund count
model_df_by_sample <- all_samps_all_sp_df %>% 
  left_join(all_samples_df) %>%
  arrange(collectDate)

# Convert sp_abund NAs to 0s
model_df_by_sample$sp_abund[is.na(model_df_by_sample$sp_abund)] <- 0
unique(model_df_by_sample$sp_abund)






### Save df ###
  
write.csv(model_df_by_ind, file="data_derived/model_df_by_individual_beetle.csv", row.names = FALSE) #df by individuals
write.csv(model_df_by_sample, file="data_derived/model_df_by_species_in_sample.csv", row.names = FALSE) #df by para_sciname in sample


