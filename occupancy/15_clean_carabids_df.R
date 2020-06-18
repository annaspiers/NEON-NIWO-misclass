

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

# I. Create df where each row is an individual beetle ------------------------

model_df_by_ind <- bet_parataxonomistID %>%
  dplyr::select(individualID, siteID, plotID, trapID, collectDate, taxonRank, morphospeciesID, "para_sciname" = scientificName) %>% #select desired columns
  mutate(para_sciname = ifelse(para_sciname == "Carabidae sp.", "Carabidae spp.", para_sciname)) %>%
  #roll similar parataxonomist IDs into one: Carabidae sp. (1) and Carabidae spp. (215) for 2015-2018
  mutate(para_sciname = ifelse(para_sciname == "Amara quenseli quenseli", "Amara quenseli", para_sciname)) %>%
  #change para_sciname subspecies ID to species-level
  mutate(para_morph = ifelse(is.na(morphospeciesID), #combine parataxonomist ID and morphospecies
                             as.character(para_sciname), 
                             paste0(as.character(para_sciname),"/", morphospeciesID,sep=""))) %>%
  mutate(col_year = lubridate::year(collectDate), 
         col_month = lubridate::month(collectDate), 
         col_day = lubridate::day(collectDate),
         dayofyear = as.numeric(strftime(collectDate, format = "%j"))) %>%
  left_join(distinct(bet_expertTaxonomistIDProcessed %>% 
                       mutate(expert_sciname = paste(genus, specificEpithet)) %>%
                       dplyr::select(individualID, expert_sciname, taxonID))) 

model_df_by_ind <- model_df_by_ind %>%
  left_join(model_df_by_ind %>% 
              distinct(col_year, collectDate) %>%
              group_by(col_year) %>%
              mutate(col_index = 1:n())) #index the collection date within each year

# # Check for parataxonomist IDs that are at the subspecies level 
# model_df_by_ind %>% select(para_sciname) %>%
#    filter(sapply(strsplit(as.character(para_sciname), " "), length) > 2)
# # change para_sciname "Amara quenseli quenseli" to "Amara quenseli" in model_df_by_ind above

# II. Create df where each row is a species in a trap from a sing --------

# Initialize df organized by each species by collectionDate by trap by plot
model_df_by_sample <- model_df_by_ind %>%
  distinct(plotID, trapID, collectDate, para_morph) %>%
  arrange(plotID, trapID, collectDate, para_morph) %>%
  complete(plotID, trapID, collectDate, para_morph) %>%
  left_join(model_df_by_ind %>% 
              group_by(plotID,trapID,collectDate,para_morph) %>%
              summarize(sp_abund = n()) ) %>% # add in non-zero species abundance
  mutate(sp_abund = ifelse(is.na(sp_abund),0,sp_abund)) %>% # turn NA abundance into 0
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


# Create df for samples but only with data 2015-2018 ----------------------

model_df_by_sample_pre2019 <- model_df_by_ind %>%
  #mutate(collectDate = as.Date(collectDate)) %>%
  #filter(collectDate <= as.Date("2018-12-31")) %>%
  filter(collectDate <= "2018-12-31 GMT") %>%
  distinct(plotID, trapID, collectDate, para_morph) %>%
  arrange(plotID, trapID, collectDate, para_morph) %>%
  complete(plotID, trapID, collectDate, para_morph) %>%
  left_join(model_df_by_ind %>% 
              group_by(plotID,trapID,collectDate,para_morph) %>%
              summarize(sp_abund = n()) ) %>% # add in non-zero species abundance
  mutate(sp_abund = ifelse(is.na(sp_abund),0,sp_abund)) %>% # turn NA abundance into 0
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


# Save dataframes to csv --------------------------------------------------

write.csv(model_df_by_ind, file="occupancy/model_df_by_individual_beetle.csv") 
write.csv(model_df_by_sample, file="occupancy/model_df_by_species_in_sample_alldates.csv") 
write.csv(model_df_by_sample_pre2019, file="occupancy/model_df_by_species_in_sample_2015-2018.csv") #df by 
