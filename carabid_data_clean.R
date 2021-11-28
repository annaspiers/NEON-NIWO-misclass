# This script downloads and cleans raw NEON carabid data. Here we compile
# carabid abundance and richness data with variables describing layers of data
# aggregation.

library(neonUtilities)
library(dplyr)
library(stringr) #word()
library(ggplot2)
library(forcats) #fct_reorder
library(tidyr) #separate()


# Load data ---------------------------------------------------------------

# Download data from NEON
if (!file.exists("data/carabids_NIWO_raw.rds")) {
  carabid_abund <- loadByProduct(dpID = "DP1.10022.001", site = "NIWO", check.size = F)
  saveRDS(carabid_abund, file="data/carabids_NIWO_raw.rds")
}

# Load data
carabid_abund <- readRDS("data/carabids_NIWO_raw.rds")
list2env(carabid_abund, .GlobalEnv)
rm(carabid_abund)


# Generate dataframes to be used in dynamic occupancy model ---------------

all_paratax_df <- bet_sorting %>%
  filter(sampleType == "carabid",
         # use data through 2019 since this is the last season that the expert taxonomist has ID'ed 
         collectDate <= "2019-12-31 GMT") %>%
  mutate(#change subspecies to species-level
         scientificName = ifelse(scientificName == "Amara quenseli quenseli", 
                                    "Amara quenseli", scientificName),
         #roll similar parataxonomist IDs into one: Carabidae sp. (1) and Carabidae spp. (215) for 2015-2018
         scientificName = ifelse(scientificName == "Carabidae sp.", "Carabidae spp.", scientificName),
         #create variable that combines parataxonomist ID and morphospecies ID (when present)
         scimorph_combo = ifelse(is.na(morphospeciesID), 
                                 scientificName, 
                                 paste0(scientificName,"/",morphospeciesID,sep="")),
         col_year = lubridate::year(collectDate))
all_paratax_df <- all_paratax_df %>%  #index the collection date within each year
  left_join(all_paratax_df %>% 
              distinct(col_year, collectDate) %>%
              group_by(col_year) %>%
              mutate(col_index = 1:n()))

pinned_df <- bet_parataxonomistID %>%
  filter( collectDate <= "2019-12-31 GMT") %>%
  mutate(#change subspecies to species-level
    scientificName = ifelse(scientificName == "Amara quenseli quenseli", 
                            "Amara quenseli", scientificName),
    #roll similar parataxonomist IDs into one: Carabidae sp. (1) and Carabidae spp. (215) for 2015-2018
    scientificName = ifelse(scientificName == "Carabidae sp.", "Carabidae spp.", scientificName),
    #create variable that combines parataxonomist ID and morphospecies ID (when present)
    scimorph_combo = ifelse(is.na(morphospeciesID), 
                            scientificName, 
                            paste0(scientificName,"/",morphospeciesID,sep="")),
    col_year = lubridate::year(collectDate))
pinned_df <- pinned_df %>%  #index the collection date within each year
  left_join(pinned_df %>% 
              distinct(col_year, collectDate) %>%
              group_by(col_year) %>%
              mutate(col_index = 1:n()))

expert_df <- bet_expertTaxonomistIDProcessed %>%
  filter(collectDate <= "2019-12-31 GMT") %>%
  mutate(scientificName = ifelse(taxonRank=="subspecies",
                                    paste0(genus," ",specificEpithet,sep=""), scientificName),
         col_year = lubridate::year(collectDate)) 
expert_df <- expert_df %>% #index the collection date within each year
  left_join(expert_df %>% 
              distinct(col_year, collectDate) %>%
              group_by(col_year) %>%
              mutate(col_index = 1:n()))


# Save dataframes locally -------------------------------------------------
saveRDS(all_paratax_df, "data/all_paratax_df.rds")
saveRDS(pinned_df, "data/pinned_df.rds")
saveRDS(expert_df, "data/expert_df.rds")
