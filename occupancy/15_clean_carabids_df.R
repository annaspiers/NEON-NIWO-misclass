# This script compiles carabid abundance and richness data with variables
# describing layers of data aggregation as well as environmental variables. 

library(dplyr)
library(stringr) #word()
library(ggplot2)
library(forcats) #fct_reorder
library(tidyr) #separate()

### Load data ####
load("data_derived/carabids_NIWO.Rdata")
list2env(carabid_abund, .GlobalEnv)
rm(carabid_abund, carabid_barcode)


# EDA ---------------------------------------------------------------------

# Do any scientific names conflict between bet_parataxonomist and bet_sorting?
bet_parataxonomistID %>%
  select(subsampleID, scientificname_pin = scientificName) %>%
  left_join(bet_sorting %>%
              select(subsampleID, scientificname_sort = scientificName)) %>%
  filter(scientificname_pin != scientificname_sort)
# No conflicts

# Do any morphospecies conflict between bet_parataxonomist and bet_sorting?
bet_parataxonomistID %>%
  select(subsampleID, morpho_pin = morphospeciesID) %>%
  left_join(bet_sorting %>%
              select(subsampleID, morpho_sort = morphospeciesID)) %>%
  filter(morpho_pin != morpho_sort)
# No conflicts


# In bet_sorting, carabids in a subsample have the same parataxonomist scientific name. 
# Does this hold in the pinned dataset, bet_parataxonomistID?
bet_parataxonomistID %>%
  count(subsampleID, sort=T) %>%
  filter(n>1)
# 400 rows. If I sort by scientific name, too, then there will be 400 rows
# again if carabids from the same subsample have the same scientific name
bet_parataxonomistID %>%
  count(subsampleID, scientificName, sort=T) %>%
  filter(n>1)
# Yes! nice. Let's visualize the scientific names for one subsample, as a sanity check
bet_parataxonomistID %>%
  filter(subsampleID=="/HCA//l5a79rKNXuRXzZyOJiS61iFZFxAPsVAzTxuGc=") %>%
  select(scientificName)


# Generate dataframes to be used in dynamic occupancy model ---------------

all_paratax_df <- bet_sorting %>%
  filter(sampleType == "carabid")  %>%
  mutate(scimorph_combo = ifelse(is.na(morphospeciesID), 
                                 scientificName, 
                                 paste0(scientificName,"/",morphospeciesID,sep="")))
pinned_df   <- bet_parataxonomistID
expert_df   <- bet_expertTaxonomistIDProcessed

# Check for scientific names that are at the subspecies level 
all_paratax_df %>% 
  #filter(sapply(strsplit(as.character(scientificName), " "), length) > 2) %>%
  filter(taxonRank=="subspecies")  %>%
  select(scientificName) %>%
  distinct()
pinned_df %>% 
  filter(taxonRank=="subspecies")  %>%
  select(scientificName) %>%
  distinct()
expert_df %>% 
  filter(taxonRank=="subspecies") %>%
  select(scientificName) %>%
  distinct()


# More cleaning
all_paratax_df <- bet_sorting %>%
  filter(sampleType == "carabid",
         #use data through 2018, since 2019 dataset is incomplete
         collectDate <= "2018-12-31 GMT") %>%
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
  filter(#use data through 2018, since 2019 dataset is incomplete
         collectDate <= "2018-12-31 GMT") %>%
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
  filter(#use data through 2018, since 2019 dataset is incomplete
    collectDate <= "2018-12-31 GMT") %>%
  mutate(scientificName = ifelse(taxonRank=="subspecies",
                                    paste0(genus," ",specificEpithet,sep=""), scientificName),
         col_year = lubridate::year(collectDate)) 
expert_df <- expert_df %>% #index the collection date within each year
  left_join(expert_df %>% 
              distinct(col_year, collectDate) %>%
              group_by(col_year) %>%
              mutate(col_index = 1:n()))


# EDA ---------------------------------------------------------------------

# How many individuals for each expert species ID?
expert_df %>%
  count(scientificName, sort=T) %>%
  data.frame()
# 38

# How many individuals for each parataxonomist species ID?
all_paratax_df %>%
  count(scientificName, sort=T) %>%
  data.frame()
# 44 unique IDs

# How many individuals for each morphospecies ID?
all_paratax_df %>%
  count(morphospeciesID, sort=T) %>% 
  data.frame()
# 24 unique IDs, including NAs

# Look at para_scinames associated with morphospecies IDs.
all_paratax_df %>%
  filter(!is.na(morphospeciesID)) %>%
  count(scientificName, morphospeciesID, sort=T) %>%
  data.frame()

# Do all records not identified to species-level have a morphospecies assigned?
all_paratax_df %>%
  filter(taxonRank != "species" & taxonRank != "subspecies" ) %>%
  count(scientificName, morphospeciesID, sort=T) %>%
  filter(is.na(morphospeciesID))
# 3 records do not identify to species level and also have no morphospecies ID.
# AIS - what to do with this?

# How many unique combinations of para_sciname and morphospecies (when a morph exists)?
all_paratax_df %>%
  count(scimorph_combo, sort=T) %>%
  data.frame()
# 62 unique IDs

# Are there IDs that the parataxonomist ID’d, but that the expert never ID’d?
all_paratax_df %>%
  select(scientificName) %>%
  filter(!scientificName %in% unique(expert_df$scientificName)) %>% #we want rows where the paratax ID is not an existing expert ID
  distinct()
# 14 parataxonomist IDs that the expert never ID'd

# Are there species that the expert taxonomist ID’d, but that the parataxonomist never ID’d?
expert_df %>%
  select(scientificName) %>%
  filter(!scientificName %in% unique(all_paratax_df$scientificName)) %>% #we want rows where the paratax ID is not an existing expert ID
  distinct()
# 8 expert tax IDs that the paratax never ID'd


# Save dataframes locally -------------------------------------------------
saveRDS(all_paratax_df, "occupancy/all_paratax_df.rds")
saveRDS(pinned_df, "occupancy/pinned_df.rds")
saveRDS(expert_df, "occupancy/expert_df.rds")
