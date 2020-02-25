# This script compiles carabid abundance and richness data with variables
# describing layers of data aggregation as well as environmental variables. Each
# row is a species at a trap on a collection date

library(dplyr)
library(stringr) #word()
library(ggplot2)
library(forcats) #fct_reorder

### Load data ####
load("data_derived/carabids_NIWO.Rdata")
#load("data_derived/soil_wc_NIWO.Rdata")
load("data_derived/woody_veg_NIWO.Rdata")
load("data_derived/litter_woodfall_NIWO.Rdata")
load("data_derived/ir_bio_temp_NIWO.Rdata")
#load("data_derived/soil_temp_NIWO.Rdata")
load("data_derived/precip_NIWO.Rdata")
load("data_derived/rad_net_NIWO.Rdata")
load("data_derived/rad_short_dir_diff_NIWO.Rdata")

### Combine into one df ###

model_df <- 
    

bet_parataxonomistID %>%
    as_tibble %>% 
    select(individualID, siteID, plotID, trapID, setDate, collectDate, taxonID, taxonRank, morphospeciesID, scientificName) %>%
    mutate(para_sciname = scientificName)

    

    
utate(para_sciname = scientificName) %>%
            select(individualID, morphospeciesID, taxonRank, para_sciname) %>%
            as_tibble


# site-level variables to include: 
# plot-level variables to include: nlcdcover, soil type
# trap-level variables to include: 

# Plot spatial arrangement of plots, labeled with habitat
# add lat and long poitns
bet_fielddata %>% 
  select(plotID, nlcdClass, decimalLatitude, decimalLongitude) %>%
  ggplot() +
  geom_point(aes(x = decimalLongitude, y = decimalLatitude, colour = nlcdClass))

# cool variables to look at: bet_fielddata$nativestatuscode
bet_parataxonomistID %>% 
  filter(scientificName %in% select_spp) %>%
  select(scientificName, nativeStatusCode) %>%
  summarize(mean(nativeStatusCode == "N") )
# All species are classified as native

taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    select(individualID, plotID , trapID , collectDate)) %>% 
  filter(para_sciname %in% select_spp) %>%
  mutate(year = lubridate::year(collectDate), 
         month = lubridate::month(collectDate), 
         day = lubridate::day(collectDate)) %>% 
  group_by(para_sciname, year, month, day) %>%
  summarize(n=n()) 



#get soil type order
carabid_spat <- def.extr.geo.os(data = carabid_abund$bet_fielddata, 'namedLocation', locOnly=T) %>%
  dplyr::select(api.decimalLatitude, api.decimalLongitude, Value.for.Plot.ID, api.soilTypeOrder) %>%
  cbind('data_type' = rep('carabid'))
#########


### Save df ###







