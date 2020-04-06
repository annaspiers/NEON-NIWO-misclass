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
#load("data_derived/precip_NIWO.Rdata")
#load("data_derived/rad_net_NIWO.Rdata")
#load("data_derived/rad_short_dir_diff_NIWO.Rdata")
load("data_derived/trap_LAI_sp.Rdata")

### Combine into one df ###
select_spp <- c("Amara alpina", "Amara quenseli", "Calathus advena", "Carabus taedatus", "Cymindis unicolor", "Harpalus nigritarsis", "Pterostichus restrictus")


# I Create df where each row is an individual beetle ------------------------

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



# II Create df where each row is a species in a trap from a sing --------

#AIS edit this to use complete()

# How many rows should this df have?
length(unique(model_df_by_ind$collectDate)) #24 unique collectDates
# 960 rows per species = (24 collectDates)*(4 traps)*(10 plots)

# initialize 960 row df 
#initialize df with plot 13 in 2018
dates2018 <- c("2018-06-05","2018-06-19","2018-07-03","2018-07-17","2018-07-31","2018-08-14","2018-08-28")
all_samples_df <- data.frame(plotID = rep("NIWO_013", 28), 
                             trapID = c(rep("N",7),rep("E",7),rep("S",7),rep("W",7)), 
                             collectDate = rep(as.factor(dates2018), 4) )
# then add on each unique collectDate by trap by plot
for ( i in 1:(length(unique(model_df_by_ind$plotID))-1) ) {
  for ( j in 1:length(unique(model_df_by_ind$trapID)) ) {
    plot_i_trap_j_df <- data.frame(plotID = rep(unique(model_df_by_ind$plotID)[i], 24), 
                                   trapID = rep(unique(model_df_by_ind$trapID)[j], 24), 
                                   collectDate = rep(unique(as.factor(model_df_by_ind$collectDate)), 1) )
    all_samples_df <- rbind(all_samples_df, plot_i_trap_j_df)
  }
}
# remove rows for  plot 4 in 2018
all_samples_df <- all_samples_df %>% 
  filter( !(plotID == "NIWO_004" & collectDate %in% dates2018) )
nrow(all_samples_df)# Now we have 960 rows

# Create sampleID column
all_samples_df <- all_samples_df %>%
  mutate(sampleID = paste(plotID, trapID, gsub("-","",collectDate), sep="."))
length(unique(all_samples_df$sampleID)) #check that there are 960 unique sampleIDs

# Add on descriptor columns (not species and abundance)
descriptors_df <- model_df_by_ind %>%
    group_by(siteID, plotID, decimalLatitude, decimalLongitude, elevation, nlcdClass, soilOrder) %>%
  summarize()
all_samp_desc_df <- all_samples_df %>%
  mutate(col_year = lubridate::year(collectDate),
         col_month = lubridate::month(collectDate),
         col_day = lubridate::day(collectDate),
         dayofyear = as.numeric(strftime(collectDate, format = "%j"))) %>%
  mutate_at(c("plotID", "trapID"), funs(as.character(.))) %>%
  left_join(descriptors_df) %>%
  arrange(collectDate)
nrow(all_samp_desc_df %>% filter(is.na(siteID)==TRUE))

# Multiply out df for number of species
temp_df <- all_samp_desc_df %>% 
  mutate(para_sciname = NA) 
for(i in 1:length(select_spp)) {
    df_i <- all_samp_desc_df %>% 
      mutate(para_sciname = select_spp[i]) 
    temp_df <- rbind(temp_df, df_i)
    rm(df_i)
}
# remove initial df chunk with para_sciname values as NA
temp_df <- temp_df %>% filter(is.na(para_sciname)==FALSE) 
nrow(temp_df) #df is 960*(7 species)=6720 rows

# Create df summarizing all samples that have >0 abundance
all_samps_w_abund_df <- model_df_by_ind %>%
    group_by(sampleID, para_sciname) %>%
  summarize(sp_abund=n())

# Merge temp_df and all_samps_w_abund_df so each possible species sample is joined with its abundance
model_df_by_sample <- temp_df %>%
  left_join(all_samps_w_abund_df)
nrow(model_df_by_sample) #6270 rows

# Convert sp_abund NAs to 0s
model_df_by_sample$sp_abund[is.na(model_df_by_sample$sp_abund)] <- 0
unique(model_df_by_sample$sp_abund) #no NAs

# Index collection date within each year
secondary_period_df <- model_df_by_sample %>%
    distinct(col_year, collectDate) %>%
    filter(!is.na(col_year)) %>%
    group_by(col_year) %>%
    mutate(col_index = 1:n()) 
model_df_by_sample <- model_df_by_sample %>%
  left_join(secondary_period_df)
  
# Add in more predictor variables -----------------------------------------

#LAI
trap_LAI_sp_edited <- trap_LAI_sp@data %>%
  separate(locationNa, into = c('plotID','junk','junk2','trapID'), sep="\\.", remove=FALSE) %>%
  dplyr::select(-c(locationNa, junk, junk2), trap.Easting=x, trap.Northing=y)

model_df_by_sample <- model_df_by_sample %>%
  left_join(trap_LAI_sp_edited)

### Save df ###
  
write.csv(model_df_by_ind, file="data_derived/model_df_by_individual_beetle.csv") #df by individuals
write.csv(model_df_by_sample, file="data_derived/model_df_by_species_in_sample.csv") #df by para_sciname in sample


