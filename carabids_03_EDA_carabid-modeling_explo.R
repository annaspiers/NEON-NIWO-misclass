# This script explores the differences in species ED between the parataxonomist, expert taxonomist, and barcode datasets on NEON ground beetles, then visualizes the data
# Para vs expert comparison written by Max

library(dplyr)
#library(stringr) #word()
library(ggplot2)
#library(forcats) #fct_reorder
#library(arm) #display

# Load carabid data
load(file="data_derived/carabids_NIWO.Rdata")
list2env(carabid_abund, .GlobalEnv)  


# Join dfs
expert <- bet_expertTaxonomistIDProcessed %>%
  mutate(expert_sciname = paste(genus, specificEpithet)) %>%
  dplyr::select(plotID, collectDate, individualID, expert_sciname, sex, nativeStatusCode, exp_sampleCondition = sampleCondition) 
para <- bet_parataxonomistID %>%
  mutate(para_sciname = scientificName) %>%
  dplyr::select(individualID, plotID, trapID,  collectDate, morphospeciesID, taxonRank, para_sciname) 
taxon_df <- left_join(para, expert) 

# Specify species with 100% match between para and expert IDs and with n>20
select_spp <- taxon_df %>%
  mutate(same_species = para_sciname == expert_sciname) %>%
  group_by(para_sciname) %>%
  summarize(pct_correct = mean(same_species, na.rm = TRUE), 
            tot_correct = sum(same_species, na.rm = TRUE), #number of individuals for this species where para = expert where expert exists
            n = n()) %>% #number of individuals with para_sciname
  filter(pct_correct == 1, n > 20) %>%
  pull(para_sciname)

carabid_df <- taxon_df %>% 
  filter(para_sciname %in% select_spp) %>%
  left_join(bet_fielddata %>%
              dplyr::select(plotID, trapID, collectDate, nlcdClass, decimalLatitude, decimalLongitude, elevation, sampleID, sampleCollected, cupStatus, lidStatus, fluidLevel)) %>%
  mutate_at(c("morphospeciesID", "taxonRank", "para_sciname", "expert_sciname", "plotID", "trapID", "nlcdClass"), funs(factor(.))) %>%
  mutate(year = lubridate::year(collectDate), 
         month = lubridate::month(collectDate), 
         day = lubridate::day(collectDate))

# Why are there 11 plots instead of 10?
unique(carabid_df$plotID)
png("output/plots_by_year.png")
carabid_df %>%
  ggplot() +
  geom_point(aes(x= year, y = plotID))
dev.off()
# In 2015-2017, plots 1-10 were sampled. In 2018, plot 4 was replaced by plot 13

# What traps couldn't be sampled and when?
unique(carabid_df$sampleCollected) # All samples were collected - that's almost hard to believe
unique(carabid_df$cupStatus)
nrow(carabid_df %>% filter(cupStatus == "Disturbed")) # Only 39 rows marked distrubed cups
unique(carabid_df$lidStatus)
unique(carabid_df$fluidLevel)
nrow(carabid_df %>% filter(fluidLevel == "None")) # Only 14 rows marked fluid level as "none"

### Modeling EDA
# What variables might be related to carabid abundance?

# Group data by plot and species and summarize by species abundance
model_test_abund <- carabid_df %>%
  group_by(plotID, trapID, para_sciname, nlcdClass) %>% #para_sciname
  summarize(sp_abund = n() ) %>%
  arrange(plotID) %>%
  data.frame()

model_test_comp <- carabid_df %>%
  group_by(plotID, nlcdClass) %>%
  summarize(sp_comp = n_distinct(para_sciname, na.rm = TRUE) ) %>%
  arrange(plotID) %>%
  data.frame()
# AIS later, try grouping by collection date and trapID

# Is habitat type a good predictor for species abundance?
# Constant
fit_abund_0 <- glm(formula = sp_abund ~ 1, 
    family = poisson(link="log"),
      data = model_test_abund)
# With habitat type
fit_abund_1 <- glm(formula = sp_abund ~ nlcdClass, 
    family = poisson(link="log"),
      data = model_test_abund)
fit_abund_2 <- glmer(formula = sp_abund ~ nlcdClass + para_sciname + (1|plotID) , 
    family = poisson(link="log"),
      data = model_test_abund)
display(fit_abund_0)
display(fit_abund_1)
display(fit_abund_2)
anova(fit_abund_0, fit_abund_1, fit_abund_2)

# Is habitat type a good predictor for species composition?
# Constant
fit_comp_0 <- glm(formula = sp_comp ~ 1, 
    family = poisson(link="log"),
      data = model_test_comp)
# With habitat type
fit_comp_1 <- glm(formula = sp_comp ~ nlcdClass, 
    family = poisson(link="log"),
      data = model_test_comp)
display(fit_comp_0)
display(fit_comp_1)


