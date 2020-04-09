# This script explores the differences in species ED between the parataxonomist, expert taxonomist, and barcode datasets on NEON ground beetles, then visualizes the data
# Para vs expert comparison written by Max

library(dplyr)
library(stringr) #word()
library(ggplot2)
library(forcats) #fct_reorder
library(arm) #display

# Load carabid data
load(file="data_derived/carabids_NIWO.Rdata")
list2env(carabid_abund, .GlobalEnv)  
list2env(carabid_barcode, .GlobalEnv)  


### First compare barcode and expert taxonomist

barcode <- bet_BOLDtaxonomy %>%
            mutate(individualID = sampleID, bc_sciname = NA) %>%
            as_tibble

expert <- bet_expertTaxonomistIDProcessed %>%
            mutate(expert_sciname = paste(genus, specificEpithet)) %>%
            as_tibble

# Some barcode species are subspecies. Remove subspecies, leaving genus/species
for (i in 1:nrow(barcode)) {
    if (sapply(strsplit(barcode$species[i]," "),length) == 2) {
        barcode$bc_sciname[i] <- barcode$species[i]
    } else {
        barcode$bc_sciname[i] <- paste(word(barcode$species[i], 1:2),collapse=" ")
    }
}
# AIS couldn't get this to work in simple ifelse format: barcode %>% mutate(bc_sciname = ifelse( sapply(strsplit(species," "),length) == 2, species, paste(word(species, 1:2),collapse=" ")))

# Join barcode and expert data
bc_exp_taxon_df <- expert %>%
            dplyr::select(individualID, expert_sciname) %>%
            left_join(distinct(barcode, individualID, bc_sciname))

# What percentage of barcode IDs match expert IDs 
bc_exp_taxon_df %>%
  filter(!is.na(bc_sciname)) %>% #samples that expert IDed and where para identified to species
  mutate(same_species = expert_sciname == bc_sciname) %>%
  group_by(expert_sciname) %>%
  summarize(pct_correct = mean(same_species, na.rm = TRUE), 
            tot_correct = sum(same_species, na.rm = TRUE), 
            n = n()) %>%
  arrange(-n) %>%
  data.frame
#100%!


### Now compare barcode/expert to parataxonomist

para <- bet_parataxonomistID %>%
            mutate(para_sciname = scientificName) %>%
            dplyr::select(individualID, morphospeciesID, taxonRank, para_sciname) %>%
            as_tibble

# Join para and expert/barcode data, and evaluate discrepancies in species ids 
taxon_df <- left_join(para, bc_exp_taxon_df)

# para_sciname and expert_sciname don't always match
taxon_df %>%
  filter(taxonRank == "species") %>%
  summarize(mean(para_sciname == expert_sciname, na.rm = TRUE)) 
#0.97 match between specimen identified by expert taxonomist with parataxonomist ID

# Are any species perfectly identified?
taxon_df %>%
  filter(taxonRank == "species", !is.na(expert_sciname)) %>% #samples that expert IDed and where para identified to species 
  mutate(same_species = para_sciname == expert_sciname) %>%
  group_by(para_sciname) %>%
  summarize(pct_correct = mean(same_species, na.rm = TRUE), 
            tot_correct = sum(same_species, na.rm = TRUE), 
            n = n()) %>%
  arrange(-n) %>%
  data.frame
# In all but one case (Cymindis cribricollis), the para and expert IDs were either 0% or 100% matching

# look at some of the discrepancies
taxon_df %>%
  filter(taxonRank == "species", 
         para_sciname != expert_sciname) %>%
  select(individualID, para_sciname, expert_sciname) %>%
  arrange(para_sciname) %>%
  data.frame
# Three mismatches checked with taxize pckg had different taxon serial numbers

# Plot discrepancies for samples identified to the species level
taxon_df %>%
  filter(taxonRank == "species", !is.na(expert_sciname)) %>%
  count(para_sciname, expert_sciname) %>%
  arrange(n) %>%
  mutate(discrepancy = para_sciname != expert_sciname) %>%
  ggplot(aes(x = para_sciname, 
             y = expert_sciname, 
             color = discrepancy)) + 
  geom_point(aes(size = n)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Species ID from parataxonomist") + 
  ylab("Species ID from expert taxonomist") + 
  scale_color_manual(values = c("black", "red"))


# Plot morphospecies IDs with expert IDs: 
taxon_df %>%
  filter(!is.na(morphospeciesID), !is.na(expert_sciname)) %>% #samples where morphospeciesID and expert ID exists,
  count(morphospeciesID, expert_sciname) %>%
  ggplot(aes(x = fct_reorder(morphospeciesID, n), 
             y = fct_reorder(expert_sciname, n))) + 
  geom_point(aes(size = n)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Morphospecies from parataxonomist") + 
  ylab("Species ID from expert taxonomist")
# all but one morphospecies map to one expert species ID.
# the D13.2016.MorphBU morphospecies is a mixture of two species

# How many morphospecies occur with each parataxonomist ID?
# note: if the parataxonomist IDed to species, then no morphospecies was assigned
taxon_df %>%
  group_by(para_sciname) %>%
  summarize(morpho_n = n_distinct(morphospeciesID, na.rm = TRUE), 
            n = n()) %>%
  arrange(-n) %>%
  data.frame
# 3 para ID's had more than one morphospecies
# Harpalus sp. (para ID) had 6 morphospecies and Amara sp. (para ID) had 9
# Doesn't seem like using the morphospeciesID is helpful

# Plot discrepancies for all samples, excluding para ID == expert ID
taxon_df %>%
  count(para_sciname, expert_sciname) %>%
  arrange(n) %>%
  mutate(discrepancy = para_sciname != expert_sciname) %>%
  filter(discrepancy == TRUE) %>%
  ggplot(aes(x = para_sciname, 
             y = expert_sciname)) + 
  geom_point(aes(size = n)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Species ID from parataxonomist") + 
  ylab("Species ID from expert taxonomist") + 
  scale_color_manual(values = c("red"))
  

# Visualize those species with 100% match between para and expert IDs and with n>20
select_spp <- taxon_df %>%
  mutate(same_species = para_sciname == expert_sciname) %>%
  group_by(para_sciname) %>%
  summarize(pct_correct = mean(same_species, na.rm = TRUE), 
            tot_correct = sum(same_species, na.rm = TRUE), #number of individuals for this species where para = expert where expert exists
            n = n()) %>% #number of individuals with para_sciname
  filter(pct_correct == 1, n > 20) %>%
  pull(para_sciname)


# Plot abundance of selected Niwot species over time
taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    dplyr::select(individualID, plotID , trapID , collectDate)) %>% 
  filter(para_sciname %in% select_spp) %>%
  ggplot() +
  geom_bar(aes(x = collectDate, fill = para_sciname)) +
  facet_wrap(. ~ para_sciname)

# Plot taxa by year (thanks SN)
taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    select(individualID, plotID , trapID , collectDate)) %>% 
  filter(para_sciname %in% select_spp) %>%
  mutate(year = lubridate::year(collectDate), 
         month = lubridate::month(collectDate), 
         day = lubridate::day(collectDate)) %>% 
  group_by(para_sciname, year) %>%
  summarize(n=n())  %>%
  ggplot() +
  geom_line(aes(x = year, y = n, colour = para_sciname, group = para_sciname)) +
  geom_point(aes(x = year, y = n, colour = para_sciname)) +  
  xlab("Collection Year") 
# Pterostichus restrictus came out of nowhere!

# Plot taxa by month through the years
taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    select(individualID, plotID , trapID , collectDate)) %>% 
  filter(para_sciname %in% select_spp) %>%
  mutate(year = lubridate::year(collectDate), 
         month = lubridate::month(collectDate), 
         day = lubridate::day(collectDate)) %>% 
  group_by(para_sciname, year, month) %>%
  summarize(n=n())  %>%
  ggplot() +
  geom_line(aes(x = month, y = n, colour = para_sciname, group = para_sciname)) +
  geom_point(aes(x = month, y = n, colour = para_sciname))  +  
  xlab("Collection Month") +
  facet_grid(year ~ .)
# July has the most abundance, but may have the most collection dates, too

# Plot taxa through time by collection day
taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    select(individualID, plotID , trapID , collectDate)) %>% 
  filter(para_sciname %in% select_spp) %>%
  mutate(year = lubridate::year(collectDate), 
         month = lubridate::month(collectDate), 
         day = lubridate::day(collectDate)) %>% 
  group_by(para_sciname, year, month, day) %>%
  summarize(n=n())  %>%
  ggplot() +
  #geom_line(aes(x = day, y = n, colour = para_sciname, group = para_sciname)) +
  geom_point(aes(x = day, y = n, colour = para_sciname)) +  
  xlab("Collection Day") +
  facet_grid(year ~ month)

# Plot spatial arrangement of plots, labeled with habitat
# add lat and long poitns
bet_fielddata %>% 
  select(plotID, nlcdClass, decimalLatitude, decimalLongitude) %>%
  ggplot() +
  geom_point(aes(x = decimalLongitude, y = decimalLatitude, colour = nlcdClass))

# Visualize what species are in which nlcdClass
taxon_df %>% 
  full_join(bet_parataxonomistID %>% 
    dplyr::select(individualID, plotID , trapID , collectDate)) %>% 
  left_join(bet_fielddata %>%
             dplyr::select(plotID , nlcdClass)) %>% 
  filter(para_sciname %in% select_spp) %>%
  ggplot() +
  geom_bar(aes(x = collectDate, fill = para_sciname)) +
  facet_wrap(para_sciname ~ nlcdClass)

# cool variables to look at: bet_fielddata$nativestatuscode
bet_parataxonomistID %>% 
  filter(scientificName %in% select_spp) %>%
  select(scientificName, nativeStatusCode) %>%
  summarize(mean(nativeStatusCode == "N") )
# All species are classified as native


