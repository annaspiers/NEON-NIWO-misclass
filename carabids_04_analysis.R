# We try different models of predicting carabid abundance and composition

library(arm)
library(dplyr)

dat <- read.csv("data_derived/model_df_by_species_in_sample.csv")

# Pick one species for now
dat %>% group_by(para_sciname) %>%
    summarize(total = sum(sp_abund)) %>%
    arrange(-total)
# Calathus advena is most abundant
dat_caladv <- dat %>% 
    filter(para_sciname == "Calathus advena") %>%
    arrange(collectDate)
    


# Try modeling Calathus advena

# Harris' baseline average model
# Investivating tthe scales of variance in the datta (aka variance components model). In Dietze book, he uses the average model as an opportunity of exploring the scales of variance in the data. 
# Between trap variation and between plot variation.
# then you have multiple samples per trap through time
# Where is most of the variation? Is it between plots or between traps?
# This average model with a random effects structure is farily sophisticated at different scales.

# See radon example

# Harris' baseline naive model
