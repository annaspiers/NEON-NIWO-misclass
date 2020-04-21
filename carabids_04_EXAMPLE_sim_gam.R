#### Example simulations for a GAM ####

library(gratia) 
library(reshape2)
library(dplyr)
library(mgcv)
library(ggplot2) #required in plot_sim_gam function

#the dataset
all7sp_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv") %>%
    mutate(sc_DOY = scale(DOY, center = TRUE, scale = TRUE),
           col_year_fact = as.factor(col_year))

# an example gam
# NEED `control = list(keepData=TRUE)` option in the gam function call so that the data gets saved in the model object
gamm1_cartae_ctrl <- gam(sp_abund ~ s(LAI_1718avg) + s(trap_CHM) + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                        s(plotID, bs="re") + s(plot_trap, bs="re") + 
                        s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                    family=poisson(), method="REML",
                    data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"),
                    control = list(keepData=TRUE)) # this bit keeps the data in the model object



numsims = 50
mod_obj = gamm1_cartae_ctrl
sims = simulate(mod_obj, nsim = numsims, seed = 32)
#columns represent separate simulations, rows are each row in the data.frame,
#values are simulated abundance


# plot to examine - using a function I wrote to avoid lots of lines
source("carabids_04_fxn_plot_sim_gam.R")
plot_sim_gam(sims, mod_obj = mod_obj) 
# simulated values jittered since simulate() produces integer simulated abundances




####### NOTES #######
# gratia package contains the simulate() function for gam objects that is the
# foundation of this script
# https://rdrr.io/cran/gratia/man/simulate.html

# can also add a newdata argument to simulate values for a test data set,
# otherwise it uses the data used for the model

# can access the fitted values, and observed values from the model object, as
# long as the control=list(keepData=TRUE) argument is included



