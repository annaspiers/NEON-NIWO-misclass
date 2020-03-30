# This script models the dynamic occupancy of a single species over a few seasons at Niwot Ridge from the NEON carabid dataset. Following example from Bayesian Population Analysis (Kery & Schaub 2012) section 13.5.

library(jagsUI)
library(tidyverse)
library(reshape2)

# Read in the data 
dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE)

# Collect the data into suitable structures
CAAD_dat <- dat %>% 
    filter(para_sciname == "Calathus advena" ) %>% #choose 1 sp 
    mutate(occ = ifelse(sp_abund > 0, 1, 0),
           scaled_doy = scale(dayofyear, center = TRUE, scale = TRUE),
           plot_trap = paste(plotID,trapID, sep=""))

secondary_period_df <- CAAD_dat %>%
    distinct(col_year, collectDate) %>%
    filter(!is.na(col_year)) %>%
    group_by(col_year) %>%
    mutate(secondary_index = 1:n()) %>%
    ungroup

# dimensions of y: 
# 1. plot_trap
# 2. number sampling occasions per year
# 3. number of years
y <- CAAD_dat %>%
    left_join(secondary_period_df) %>%
    distinct(plot_trap, secondary_index, col_year, occ) %>%
    filter(!is.na(col_year)) %>%
    reshape2::acast(plot_trap ~ secondary_index ~ col_year, 
                    value.var = "occ")
str(y) # dim is [1:44, 1:7, 1:4] 44 sites (plot_trap), up to 7 surveys, 4 years

# Look at the number of sites with detections for each year     
tmp <- apply(y, c(1,3), max, na.rm = TRUE) 
tmp[tmp == "-Inf"] <- NA 
apply(tmp, 2, sum, na.rm = TRUE) #[1] 0031017176


# # Grab covariates
# nlcd <- CAAD_dat  %>%
#     select(plot_trap, nlcdClass) %>%
#     distinct() %>%
#     arrange(plot_trap) %>%
#     mutate(nlcdCat = ifelse(nlcdClass == "evergreenForest",1,0)) %>%
#     select(-c(plot_trap,nlcdClass)) #remove plot_trap IDs so they are not part of the analysis
# #AIS nlcdClass would not work if nlcdClass was a factor, had to convert to 0/1. Is this a rule of thumbs with JAGS (no factors, only numeric)?
# sc_doy <- CAAD_dat %>%
#     select(plot_trap, collectDate, scaled_doy) %>%
#     pivot_wider(names_from=collectDate, values_from=scaled_doy) %>%
#     arrange(plot_trap) %>%
#     select(-plot_trap) #remove plot_trap IDs so they are not part of the analysis


# Run JAGS model
str(JAGSdata <- list(y = y,       #dim(y) is 44 7 4
                     nsite = dim(y)[1], 
                     nrep = dim(y)[2], 
                     nyear = dim (y)[3]))
# Initial values 
zst <- apply(y, c(1,3), max, na.rm = TRUE)
zst[zst == "-Inf"] <- NA 
# Observed occurrence as inits for z 
JAGSinits <- function() { list(z = zst) } #dim(z) is 44 4
JAGSparams <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")

# MCMC settings
nc <- 3     #MCMC chains
ni <- 5000  #MCMC iterations
nb <- 1000  #MCMC burnin
nt <- 4     #MCMC thin

# JAGS model
JAGSout <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/04_neon_dynamic_single-sp_JAGS.txt",
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(JAGSout, dig=2)
names(JAGSout)
plot(JAGSout) #AIS why is plot of density of occ.fs not a smooth line?
