# This script is models the static occupancy of NEON carabid speciesat Niwot Ridge. Data include 7 most abundant and accurately-identified species

library(jagsUI)
library(tidyverse)
library(reshape2)


# Read in and clean data  -------------------------------------------------------

dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE) %>%
    mutate(occ = ifelse(sp_abund > 0, 1, 0),
           plot_trap = paste(plotID,trapID, sep=""))


# Create detection array --------------------------------------------------

y <- dat %>%
    filter(col_year == 2018) %>%
    reshape2::acast(plot_trap ~ para_sciname ~ col_index, 
                    value.var = "occ")
str(y) # dim is [1:40, 1:7, 1:7] 40 traps, 7 species, 7 surveys
image(t(y[,,3]))
#table(y[,,3])


# Fit model ---------------------------------------------------------------

str(JAGSdata <- list(y = y, 
                     ntrap = dim(y)[1], 
                     nspec = dim(y)[2], 
                     nsurv = dim(y)[3])) #bundle data
zst <- apply(y, c(1,2), max, na.rm = TRUE) 
JAGSinits <- function(){list(z = zst) } 
JAGSparams <- c("psi", "p") #params monitored
nc <- 4 #MCMC chains #AIS more chains decreases n.eff
ni <- 5000 #MCMC iterations
nb <- 1200 #MCMC burnin - increasing from 500 to 1200 got rid of the error "Slicer stuck at value with infinite density"
nt <- 2  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/07_neon_static_multisp_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)
names(jags_out)
plot(jags_out) 

