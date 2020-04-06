# This script is models the static occupancy of NEON carabid speciesat Niwot Ridge. Data include 7 most abundant and accurately-identified species

library(jagsUI)
library(tidyverse)
library(reshape2)


# Read in and clean data  -------------------------------------------------------

raw_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE)

clean_dat <- raw_dat %>% 
    mutate(occ = ifelse(sp_abund > 0, 1, 0),
           plot_trap = paste(plotID,trapID, sep=""))

# Choose one year to model
dat_2018 <- clean_dat %>%
    filter(col_year == 2018)


# Create detection array --------------------------------------------------

y <- dat_2018 %>%
    #distinct(plot_trap, col_index, occ) %>%
    reshape2::acast(plot_trap ~ para_sciname ~ col_index, 
                    value.var = "occ")
str(y) # dim is [1:40, 1:7, 1:7] 40 traps, 7 species, 7 surveys
#image(t(y[,,3]))
#table(y[,,3])


# Fit model ---------------------------------------------------------------

str(JAGSdata <- list(y = y, 
                     ntrap = dim(y)[1], 
                     nspec = dim(y)[2], 
                     nsurv = dim (y)[3])) #bundle data
zst <- apply(y, c(1,2), max, na.rm = TRUE) #AIS what does zst mean, intuitively? z start?
zst[zst == "-Inf"] <- NA #AIS why would a value of -Inf be generated?
JAGSinits <- function(){list(z = zst) } 
JAGSparams <- c("psi", "p") #params monitored
nc <- 3 #MCMC chains
ni <- 2500 #MCMC iterations
nb <- 500 #MCMC burnin
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
# AIS 

# return in the morning - make a 08 script that is dynamic, multispecies. or try incorporating all samples, even rare ones. Also, try incorporating covariates!!
# AIS this model assumes closure within a season. Allow z to vary by survey to incorporate seasonality
# Next, allow parameters to vary by sampling occassion (bootom of bpa ch 14 p.456)
