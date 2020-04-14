# This script is models the dynamic occupancy of NEON carabid speciesat Niwot Ridge 2015-2018. Data include 7 most abundant and accurately-identified species

library(jagsUI)
library(dplyr)
library(reshape2)


# Read in and clean data  -------------------------------------------------------

dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE) %>%
    mutate(occ = ifelse(sp_abund > 0, 1, 0),
           plot_trap = paste(plotID,trapID, sep=""))


# Create detection array --------------------------------------------------

# dimensions of y: [plot_trap, species, surveys within year, year]
y <- dat %>%
    reshape2::acast(plot_trap ~ para_sciname ~ col_index ~ col_year, 
                    value.var = "occ")
str(y) # dim is [1:40, 1:7, 1:7] 40 traps, 7 species, 7 surveys
image(t(y[,,1,1])) #AIS why is top row blank?
image(t(y[,,2,1]))
image(t(y[,,3,1]))
image(t(y[,,4,1]))
image(t(y[,,5,1]))
#table(y[,,1,1])
#table(y[,,5,1])


# Fit model ---------------------------------------------------------------

str(JAGSdata <- list(y = y, 
                     ntrap = dim(y)[1], 
                     nspec = dim(y)[2], 
                     nsurv = dim(y)[3],
                     nyear = dim(y)[4])) 
zst <- apply(y, c(1,2,4), max, na.rm = TRUE) # closure across surveys. varies by trap, species, and year.
zst[zst == "-Inf"] <- NA #this happens where a plot isn't sampled across years 
# AIS how does JAGS handle inits (i.e. compared to handling data)?
JAGSinits <- function(){list(z = zst) } 
JAGSparams <- c("psi", "p", "phi", "gamma", "n.occ", "growth", "turnover") 
nc <- 4     #MCMC chains 
ni <- 10000  #MCMC iterations
nb <- 1800  #MCMC burnin
nt <- 2     #MCMC thin

# JAGS model
JAGSout <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/08_neon_dynamic_multisp_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(JAGSout, dig=2)
names(JAGSout) 
plot(JAGSout) # horizontal lines through the traceplots show convergence visually


# Notes on comparing data to model output  --------------------------------

dimnames(y)

#compare to "data_derived/species_abund_by_year_by_trap.png"

JAGSout$mean$psi
# Calathus advena occupancy consistent through years
# Species 2 not present in survey data in 2015
# Species 7 not present in survey data in 2015-2016

JAGSout$mean$p

JAGSout$mean$phi

JAGSout$mean$gamma
# High colonization for gamma[4] - doesn't make sense with data

JAGSout$mean$n.occ
# Species 2 and 7 have low n.occ in first year(s) - consistent with data

JAGSout$mean$growth
# growth[,2] for species 2 and 7 are unrealistically large due to 0 detection in first year(s)
# Species 6 vlaues seem too high

JAGSout$mean$turnover


