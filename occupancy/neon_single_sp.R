# This script models the occupancy of a single species at Niwot Ridge from the NEON carabid dataset. Following example from Bayesian Population Analysis (Kery & Schaub 2012) section 13.4.

library(jagsUI)
library(tidyverse)

# Read in the data 
dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE)

# Which species to model?
dat %>%
    filter(para_sciname == "Calathus advena" | para_sciname == "Carabus taedatus" ) %>%
    group_by(para_sciname, plotID, col_year) %>%
    summarize(n=sum(sp_abund)) %>%
    ggplot() +
    geom_line(aes(x = col_year, y = n, colour = para_sciname, group = para_sciname)) +
    geom_point(aes(x = col_year, y = n, colour = para_sciname)) +  
    facet_grid(plotID ~ .) +
    xlab("Collection Year") 
#neither species is present in all plots

# Collect the data into suitable structures
# Start with 2015
CAAD_dat_2015 <- dat %>% 
    filter(para_sciname == "Calathus advena" & col_year == "2015") %>% #choose 1 sp and 1 season
    mutate(occ = ifelse(sp_abund > 0, 1, 0)) %>% # Reduce counts to 0/1
    mutate(scaled_doy = scale(dayofyear, center = TRUE, scale = TRUE)) %>% # Standardise day of year
    mutate(plot_trap = paste(plotID,trapID, sep=""))

y <- CAAD_dat_2015 %>%
    select(plot_trap, collectDate, occ) %>%
    pivot_wider(names_from=collectDate, values_from=occ) %>% 
    arrange(plot_trap) %>%
    select(-plot_trap) #remove plot_trap IDs so they are not part of the analysis

# Grab covariates
nlcd <- CAAD_dat_2015  %>%
    select(plot_trap, nlcdClass) %>%
    distinct() %>% 
    arrange(plot_trap) %>%
    mutate(nlcdCat = ifelse(nlcdClass == "evergreenForest",1,0)) %>%
    select(-c(plot_trap,nlcdClass)) #remove plot_trap IDs so they are not part of the analysis
#AIS nlcdClass would not work if nlcdClass was a factor, had to convert to 0/1. Is this a rule of thumbs with JAGS (no factors, only numeric)?
sc_doy <- CAAD_dat_2015 %>% 
    select(plot_trap, collectDate, scaled_doy) %>%
    pivot_wider(names_from=collectDate, values_from=scaled_doy) %>% 
    arrange(plot_trap) %>%
    select(-plot_trap) #remove plot_trap IDs so they are not part of the analysis


# Run JAGS model
str(JAGSdata <- list(y=y, T=ncol(y), R=nrow(y), nlcd=nlcd, sc_doy=sc_doy)) #bundle data
zst <- apply(y, 1, max) # Good starting values crucial
JAGSinits <- function(){
    list(z = zst, 
         alpha.psi=runif(1, -3, 3), #AIS why choose this dist?
         alpha.p = runif(1, -3, 3))}
JAGSparams <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p")
nc <- 3 #MCMC chains        #AIS, jags required the same number of chains as init values. only one init variable
ni <- 3000 #MCMC iterations
nb <- 2000 #MCMC burnin #AIS (ni=3000, nb=20000) got an error message saying that number of iterations must be larger than burnin - why?
nt <- 10  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/neon_single_sp_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)
names(jags_out)
plot(jags_out) #AIS why is plot of density of occ.fs not a smooth line?

