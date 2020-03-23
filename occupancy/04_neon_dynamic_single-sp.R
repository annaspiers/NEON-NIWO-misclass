# This script models the dynamic occupancy of a single species over a few seasons at Niwot Ridge from the NEON carabid dataset. Following example from Bayesian Population Analysis (Kery & Schaub 2012) section 13.5.

library(jagsUI)
library(tidyverse)
library(unmarked)

# Read in the data 
dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE)

# Collect the data into suitable structures
CAAD_dat <- dat %>% 
    filter(para_sciname == "Calathus advena" ) %>% #choose 1 sp 
    complete(collectDate, plotID, trapID) %>% #AIS added in this script so that NA's show up for plots not sampled in a given year
    mutate(occ = ifelse(sp_abund > 0, 1, 0)) %>% # Reduce counts to 0/1
    mutate(scaled_doy = scale(dayofyear, center = TRUE, scale = TRUE)) %>% # Standardise day of year
    mutate(plot_trap = paste(plotID,trapID, sep=""))

# Generate detection data for each season
y <- array(NA, dim = c(40, 7, 4) ) # 40 sites per year, max of 7 reps per year, 4 years
# First try 7 surveys per season
for (k in 1:4) {
    year_subset <- CAAD_dat %>%
        filter(col_year == c(2015, 2016, 2017, 2018)[k]) %>%
        select(plot_trap, collectDate, occ) %>%
        pivot_wider(names_from=collectDate, values_from=occ) %>% 
        arrange(plot_trap) %>% 
        select(-plot_trap) 
    y[,,k] <- as.matrix(year_subset)
}

# Next, try the accurate number of surveys per season
y <- list(NA, dim = 4 ) 
for (k in 1:4) {
    sel.rows <- CAAD_dat %>%
        expand(nesting(plotID, trapID)) %>%
        filter(col_year == c(2015, 2016, 2017, 2018)[k]) %>%
        select(plot_trap, collectDate, occ) %>%
        pivot_wider(names_from=collectDate, values_from=occ) %>% 
        arrange(plot_trap) %>% 
        select(-plot_trap) 
    y[[k]] <- as.matrix(sel.rows)
}
names(y) <- c("2015", "2016", "2017", "2018")
str(y)

df <- tibble(
  group = c(1:2, 1),
  item_id = c(1:2, 2),
  item_name = c("a", "b", "b"),
  value1 = 1:3,
  value2 = 4:6
)
df %>% complete(group, nesting(item_id, item_name))




    for (j in 1:length(y[1,,1]) ) {
        for (i in 1:length(y[,1,1]) ) {
            
        }
    }
# AIS how to deal with irregular plot surveying


for (k in 2015:2018)  {
    y[,,k] <- array(NA, dim = c(nrow(CAAD_dat %>%
                                       filter(col_year == k) %>%
                                       distinct(plot_trap) ),
                            nrow(CAAD_dat %>%
                                       filter(col_year == k) %>%
                                       distinct(collectDate) ) ) )
}
    
y_2015 <- array(NA, dim = c(nrow(CAAD_dat %>%
                                       filter(col_year == 2015) %>%
                                       distinct(plot_trap) ),
                            nrow(CAAD_dat %>%
                                       filter(col_year == 2015) %>%
                                       distinct(collectDate) ) ) )


# y <- CAAD_dat %>%
#     select(plot_trap, collectDate, occ) %>%
#     pivot_wider(names_from=collectDate, values_from=occ) %>% 
#     arrange(plot_trap) %>%
#     select(-plot_trap) #remove plot_trap IDs so they are not part of the analysis
# 
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
# 
# 
# # Run JAGS model
# str(JAGSdata <- list(y=y, 
#                      T=ncol(y), 
#                      R=nrow(y), 
#                      nlcd=unlist(nlcd), 
#                      sc_doy=sc_doy)) #bundle data
# zst <- apply(y, 1, max) # Good starting values crucial
# #JAGSinits <- function(){list(z = zst)}
# JAGSinits <- function(){
#     list(z = zst, 
#          alpha.psi=runif(1, -3, 3), #AIS why choose this dist?
#          alpha.p = runif(1, -3, 3))}
# JAGSparams <- c("psi", "p", "occ.fs") 
# nc <- 3 #MCMC chains        #AIS, jags required the same number of chains as init values. only one init variable
# ni <- 1500 #MCMC iterations
# nb <- 200 #MCMC burnin #AIS (ni=3000, nb=20000) got an error message saying that number of iterations must be larger than burnin - why?
# nt <- 5  #MCMC thin
# 
# # JAGS model
# jags_out <- jags(data = JAGSdata,
#                  inits = JAGSinits,
#                  parameters.to.save = JAGSparams,
#                  model.file = "occupancy/neon_single-season_single-sp_JAGS.txt", 
#                  n.chains = nc,
#                  n.iter = ni,
#                  n.burnin = nb,
#                  n.thin = nt)
# 
# print(jags_out, dig=2)
# names(jags_out)
# plot(jags_out) #AIS why is plot of density of occ.fs not a smooth line?
# 

# Following Kery & Chandler 2016 section 4 example
data(crossbill)
DATE <- as.matrix(crossbill[,32:58])
y.cross <- as.matrix(crossbill[,5:31])
y.cross[is.na(DATE) != is.na(y.cross)] <- NA
sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE
years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross,
                         siteCovs=crossbill[,2:3], yearlySiteCovs=list(year=years),
                         obsCovs=list(date=DATE),
                         numPrimary=9)
 