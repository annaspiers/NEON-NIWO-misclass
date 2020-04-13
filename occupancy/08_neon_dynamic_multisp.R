# This script is models the dynamic occupancy of NEON carabid speciesat Niwot Ridge 2015-2018. Data include 7 most abundant and accurately-identified species

library(jagsUI)
library(dplyr)
library(reshape2)


# Read in and clean data  -------------------------------------------------------

dat <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE) %>%
    mutate(occ = ifelse(sp_abund > 0, 1, 0),
           plot_trap = paste(plotID,trapID, sep=""))


# Create detection array --------------------------------------------------

# dimensions of y: 
# 1. plot_trap
# 2. species
# 3. number sampling occasions per year
# 4. number of years
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
# with autojags, I could see that growth[2,2], gamma[7,1],gamma[1,1], gamma[1,3], gamma[2,3], gamma had trouble converging
# AIS species 2 and 7 were not present at beginning of survey 92015 and 2015-2016, respectively) - see output photo

print(JAGSout, dig=2)
names(JAGSout) #JAGout$mean, see growth[2,2] and growth[7,3] - too big!
#JAGSout$mean - most p values are larger than corresponding psi. Is p conditional on phi?
plot(JAGSout) #nonconverged: psi[2,1], psi[7,1], psi[7,2]. many of the phi's look wonky. most of the gamma's look wonky too. growth[2,2], growth[7,2], growth[7,3] DNE. bimodal turnover[7,1]
#AIS what are the horizontal lines through the trace plots?
View(JAGSout)
whiskerplot(JAGSout, c("psi"), zeroline=T)
pp.check(JAGSout, observed, simulated, xlab=NULL, ylab=NULL, main=NULL)


# IN THE MORNING - how to visualize these results nicely?

# AIS Next, allow parameters to vary by sampling occassion (bootom of bpa ch 14 p.456)
