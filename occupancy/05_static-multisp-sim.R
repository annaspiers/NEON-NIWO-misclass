# This script models the static occupancy of two species by simulation. Following example from Bayesian Population Analysis (Kery & Schaub 2012) section 13.6 and Mackenzie et al. 2004a. This is a multi-state model with possible states: “occupied by both species”, “occupied by species A only”, “occupied by species B only”, "occupied by neither species". We use a multinomial parameterization to describe the probability structure for the occupancy state of a unit

library(jagsUI)

# Select sample sizes (spatial and temporal replication) 
R <- 200 
T <- 3
nspec <- 10


# Generate parameter values -----------------------------------------------

psi <- rbeta(nspec, 4, 1)
psi

p <- rbeta(nspec, 1, 2)
p

# Simulate occupancy states
Z <- matrix(nrow = R, ncol = nspec)
for (i in 1:R) {
    Z[i, ] <- rbinom(nspec, 1, psi)
}



# Simulate data -----------------------------------------------------------

Y <- matrix(nrow = R, ncol = nspec)
for (i in 1:R) {
    Y[i, ] <- rbinom(nspec, T, p * Z[i, ])
}



# Fit model ---------------------------------------------------------------


str(JAGSdata <- list(Y = Y, 
                     T = T, 
                     R = R, 
                     nspec = nspec)) #bundle data
z_init <- (Y > 0) * 1
JAGSinits <- function(){list(Z = z_init) } #rep(JAGSinits, nc) too
JAGSparams <- c("psi", "p") #params monitored
nc <- 3 #MCMC chains
ni <- 2000 #MCMC iterations
nb <- 500 #MCMC burnin
nt <- 2  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/05_static_multisp_sim_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)
names(jags_out)
plot(jags_out) #AIS why is plot of density of occ.fs not a smooth line?

# Next, allow parameters to vary by sampling occassion (bootom of bpa ch 14 p.456)
# how to expand this to three species? how to do this not as a multistate model, but simply as a multispecies model (i.e. columns are species)?
