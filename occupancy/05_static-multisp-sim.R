# This script models the static occupancy of two species by simulation. Following example from Bayesian Population Analysis (Kery & Schaub 2012) section 13.6 and Mackenzie et al. 2004a. This is a multi-state model with possible states: “occupied by both species”, “occupied by species A only”, “occupied by species B only”, "occupied by neither species". We use a multinomial parameterization to describe the probability structure for the occupancy state of a unit

library(jagsUI)

# Select sample sizes (spatial and temporal replication) 
R <- 200 
T <- 3

# Determine process parameters 
# Species given equal probabilities for occupying sites. The following scenarios outlined by Mackenzie et al. 2004
# Scenario 1: {0·4, 0·4, 0·08}, low prob of species occupying a site, species exhibit low association 
# Scenario 2: {0·4, 0·4, 0·24}, low prob of species occupying a site, species exhibit high association 
# Scenario 3: {0·7, 0·7, 0·4}, high prob of species occupying a site, species exhibit low association
# Scenario 2: {0·7, 0·7, 0·6125}, high prob of species occupying a site, species exhibit high association
psi_A <- psi_B <- 0.4
psi_AB <- 0.08
# Assume species to have equal probability of detection, regardless of true state
p_A <- p_B <- r_AB <- r_aB <- r_Ab <- 0.214

# State vector
phi <- array(NA, dim=4)   # initialize
phi[1] <- psi_AB          # Prob species A and B present
phi[2] <- psi_A           # Prob species A present
phi[3] <- psi_B           # Prob species B present
phi[4] <- 1 - psi_AB - psi_A - psi_B     # Prob no species present

# AIS possible states for z are 0,1,2,3

# Ecological process: Sample true (latent) occurrence from a binomial dist
z <- sample(c(0,1,2,3), size=R, prob=phi, replace=TRUE)

# Observation matrix
p <- array(NA, dim=c(4,4) ) # initialize
p[1,1] <- r_AB
p[1,2] <- r_Ab
p[1,3] <- r_aB
p[1,4] <- 1 - r_AB - r_Ab - r_aB
p[2,1] <- 0
p[2,2] <- p_A
p[2,3] <- 0
p[2,4] <- 1 - p_A
p[3,1] <- 0
p[3,2] <- 0
p[3,3] <- p_B
p[3,4] <- 1 - p_B
p[4,1] <- 0
p[4,2] <- 0
p[4,3] <- 0
p[4,4] <- 1

# Observation process: Sample (non)detection 
y <- matrix(NA, nrow = R, ncol = T) # s sites by t sampling occassions
for (s in 1:R) {
    y[s,] <- sample(c(0,1,2,3), size=3, prob=p[z[s],], replace=TRUE)
}

# Run JAGS model
str(JAGSdata <- list(y=y, T=T, R=R)) #bundle data
zst <- apply(y, 1, max) #observed occurrence as starting values for z
JAGSinits <- function(){list(z=zst) } #rep(JAGSinits, nc) too
JAGSparams <- c("psi_A", "psi_B", "psi_AB", "p", "p_A", "p_B", "r_AB", "r_aB", "r_Ab", "n.occ") #params monitored
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
