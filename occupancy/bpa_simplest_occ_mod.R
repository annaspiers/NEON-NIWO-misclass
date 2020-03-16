# the simplest possible site-occupancy model

library(jagsUI)

# Select sample sizes (spatial and temporal replication) 
R <- 200 
T <- 3

# Determine process parameters 
psi <- 0.8 # Occupancy probability
p <- 0.5 # Detection probability

# Create structure to contain counts 
y <- matrix(NA, nrow = R, ncol = T) 

# Ecological process: Sample true occurrence (z, yes/no) from a Bernoulli (occurrence probability = psi)
z <- rbinom(n = R, size = 1, prob = psi) # Latent occurrence state

# Observation process: Sample detection/nondetection observations from a Bernoulli(with p) if z=1
for ( j in 1:T ) { 
    y[,j] <- rbinom(n = R, size = 1, prob = z * p) 
}

# Look at truth and at our imperfect observations 
sum(z) # Realized occupancy among 200 surveyed sites 
sum(apply(y, 1, max)) # Observed occupancy 

# Run JAGS model
str(JAGSdata <- list(y=y, T=T, R=R)) #bundle data
zst <- apply(y, 1, max) #observed occurrence as starting values for z
JAGSinits <- list(z=zst) 
JAGSparams <- c("psi", "p", "occ.fs") #params monitored
nc <- 1 #MCMC chains        #AIS, jags required the same number of chains as init values. only one init variable
ni <- 1200 #MCMC iterations
nb <- 200 #MCMC burnin
nt <- 2  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "jags_model.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)
names(jags_out)
plot(jags_out) #AIS why is plot of density of occ.fs not a smooth line?

# Resources
# http://www.mikemeredith.net/blog/2017/Occupancy_RNmodel.htm
