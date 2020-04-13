# This script is based off 05_static_multisp-y-array. In 05_... the detection matrix y sums detection values across all T surveys into one matrix. In this script, y is a 3D array with a detection matrix for each T survey.

library(jagsUI)

# Select sample sizes (spatial and temporal replication) 
R <- 200 
T <- 3
nspec <- 10


# Generate parameter values -----------------------------------------------

psi <- rbeta(nspec, 4, 1)
#hist(rbeta(10000, 4, 1))
psi

p <- rbeta(nspec, 1, 2)
#hist(rbeta(10000, 1, 2))
p

# Simulate true occupancy states
z <- matrix(nrow = R, ncol = nspec)
for (i in 1:R) {
    z[i, ] <- rbinom(nspec, 1, psi)
}
# z has closure throughout the season (does not vary across surveys)
image(t(z))
table(z)

# Simulate data -----------------------------------------------------------

# Make y a 3D array
y <- array(NA, dim=c(nrow=R, ncol=nspec, nmat=T))
for (k in 1:T) {
    for (i in 1:R) {
        y[i,,k] <- rbinom(nspec, 1, p * z[i, ])
    }
}
image(t(y[,,1]))
table(y[,,1])

# Fit model ---------------------------------------------------------------

str(JAGSdata <- list(y = y, 
                     T = T, 
                     R = R, 
                     nspec = nspec)) #bundle data
#z_init_test <- matrix((y > 0) * 1, nrow=R, ncol=nspec) #this only makes sense if y is a matrix (summed detection through surveys)
z_init <- apply(y, c(1,2), max, na.rm = TRUE)
z_init[z_init == "-Inf"] <- NA 
JAGSinits <- function(){list(z = z_init) } #rep(JAGSinits, nc) too
JAGSparams <- c("psi", "p") #params monitored
nc <- 3 #MCMC chains
ni <- 2000 #MCMC iterations
nb <- 500 #MCMC burnin
nt <- 2  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/06_static_multisp_y-array_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)
names(jags_out)
plot(jags_out) 


# Next, allow parameters to vary by sampling occassion (bootom of bpa ch 14 p.456)
# visualize these data across surveys(image or table or michael's sim visualizations)
