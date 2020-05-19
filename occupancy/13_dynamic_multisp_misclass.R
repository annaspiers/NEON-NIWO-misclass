# In building up to a full misclassification model using NEON carabid data, we first build a simpler misclassificaiton models. See this resource for a full description: 
    # https://www.overleaf.com/project/5e938b1d9eb3590001f3e5e5
# Here, we extend the static occupancy model with misclassification in script 12_ with simulated data. 

library(MCMCpack) #rdirchlet
library(jagsUI)

# Define function to run JAGS model
jags_misclass_fn <- function(){
    str(JAGSdata <- list(nspec = nspec,
                         nsite = nsite, 
                         nsurv = nsurv, 
                         nyear = nyear,
                         c_obs = c_obs,
                         alpha = alpha,
                         n = n,
                         M = M )) #bundle data
    JAGSinits <- function(){list(Z = Z, 
                                 theta = theta) }
    JAGSparams <- c("psi", "lambda", "theta", "Z", "p", "phi", "gamma", "n.occ", "growth", "turnover") #params monitored
    nc <- 4 #MCMC chains
    ni <- 8000 #MCMC iterations
    nb <- 4000 #MCMC burnin
    nt <- 1  #MCMC thin
    
    # JAGS model
    jags_out <- jags(data = JAGSdata,
                     inits = JAGSinits,
                     parameters.to.save = JAGSparams,
                     model.file = "occupancy/13_dynamic_multisp_misclass_JAGS.txt", 
                     n.chains = nc,
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt)
    return(jags_out)
}


# Simulate data for a static occupancy misclassification model with the following assumptions:
# There are K species. Theta and M_k are KxK matrices for classification probability and observed counts, respectively. Rows of theta are distributed as dirichlet random variables and rows of M_k are distributed as multinomial random variables. 
# Choose reasonable priors for theta that place higher probability mass on the parataxonomist getting the correct species ID than the incorrect one
# Use  realistic values for n (number of species k identified by expert taxonomist) and include cases where no expert data are available (n[k]=0 for species k)
# AIS This model assumes closure for Z and M throughout the season. 

nsite <- 8 
nsurv <- 5
nspec <- 2
nyear <- 4


# Dynamic occupancy model --------------------------------------------------
# The parameters below assume only a parataxonomist's classification is available

# Set initial occupancy and demographic parameters
# Psi: occupancy probability. dim: nspec 
# P:  detection probability. dim: nspec
# Phi:  survival probability. dim: nspec 
# Gamma:  colonization probability. dim: nspec 
psi     <- rbeta(nspec, 3, 4) 
p       <- rbeta(nspec, 3, 4)  #chose to be somewhere in the middle %
phi     <- rbeta(nspec, 8, 2)  #chose to be high %
gamma   <- rbeta(nspec, 1, 15) #chose to be low %


# Z: true occupancy at site i of species k (latent). dim = nsite, nspec, nyear
Z <- array(NA, dim = c(nsite=nsite, nspec=nspec, nyear=nyear))
for (i in 1:nsite) {
    Z[i, ,1] <- rbinom(nspec, 1, psi)
    for (l in 2:nyear) {
        Z[i, ,l] <- rbinom(nspec, 1, Z[i, ,l-1]*phi + (1-Z[i, ,l-1])*gamma)
	}
}

# Lambda: expected number of detections at an occupied site (latent). dim = nsite, nsurv, nspec, nyear
# Y: actual number of detections at site i for species k on visit j, conditional on Z (latent)
lambda  <- array(NA, dim = c(nsite=nsite, nsurv=nsurv, nspec=nspec, nyear=nyear)) 
Y       <- array(NA, dim = dim(lambda))
for (i in 1:nsite) {
    for (j in 1:nsurv) {
        for (l in 1:nyear) {
            lambda[i,j, ,l] <- rgamma(nspec, 5, 1) #arbitrarily chosen
    	    Y[i,j, ,l]      <- rpois(nspec, Z[i, ,l] * lambda[i,j, ,l]) 
        }
	}
}

# Misclassification model -------------------------------------------------
# The parameters below assume both a parataxonomist's and an expert taxonomist's classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n <- round(runif(nspec,min=1,max=120)) # n values are need to be whole for multinom dist
#n[nspec] <- 0     #make some species unobserved by the expert taxonomist


# Alpha: dirichlet concentration parameter
alpha <- matrix(1, nrow=nspec, ncol=nspec)
diag(alpha) <- 15 #place higher probability mass on parataxonomist getting the correct classification

# Theta: misclassification probability matrix [KxK]
# M_k: M[k,k'] is the number of individuals from species k (according to expert) that were identified as species k' by parataxonomist [KxK]
M <- theta <- matrix(NA, nrow=nspec, ncol=nspec)
for (k in 1:nspec) {
    theta[k, ] <- MCMCpack::rdirichlet(1, alpha[k, ])
    M[k, ]     <- rmultinom(1, size=n[k], prob=theta[k, ])
}


# Combine occupancy and misclassification models to simulate observed data --------

# C: C[i,j,k,l] is a vector of counts for species classified as species 1,...,K whose elements sum to Y[i,j,k]. dim = [KxK], rows are expert species ID's, columns are parataxonomist species ID's
# c_obs: c_obs[i,j,k',l] are the elements of vector C, and represent the number of individuals that were classified as k'
c_obs   <- array(NA, dim = c(nsite=nsite, nsurv=nsurv, nspec=nspec, nyear=nyear))
C       <- array(NA, dim = c(dim(c_obs), nspec=nspec))
for (i in 1:nsite) {
    for (j in 1:nsurv) {
        for (l in 1:nyear) {
            for (k in 1:nspec) { # actual species
        	    C[i,j,k,l, ] <- rmultinom(1, Y[i,j,k,l], theta[k, ])
    	    }
            for (k_prime in 1:nspec) {
                c_obs[i,j,k_prime,l] <- sum(C[i,j, ,l,k_prime]) #observed data
            }
        }
	}
}


# JAGS model --------------------------------------------------------------

# Run model in JAGS. 
out <- jags_misclass_fn()

# How well does the model estimate specified parameters?
print(out, dig=2)
traceplot(out, parameters="Z")

# Compare true and recaptured psi 
psi
out$mean$psi

# Compare true and recaptured lambda 
plot(NA,xlim=c(0,max(lambda)),ylim=c(0,max(lambda)),
     xlab="Observed",ylab="Predicted",main="Lambda comparison")
abline(0,1)
for (i in 1:dim(Z)[1]) {
    for (k in 1:dim(Z)[2]) {
        for (l in 1:dim(Z)[3])
            points(c(lambda[i,,k,l]), c(out$mean$lambda[i,,k,l]), main="Lambda", 
                   col=ifelse( round(out$mean$Z[i,k,l])==0,"red","black")) 
    }
}


# Compare true and recaptured theta 
plot(c(theta), c(out$mean$theta), main="Theta")
abline(0,1)

# Compare true and recaptured Z 
plot(c(Z), c(out$mean$Z), main="Z")
abline(0,1)

# Visualize predictions of species unobserved by expert taxonomist
#par(mfrow=c(1,1))
#plot(NA,xlim=c(0,1),ylim=c(0,1),
#     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
#abline(0,1)
#points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")

# Plot apparent occupancy 
    psi.app <- apply(apply(y, c(1,3), max), 2, mean) 
    lines(year, psi.app, type = "l", col = "black", lwd = 2) 
    text(0.85*K, 0.06, labels = "red solid – true occupancy\n red dashed – detection probability\n black – observed occupancy")