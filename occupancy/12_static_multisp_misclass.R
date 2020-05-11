# In building up to a full misclassification model using NEON carabid data, we first build a simpler misclassificaiton models. See this resource for a full description: 
    # https://www.overleaf.com/project/5e938b1d9eb3590001f3e5e5
# Here, we combine scripts 10_ and 11_ to create a static occupancy model estimating misclassification probabilities with simulated data. 

library(MCMCpack)
library(jagsUI)

# Define function to run JAGS model
jags_misclass_fn <- function(){
    str(JAGSdata <- list(nspec = nspec,
                         nsite = nsite, 
                         nsurv = nsurv, 
                         c_obs = c_obs,
                         alpha = alpha,
                         n = n,
                         M = M )) #bundle data
    JAGSinits <- function(){list(Z = Z, 
                                 theta = theta) }
    JAGSparams <- c("psi", "lambda", "theta", "Z") #params monitored
    nc <- 4 #MCMC chains
    #ni <- 20000 #MCMC iterations
    nb <- 4000 #MCMC burnin
    nt <- 1  #MCMC thin
    
    # JAGS model
    jags_out <- autojags(data = JAGSdata,
                     inits = JAGSinits,
                     parameters.to.save = JAGSparams,
                     model.file = "occupancy/12_static_multisp_misclass_JAGS.txt", 
                     n.chains = nc,
                     #n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     Rhat.limit = 1.05)
    return(jags_out)
}

# Simulate data for a static occupancy misclassification model with the following assumptions:
# There are K species. Theta and M_k are KxK matrices for classification probability and observed counts, respectively. Rows of theta are distributed as dirichlet random variables and rows of M_k are distributed as multinomial random variables. 
# Choose reasonable priors for theta that place higher probability mass on the parataxonomist getting the correct species ID than the incorrect one
# Use  realistic values for n (number of species k identified by expert taxonomist) and include cases where no expert data are available (n[k]=0 for species k)
# AIS This model assumes closure for Z and M throughout the season. 

nsite <- 10 
nsurv <- 4
nspec <- 3 


# Static occupancy model --------------------------------------------------
# The parameters below assume only a parataxonomist's classification is available

# Psi: occupancy probability
psi <- rbeta(nspec, 3,4) #one value per species

# Z: true occupancy at site i of species k (latent). dim = nsite rows, nspec columns
# Lambda: expected number of detections at an occupied site (latent). dim = nsite rows, nsurv columns, nspec matrices
# Y: actual number of detections at site i for species k on visit j, conditional on Z (latent)
Z       <- array(NA, dim = c(nsite=nsite, nspec=nspec))
lambda  <- array(NA, dim = c(nsite=nsite, nsurv=nsurv, nspec=nspec))
Y       <- array(NA, dim = c(nsite=nsite, nsurv=nsurv, nspec=nspec))
for (i in 1:nsite) {
    Z[i, ] <- rbinom(nspec, 1, psi)
    for (j in 1:nsurv) {
    	lambda[i,j, ] <- rgamma(nspec, 5, 1) #arbitrarily chosen
    	Y[i,j, ]      <- rpois(nspec, Z[i, ] * lambda[i,j, ]) 
	}
}

# Misclassification model -------------------------------------------------
# The parameters below assume both a parataxonomist's and an expert taxonomist's classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n <- round(runif(nspec,min=1,max=120)) # n values are need to be whole for multinom dist
#n[nspec] <- 0     #make some species unobserved by the expert taxonomist

# Alpha: dirichlet concentration parameter
alpha <- matrix(1,nrow=nspec,ncol=nspec)
diag(alpha) <- 15 #place higher probability mass on parataxonomist getting the correct classification

# Theta: misclassification probability matrix [KxK]
# M_k: M[k,k'] is the number of individuals from species k (according to expert) that were identified as species k' by parataxonomist [KxK]
M <- theta <- matrix(NA, nrow=nspec, ncol=nspec)
for (i in 1:nspec) {
    theta[i, ] <- MCMCpack::rdirichlet(1, alpha[i, ])
    M[i, ]     <- rmultinom(1, size=n[i], prob=theta[i, ])
}
image(theta, main="Simulated Theta")
image(M, main="Simulated M")


# Combine occupancy and misclassification models to simulate observed data --------

# C: C[i,j,k] is a vector of counts for species classified as species 1,...,K whose elements sum to Y[i,j,k]. dim = [KxK], rows are expert species ID's, columns are parataxonomist species ID's
# c_obs: c_obs[i,j,k'] are the elements of vector C, and represent the number of individuals that were classified as k'
c_obs   <- array(NA, dim = c(nsite=nsite, nsurv=nsurv, nspec=nspec))
C       <- array(NA, dim = c(dim(c_obs), nspec=nspec))
for (i in 1:nsite) {
    for (j in 1:nsurv) {
    	for (k in 1:nspec) { # actual species
        	C[i,j,k, ] <- rmultinom(1, Y[i,j,k], theta[k, ])
    	}
        for (k_prime in 1:nspec) {
            c_obs[i,j,k_prime] <- sum(C[i,j, ,k_prime]) #observed data
        }
	}
}


# JAGS model --------------------------------------------------------------

# Run model in JAGS. 
out <- jags_misclass_fn()

# How well does the model estimate specified parameters?
print(out, dig=2)
traceplot(out, parameters="psi")

# Compare true and recaptured psi 
psi
out$mean$psi

# Compare true and recaptured lambda 
plot(c(lambda), c(out$mean$lambda), main="Lambda") #we'll see 1-1 line 
abline(0,1)
for (i in 1:nrow(Z)) {
    for (k in 1:ncol(Z)) {
        points(c(lambda[i,,k]), c(out$mean$lambda[i,,k]), main="Lambda", 
             col=ifelse( round(out$mean$z[i,k])==0,"red","black")) 
    }
}
abline(0,1)

par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison")
abline(0,1)
for (i in 1:nspec) {
    points(x=out$mean$theta[i, ],
           y=theta[i,],
           col=i)
}

# Compare true and recaptured theta 
plot(c(theta), c(out$mean$theta), main="Theta")
abline(0,1)

# Compare true and recaptured Z 
plot(c(Z), c(out$mean$z), main="Theta")
abline(0,1)

# Visualize predictions of species unobserved by expert taxonomist
par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
abline(0,1)
points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")
