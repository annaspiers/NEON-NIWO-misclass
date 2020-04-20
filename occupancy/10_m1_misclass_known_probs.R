# In building up to a full misclassification model using NEON carabid data, we first build a simpler misclassificaiton model. See this resource for a full description: 
    # https://www.overleaf.com/project/5e938b1d9eb3590001f3e5e5


# MODEL 1: model with known misclassification probabilities

library(jagsUI)

# Problem 1 ---------------------------------------------------------------
# What will C_{i,j,k} be if a site is unoccupied (z_{i,k} = 0)?

# If occupancy is 0 (z_{i,k} = 0), then detection is 0 (y_{i,j,k} = 0). 
# The elements of the vector of counts, C_{i,j,k}, sum to y_{i,j,k}. 
# Thus, if z_{i,k} = 0, then C_{i,j,k} is a vector of 0's of length K



# Problem 2 ---------------------------------------------------------------
# Use the above model formulation to generate a simulated dataset with n = 20 sites, J = 5 visits, and K = 2 species. Assume that the confusion matrix is Θ = [[.8, .2], [.1, .9]], and choose reasonable values for the rest of the parameters ψ,z,λ

nsite <- 20 
nsurv <- 5
nspec <- 2
theta <- matrix(c(.8,.2,.1,.9), nrow=2, ncol=2, byrow=TRUE)
lambda_prob1 <- runif(1,0,1)
lambda_prob2 <- runif(1,0,1)

# one psi value per species
psi <- rbeta(nspec, 3,4)

# z dimensions: nsite rows, nspec columns
Z <- array(NA, dim=c(nrow=nsite, ncol=nspec))
for (i in 1:nsite) {
    Z[i,] <- rbinom(nspec, 1, psi)
}

# lambda dimensions: nsite rows, nsurv columns, nspec matrices
lambda <- array(NA, dim=c(nrow=nsite, ncol=nsurv, n3d=nspec ))
for (j in 1:nsurv) {
    for (i in 1:nsite) {
        lambda[i,j,] <- rbinom(nspec, 10, c(lambda_prob1,lambda_prob2)) #AIS increased size. Species 2 more abundant than species 1 when present
    }
}


# Simulate C, count of detections

# C dimensions: KxK, rows are correct species ID's, columns are (mis)classified species ID's - AIS right orientation?
C <- array(NA, dim=c(nrow=nspec, ncol=nspec, n3d=nsite, n4d=nsurv))
for (j in 1:nsurv) {
    for (i in 1:nsite) {
        for (k in 1:nspec) { # actual species
            for (k_class in 1:nspec) { # classified as 
                C[k,k_class,i,j] <- rpois(1, lambda = Z[i,k] * lambda[i,j,k] * theta[k,k_class])
            }
        }
    }
}


# Problem 3 ---------------------------------------------------------------
# Implement this model in JAGS, using your simulated data.
# How well can you recapture the true values of λ and ψ?

str(JAGSdata <- list(C=C, #AIS should we pass in C as data? if not, why do we do this in occ models?
                     nsite = nsite, 
                     nsurv = nsurv, 
                     nspec = nspec,
                     theta=theta)) #bundle data
z_init <- Z # AIS for the purposes of this model exercise, just allowed z_init to be Z
JAGSinits <- function(){list(z = z_init) } #rep(JAGSinits, nc) too
JAGSparams <- c("psi", "lambda") #params monitored
nc <- 4 #MCMC chains
ni <- 10000 #MCMC iterations
nb <- 3000 #MCMC burnin
nt <- 2  #MCMC thin

# JAGS model
jags_out <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/10_m1_misclass_known_probs_JAGS.txt", 
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(jags_out, dig=2)

# Ccompare true and recaptured psi 
jags_out$summary[1:2,]
psi

# Compare true and recaptured lambda 
hist(jags_out$summary[-c(1,2,nrow(jags_out$summary)),"mean"])
hist(lambda)

# AIS how to reconfigure model to make psi and lambda estimates better?