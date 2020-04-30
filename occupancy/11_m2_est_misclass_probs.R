# In building up to a full misclassification model using NEON carabid data, we first build a simpler misclassificaiton model. See this resource for a full description: 
    # https://www.overleaf.com/project/5e938b1d9eb3590001f3e5e5
# nice explanation of dirichlet and beta distributions: https://stats.stackexchange.com/questions/244917/what-exactly-is-the-alpha1-in-the-dirichlet-distribution

# MODEL 2: estimate misclassification probabilities

library(MCMCpack)
library(jagsUI)

# Problem 4 ---------------------------------------------------------------
# Experiment with visualizing the prior over Θ1 for some different Dirichlet distributions, including α1=  (1,1), αk=  (10,1), αk=(20,2).

vis_theta_dens <- function(alphak ) { # alpha is Dirichlet parameter vector
    thetak <- MCMCpack::rdirichlet(10000, alphak) # simulate 10,000 draws from the prior for thetak
    # plot prior probability of a correct classification
    plot(density(thetak[, 1], from = 0, to = 1), xlim = c(0, 1)) 
    lines(density(thetak[, 2], from = 0, to = 1), col = "red") # add prior probability of incorrect classif.
}

vis_theta_dens(alphak = c(1, 1))
vis_theta_dens(alphak = c(10, 10))
vis_theta_dens(alphak = c(100, 100)) 
#with equal alpha elements, larger magnitudes make a narrower density peak

vis_theta_dens(alphak = c(10, 1))
vis_theta_dens(alphak = c(50, 5))
vis_theta_dens(alphak = c(100, 10))
vis_theta_dens(alphak = c(1000, 100)) 
#element ratio of 10:1 centers peak over 0.9, larger magnitudes makes narrower peak

vis_theta_dens(alphak = c(5, 4))
vis_theta_dens(alphak = c(20, 19))
vis_theta_dens(alphak = c(100, 99))
# with element1 > element2 by 1, larger magnitudes approach equal elements

vis_theta_dens(alphak = c(10, 5))
vis_theta_dens(alphak = c(50, 25))
vis_theta_dens(alphak = c(200, 100))
# with element1 = 2 * element2, peaks center over 1/3 and 2/3

# Problem 5 ---------------------------------------------------------------
# Simulate data from this misclassification model for Theta and M_k (both dim=KxK) for K= 3 species, using the following prior distributions
# theta1 ~ Dirichlet(10,1,1)
# theta2 ~ Dirichlet(1,10,1)
# theta3 ~ Dirichlet(1,1,10)
# Assume  n1=25, n2=16, and n3=5, where expert taxonomist ID's nk individuals of species k. 

nspec <- 3 

n1 <- 25
n2 <- 16
n3 <- 5

# Define dirichlet parameter vectors
alpha1 <- c(10,1,1)
alpha2 <- c(1,10,1)
alpha3 <- c(1,1,10)

# Simulate rows of theta with dirichlet priors
theta1 <- MCMCpack::rdirichlet(1, alpha1)
theta2 <- MCMCpack::rdirichlet(1, alpha2)
theta3 <- MCMCpack::rdirichlet(1, alpha3)

# Define (mis)classification probability matrix
theta <- matrix(c(theta1,theta2,theta3), nrow=nspec, ncol=nspec, byrow=TRUE)
theta
# So there is a 75.5% chance that species 1 was classified as species 1, 20.8% it was misclassified as species 2, and 3.7% chance it was misclassified as species 3

# Simulate observation matrix M
M <- matrix(NA, nrow=nspec, ncol=nspec)
M[1,] <- rmultinom(1,size=n1,prob=theta1)
M[2,] <- rmultinom(1,size=n2,prob=theta2)
M[3,] <- rmultinom(1,size=n3,prob=theta3)


# Problem 6 ---------------------------------------------------------------
# Implement this model in JAGS. Model the priors for rows in Theta and a likelihood for observation matrix M. Sample from the posterior using the data you simulated in the last problem. How well could you estimate the true values of Θ?  Did you get better estimates for species with larger sample sizes?

jags_misclass_fn <- function(M_arg, nspec_arg, theta_arg){
    str(JAGSdata <- list(M=M_arg, 
                     nspec=nspec_arg)) #bundle data
    JAGSinits <- function(){list(theta=theta_arg) } #rep(JAGSinits, nc) too
    JAGSparams <- c("theta") #params monitored
    nc <- 4 #MCMC chains
    ni <- 7000 #MCMC iterations
    nb <- 1000 #MCMC burnin
    nt <- 2  #MCMC thin
    
    # JAGS model
    jags_out <- jags(data = JAGSdata,
                     inits = JAGSinits,
                     parameters.to.save = JAGSparams,
                     model.file = "occupancy/11_m2_est_misclass_probs_JAGS.txt", 
                     n.chains = nc,
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt)
    return(jags_out)
}
out <- jags_misclass_fn(M_arg=M, nspec_arg=nspec, theta_arg=theta)
print(out, dig=2)
plot(out)

# AIS why is deviance 0?

# Compare true and recaptured theta 
theta
out$mean$theta
# The model underestimates the true species k-to-k' probabilities and overestimates the species k-to-k' misclassifications.
# AIS can anything be done to the model to make these more accurate? Add more data?

plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison")
abline(0,1)
for (i in 1:nspec) {
    points(x=out$mean$theta[i,],
         y=theta[i,],
         col=c("red","blue","green")[i])
}

# Species 1 had the largest sample size, but did not necessarily yield the best prediction of all the species. It appears that species 2 had the best predictions (i.e. closest to 1:1 line in plot). However, increasing the sample size of each species my a couple orders of magnitude yielded more accurate estimates


# Problem 7 ---------------------------------------------------------------
# Simulate data for K= 20 species, choosing reasonable priors for Θ (that place higher probability mass on the parataxonomist getting the correct species ID than the incorrect one), and whatever values you want for n1,...,nK. Implement this model in JAGS. How well did the model estimate Θ?

nspec <- 20 

n <- runif(nspec,min=1,max=50) #AIS n values are not whole

# Define dirichlet parameter vectors
alpha <- matrix(1, nrow=nspec, ncol=nspec)
diag(alpha) <- 20

# Simulate misclassification prob matrix, theta, and observation matrix, M 
M <- theta <- matrix(NA, nrow=nspec, ncol=nspec)
for (i in 1:nspec) {
    theta[i,] <- MCMCpack::rdirichlet(1, alpha[i,])
    M[i,] <- rmultinom(1,size=n[i],prob=theta[i,])
}
image(theta)
image(M)

out <- jags_misclass_fn(M_arg=M, nspec_arg=nspec, theta_arg=theta)
print(out, dig=2)
plot(out)

# Compare true and recaptured theta 
par(mfrow=c(1,2))
image(theta, main="Observed theta")
image(out$mean$theta, main="Predicted theta")

par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison")
abline(0,1)
for (i in 1:nspec) {
    points(x=out$mean$theta[i,],
         y=theta[i,],
         col=i)
}
# Prerdicts misclassification probs pretty well, but underpredicts correct identifications



# Problem 8 ---------------------------------------------------------------
# Check your implementation to see whether it still works in caseswhere no expert data are available (nk=0 for some speciesk).
n[c(5,10)] <- 0

M <- theta <- matrix(NA, nrow=nspec, ncol=nspec)
for (i in 1:nspec) {
    theta[i,] <- MCMCpack::rdirichlet(1, alpha[i,])
    M[i,] <- rmultinom(1,size=n[i],prob=theta[i,])
}
image(theta)
image(M) # no abundance for M[5,] and M[10,]

out <- jags_misclass_fn(M_arg=M, nspec_arg=nspec, theta_arg=theta)
print(out, dig=2)
plot(out)

# Does the model correctly model theta to show uniform probability across species with no expert ID? This is what we'd expect for a species with no expert ID
par(mfrow=c(1,2))
image(theta, main="Observed theta")
image(out$mean$theta, main="Predicted theta")

par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
abline(0,1)
points(x=out$mean$theta[5,],y=theta[5,],col="red")
points(x=out$mean$theta[10,],y=theta[10,],col="blue")
# Cool!
