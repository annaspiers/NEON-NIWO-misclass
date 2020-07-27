# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a dynamic occupancy 
# model estimating misclassification probabilities with simulated data

library(reshape2) 
library(tidyverse) 
library(gtools) 
library(R2jags)
library(MCMCpack) #rdirchlet

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
Z <- array(dim = c(nsite=nsite, nspec=nspec, nyear=nyear))
for (i in 1:nsite) {
    Z[i, ,1] <- rbinom(nspec, 1, psi)
    for (y in 2:nyear) {
        Z[i,  , y] <- rbinom(nspec, 1, Z[i, ,y-1]*phi + (1-Z[i, ,y-1])*gamma)
	}
}

# Lambda: expected number of detections at an occupied site (latent). dim = nsite, nsurv, nspec, nyear
# Y: actual number of detections at site i for species k on visit j, conditional on Z (latent)
lambda  <- array(NA, dim = c(nspec=nspec, nyear=nyear)) 
for (y in 1:nyear) {
    lambda[ ,y] <- exp(rnorm(nspec, mean =1, sd = 1))
}

# Misclassification model -------------------------------------------------
# The parameters below assume both a parataxonomist's and an expert taxonomist's classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n <- array(dim = c(nsite, nsurv, nspec, nyear))
for (i in 1:nsite) {
    for (j in 1:nsurv) {
        for (k in 1:nspec) {
            for (y in 1:nyear) {
                n[i, j, k, y] <- rpois(1, Z[i, k, y] * lambda[k, y])
            }
        }
    }
}

L <- apply(n, c(1, 2, 4), sum)

# Now disaggregate the encounter counts to generate a record for each encounter.
n_df <- reshape2::melt(n, varnames = c("site", "survey", "species","year"))
# each row in k is an individual beetle. row l marks the individual index for
# each unique combo of sitexsurveyxspeciesxyear 
k_df <- n_df %>%
    group_by(site, survey, species, year) %>%
    summarize(l = list(seq_len(value))) %>%
    unnest(l) %>%
    ungroup

## Simulating imperfect species classifications Now, generate a probability
#vector `y` for each of these detections. Assume that we have a noisy
#classifier, and that the skill of the classifier might vary by species.
alpha <- matrix(1, nrow = nspec, ncol = nspec)
diag(alpha) <- 10

Theta_true <- matrix(nrow = nspec, ncol = nspec)
for (k in 1:nspec) {
    Theta_true[k, ] <- rdirichlet(1, alpha[k, ])
}
print(Theta_true)

noisy_classifier <- function(true_species, Theta_true) {
    sample(nspec, size = length(true_species), replace = TRUE, 
           prob = Theta_true[true_species, ])
}

# get probabilities from the classifier
y_df <- k_df %>%
    rowwise %>%
    mutate(y = list(noisy_classifier(species, Theta_true))) %>%
    ungroup %>%
    unnest(y) %>%
    mutate(idx = 1:n())
#AIS here, species variable is the true ID while the y variable is the 
# imperfectly classified ID?

# Simulating "ground truth" data
# Assume that we have a subset of the data with known species IDs. 
pct_known_species <- .3
y_df <- y_df %>%
    mutate(true_species_known = rbinom(n(), size = 1, prob = pct_known_species))

# How many labeled examples do we have of each species?
y_df %>%
    group_by(species) %>%
    summarize(n_known_species_IDs = sum(true_species_known))

## Model fitting
# Fit the model and check convergence.
jags_d <- list(nsite = nsite, 
               K = nspec, 
               nsurv = nsurv, 
               nyear = nyear,
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = k_df$site, 
               # if the individual was labeled by the expert, true ID is known
               k = ifelse(y_df$true_species_known, y_df$species, NA),
               # for all individuals, we get paratxonomist IDs
               y = y_df$y)

init_fn_factory <- function(nsite, nspec, nyear) {
    function(){
        list(z = matrix(1, nsite, nspec, nyear))
    } 
}

ji <- init_fn_factory(nsite = nsite, nspec = nspec, nyear)

# initialize the function
inits <- ji()

jm <- jags.parallel(
    data = jags_d, 
    inits = ji, 
    parameters.to.save = c("psi", "p", "phi", "gamma", "lambda", 
                           "Theta", "pi", "k", "y"), 
    model.file = "occupancy/13.1_dynamic_multisp_misclass_JAGS.txt", 
    n.chains = 6, 
    n.iter = 10000, 
    n.thin = 1,
    DIC = FALSE)

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

# Plot apparent occupancy 
    psi.app <- apply(apply(y, c(1,3), max), 2, mean) 
    lines(year, psi.app, type = "l", col = "black", lwd = 2) 
    text(0.85*K, 0.06, 
         labels = "red solid – true occupancy\n red dashed – detection probability\n black – observed occupancy")