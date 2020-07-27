# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a dynamic occupancy 
# model estimating misclassification probabilities with simulated data

library(reshape2) 
library(tidyverse) 
library(gtools) 
library(R2jags)
library(MCMCpack) #rdirchlet
library(MCMCvis)

set.seed(1234567)
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
               year = k_df$year,
               # if the individual was labeled by the expert, true ID is known
               k = ifelse(y_df$true_species_known, y_df$species, NA),
               # for all individuals, we get paratxonomist IDs
               y = y_df$y)

init_fn_factory <- function(nsite, nspec, nyear) {
    function(){
        list(Z = array(1, dim = c(nsite, nspec, nyear)))
    } 
}

ji <- init_fn_factory(nsite = nsite, nspec = nspec, nyear = nyear)

# initialize the function
inits <- ji()

jm <- jags.parallel(
    data = jags_d, 
    inits = ji, 
    parameters.to.save = c("psi", "p", "phi", "gamma", "lambda", 
                           "Theta", "pi"), 
    model.file = "occupancy/13.1_dynamic_multisp_misclass_JAGS.txt", 
    n.chains = 6, 
    n.iter = 10000, 
    n.thin = 1,
    DIC = FALSE)

# How well does the model estimate specified parameters?
jm_summ <- MCMCsummary(jm, dig=2)
hist(jm_summ$Rhat, breaks=40)
range(jm_summ$Rhat, na.rm=TRUE)

# Look at high Rhat values
jm_summ <- rownames_to_column(jm_summ)
jm_summ %>% filter(Rhat > 3)

# Compare true vs recaptured values
# psi
MCMCtrace(jm, params="Theta", type = 'both', ind = F, pdf=F, ISB=F, Rhat=T)

# p
# gamma
# lambda
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

# theta
MCMCtrace(jm, params="Theta", type = 'both', ind = F, pdf=F, ISB=F, Rhat=T)
plot(c(Theta_true), c(Theta_summ$mean$Theta), main="Theta",
     xlim=c(0,1), ylim=c(0,1))
abline(0,1)

# pi
pi_summ <- MCMCsummary(jm, params="pi", dig=2)
