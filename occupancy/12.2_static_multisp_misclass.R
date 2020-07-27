# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a static occupancy 
# (single season) model estimating misclassification probabilities with 
# NEON data.

library(reshape2) 
library(tidyverse) 
library(gtools) 
library(R2jags)

all_paratax_df <- readRDS("occupancy/all_paratax_df.rds")

# Create df for disaggregated individual model. 
# Modify bet_sorting such that each row is an individual beetle
all_paratax_by_ind <- all_paratax_df %>%
    uncount(individualCount) %>%
    rownames_to_column()

# Filter to one year for static model
all_paratax_by_ind_2018 <- all_paratax_by_ind %>%
    filter(col_year == "2018")
    
    
set.seed(1234567)
nsite <- 15
nspecies <- 2
noccasion <- 3

psi <- runif(nspecies)
lambda <- exp(rnorm(nspecies, mean =1, sd = 1))

z <- matrix(nrow = nsite, ncol = nspecies)
for (i in 1:nsite) {
    for (k in 1:nspecies) {
        z[i, k] <- rbinom(1, 1, psi[k])
    }
}

n <- array(dim = c(nsite, noccasion, nspecies))
for (i in 1:nsite) {
    for (j in 1:noccasion) {
        for (k in 1:nspecies) {
            n[i, j, k] <- rpois(1, z[i, k] * lambda[k])
        }
    }
}

L <- apply(n, c(1, 2), sum)

# Now disaggregate the encounter counts to generate a record for each encounter.
n_df <- reshape2::melt(n, varnames = c("site", "occasion", "species"))
#AIS ^ this is where I'll need to modify bet_sorting so each row is an
#individual. Right? but I thought n was true species, not observed species like
#bet_sorting represents
k_df <- n_df %>%
    group_by(site, occasion, species) %>%
    summarize(l = list(seq_len(value))) %>%
    unnest(l) %>%
    ungroup

## Simulating imperfect species classifications Now, generate a probability
#vector `y` for each of these detections. Assume that we have a noisy
#classifier, and that the skill of the classifier might vary by species.
alpha <- matrix(1, nrow = nspecies, ncol = nspecies)
diag(alpha) <- 10

Theta_true <- matrix(nrow = nspecies, ncol = nspecies)
for (k in 1:nspecies) {
    Theta_true[k, ] <- rdirichlet(1, alpha[k, ])
}
print(Theta_true)

noisy_classifier <- function(true_species, Theta_true) {
    sample(nspecies, size = length(true_species), replace = TRUE, 
           prob = Theta_true[true_species, ])
}

# get probabilities from the classifier
y_df <- k_df %>%
    rowwise %>%
    mutate(y = list(noisy_classifier(species, Theta_true))) %>%
    ungroup %>%
    unnest(y) %>%
    mutate(idx = 1:n())

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
               K = nspecies, 
               noccasion = noccasion, 
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = k_df$site, 
               # if the individual was labeled by the expert, true ID is known
               k = ifelse(y_df$true_species_known, y_df$species, NA),
               # for all individuals, we get paratxonomist IDs
               y = y_df$y)

init_fn_factory <- function(nsite, nspecies) {
    function(){
        list(z = matrix(1, nsite, nspecies))
    } 
}

ji <- init_fn_factory(nsite = nsite, nspecies = nspecies)

# initialize the function
inits <- ji()

jm <- jags.parallel(
    data = jags_d, 
    inits = ji, 
    parameters.to.save = c("psi", "lambda", "Theta"), 
    model.file = "occupancy/disagg-ind-mod.txt", 
    n.chains = 6, 
    n.iter = 100000, 
    n.thin = 1,
    DIC = FALSE)

# Let's look at the trace plots.
traceplot(jm, 
          ask = FALSE, 
          mfrow = c(4, 2), 
          varname = c("psi", "lambda", "Theta"))

# Check convergence
jm_summ <- print(jm, dig=2)
par(mfrow=c(1,1))
hist(jm_summ$summary[,"Rhat"])

## Compare true vs recaptured values -------
# Theta - compare true and recaptured
plot(c(Theta_true), c(jm_summ$mean$Theta), main="Theta",
     xlim=c(0,1), ylim=c(0,1))
abline(0,1)

# lambda - compare true and recaptured  
traceplot(jm, varname="lambda")
plot(c(lambda), c(jm_summ$mean$lambda), main="Lambda",
     xlim=c(0,max(lambda)), ylim=c(0,max(lambda))) 
abline(0,1)
# for (i in 1:nrow(Z)) {
#     for (k in 1:ncol(Z)) {
#         points(c(lambda[i,,k]), c(out$mean$lambda[i,,k]), main="Lambda", 
#                col=ifelse( round(out$mean$z[i,k])==0,"red","black")) 
#     }
# }
# abline(0,1)

# psi - compare true and recaptured  
traceplot(jm, varname="psi")
plot(c(psi), c(jm_summ$mean$psi), main="Psi", xlim=c(0,1),ylim=c(0,1))
abline(0,1)
