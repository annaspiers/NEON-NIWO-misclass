# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a dynamic occupancy 
# model estimating misclassification probabilities with simulated data

library(reshape2) 
library(tidyverse) 
library(gtools) 
library(R2jags)
library(MCMCpack) #rdirchlet
library(MCMCvis)

nsite <- 20 
nsurv <- 2
nspec <- 2
nyear <- 4

simulate_data <- function() {
  mu <- rnorm(4)
  Sigma_epsilon <- MCMCpack::riwish(10, diag(4))
  Sigma_alpha <- MCMCpack::riwish(10, diag(4))
  
  alpha <- t(chol(Sigma_alpha)) %*% matrix(rnorm(nspec * 4), nrow = 4)
  epsilon <- t(chol(Sigma_epsilon)) %*% matrix(rnorm(nsite * 4), nrow = 4)
  
  psi1 <- matrix(nrow = nsite, ncol = nspec)
  phi <- matrix(nrow = nsite, ncol = nspec)
  gamma <- matrix(nrow = nsite, ncol = nspec)
  lambda <- matrix(nrow = nsite, ncol = nspec)
  
  for (i in 1:nsite) {
    for (k in 1:nspec) {
      psi1[i, k] <- plogis(mu[1] + alpha[1, k] + epsilon[1, i])
      phi[i, k] <- plogis(mu[2] + alpha[2, k] + epsilon[2, i])
      gamma[i, k] <- plogis(mu[3] + alpha[3, k] + epsilon[3, i])
      lambda[i, k] <- exp(mu[4] + alpha[4, k] + epsilon[4, i])
    }
  }
  
  # z: true occupancy at site i of species k (latent). dim = nsite, nspec, nyear
  z <- array(dim = c(nsite=nsite, nspec=nspec, nyear=nyear))
  for (i in 1:nsite) {
    z[i, ,1] <- rbinom(nspec, 1, psi1[i, ])
    for (y in 2:nyear) {
      z[i,  , y] <- rbinom(nspec, 
                           1, 
                           z[i, ,y-1]*phi[i, ] + 
                             (1-z[i, ,y-1])*gamma[i, ])
    }
  }
  
  # Encounter model -------------------------------------------------
  L <- array(dim = c(nsite, nsurv, nspec, nyear))
  for (i in 1:nsite) {
    for (j in 1:nsurv) {
      for (k in 1:nspec) {
        for (t in 1:nyear) {
          L[i, j, k, t] <- rpois(1, z[i, k, t] * lambda[i, k])
        }
      }
    }
  }
  
  Ldot <- apply(L, c(1, 2, 4), sum)
  
  # Now disaggregate the encounter counts to generate a record for each encounter.
  n_df <- reshape2::melt(L, varnames = c("site", "survey", "species","year"))
  # each row in k is an individual beetle. row l marks the individual index for
  # each unique combo of sitexsurveyxspeciesxyear 
  k_df <- n_df %>%
    group_by(site, survey, species, year) %>%
    summarize(l = list(seq_len(value)), .groups = "keep") %>%
    unnest(l) %>%
    ungroup
  
  ## Simulating imperfect species classifications Now, generate a probability
  #vector `y` for each of these detections. Assume that we have a noisy
  #classifier, and that the skill of the classifier might vary by species.
  a <- matrix(1, nrow = nspec, ncol = nspec)
  diag(a) <- 20
  
  Theta_true <- matrix(nrow = nspec, ncol = nspec)
  for (k in 1:nspec) {
    Theta_true[k, ] <- rdirichlet(1, a[k, ])
  }
  print(Theta_true)

  noisy_classifier <- function(true_species, Theta_true) {
    sample(nspec, size = length(true_species), replace = TRUE, 
           prob = Theta_true[true_species, ])
  }
  
  # get imperfect classes from the classifier
  y_df <- k_df %>%
    rowwise %>%
    mutate(y = list(noisy_classifier(species, Theta_true))) %>%
    ungroup %>%
    unnest(y) %>%
    mutate(idx = 1:n())
  # here, species variable is the true ID while the y variable is the 
  # imperfectly classified ID?
  
  # Simulating "ground truth" data
  # Assume that we have a subset of the data with known species IDs. 
  pct_known_species <- .3
  y_df <- y_df %>%
    mutate(true_species_known = rbinom(n(), size = 1, prob = pct_known_species))
  
  list(
    nsite = nsite, 
    K_exp = nspec, 
    K_para = nspec,
    nsurv = nsurv, 
    nyear = nyear,
    L = Ldot, 
    a = a,
    R = diag(4),
    Ltot = sum(L), 
    site = k_df$site, 
    year = k_df$year,
    # if the individual was labeled by the expert, true ID is known
    k = ifelse(y_df$true_species_known, y_df$species, NA),
    # for all individuals, we get paratxonomist IDs
    y = y_df$y, 
    Theta_true = Theta_true, 
    Tau_spec_true = solve(Sigma_alpha), 
    Tau_site_true = solve(Sigma_epsilon), 
    eps_site_true = t(epsilon), 
    eps_spec_true = t(alpha), 
    mu_true = mu)
}


init_fn_factory <- function(nsite, nspec, nyear) {
  function(){
    list(z = array(1, dim = c(nsite, nspec, nyear)))
  } 
}

ji <- init_fn_factory(nsite = nsite, nspec = nspec, nyear = nyear)

# initialize the function
inits <- ji()

sbc_sim <- function(x) {
  d <- simulate_data()
  jm <- R2jags::jags.parallel(
    data = d, 
    inits = ji, 
    parameters.to.save = c("ranks_Theta", "ranks_Tau_site", "ranks_Tau_spec", 
                           "ranks_eps_site", "ranks_eps_spec", "ranks_mu"), 
    model.file = "neon_dynamic_disagg_multisp_sbc.txt", 
    n.chains = 4, 
    n.iter = 5000, 
    DIC = FALSE)
  # return a n_iter X n_param matrix of ranks
  jm$BUGSoutput$sims.array %>%
    apply(3, c)
}
ranks <- pbapply::pblapply(1:1000, sbc_sim)

dir.create("output", showWarnings = FALSE)
write_rds(ranks, "output/sbc_ranks.rds")


plot_sbc <- function(ranks, thin = 3, ...) {
  thinner <- seq(from = 1, to = nrow(ranks[[1]]), by = thin)
  u <- t(sapply(ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
  parameter <- as.factor(rep(colnames(u), each = nrow(u)))
  parameter <- gsub("ranks_", "", parameter)
  d <- data.frame(u = c(u), parameter)
  suppressWarnings(ggplot2::ggplot(d) + 
                     ggplot2::geom_histogram(ggplot2::aes(x = u), ...) + 
                     ggplot2::facet_wrap("parameter"))
}

plot_sbc(ranks, bins = 20) + 
  ylab("Count")
ggsave("figures/sbc.pdf", width=10, height = 7)
