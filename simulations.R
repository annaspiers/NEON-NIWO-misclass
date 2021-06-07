# Simulation varying fraction of samples validated 

# This simulation demonstrates the model behavior when varying fractions of
# samples are validated at the individual level. Data are generated using model
# priors. We implement the simulations using a single-season version of the featured
# classification-occupancy model. This simulation illustrates model behavior on
# a generic dataset where all samples (e.g., individuals, states) are
# imperfectly classified so as to reach a broader audience of ecologists that
# are generally interested in classification problems, regardless of
# application. 

# AIS revisit the following to make sure we actually do these:
# We show how this model improves ecological estimation/inference by
# illustrating how close the estimated occupancy or demographic rates are to the
# true parameter values. We also show what fraction of data need to be validated to
# achieve certain model accuracy.

# Model assumptions
# Occupancy has closure throughout the season

library(dclone) #jags.parfit()
library(MCMCvis)
library(dplyr)
library(ggmcmc)
library(tibble) #rownames_to_column
library(MASS) #mvrnorm()
library(boot) #inv.logit()
library(MCMCpack) #rdirichlet()
library(extraDistr) #rcat()
options(digits = 3)

# Simulate data from model priors
# By generating data from the model, fitting the model to these data, then computing posteriors for (occupancy, demographic rate) parameter estimates, we have theoretical guarantees about interval coverage. 

# Generate data from model priors --------------------------------------------------

set.seed(201023923)

nsite <- 10 #number of sites
nsurv <- 4  #number of suveys in the season
K_imp <- 8  #number of states identified by imperfect classifier
K_val <- 5  #number of states identified in validation

# Simulate 
#n <- 1000 #AIS what should the sample size be, maybe 10,000 in the end?
#df <- 10 #degrees of freedom
R <- diag(rep(1, 2)) # 2 dim since we're estimating psi and lambda

mu_spec <- rnorm(2, mean=0, sd=1)
mu_site <- rep(0,2) 
  
#Tau_spec <- Tau_site <- rWishart(1, df, R) #AIS what am I doing wrong where I can't use tau as an argument in mvrnorm?
eps_spec <- mvrnorm(n=K_val, mu=mu_spec, Sigma=R) 
eps_site <- mvrnorm(n=nsite, mu=mu_site, Sigma=R)
  

### Simulate data for occupancy model
# The parameters below assume only imperfect classifications are available

# Psi (occupancy probability) and lambda (encounter rate)
psi <- lambda <- matrix(NA,nrow=nsite,ncol=K_val)
for (i in 1:nsite){
    for (k in 1:K_val) {
      psi[i, k] <- inv.logit(eps_site[i, 1] + eps_spec[k, 1])
      lambda[i, k] <- exp(eps_site[i, 2] + eps_spec[k, 2])
    }
  }

# Theta (classification probability)
alpha <- matrix(2, nrow=K_val, ncol=K_imp)
diag(alpha) <- 200
Theta <- matrix(NA, nrow=K_val, ncol=K_imp)
for (k in 1:K_val) {
  Theta[k,1:K_imp] <- rdirichlet(1,alpha[k,1:K_imp])
}

# Z (expected occupancy) and zlam (expected abundance given occupancy)
z <- zlam <- matrix(NA, nrow=nsite, ncol=K_val)
for (i in 1:nsite) {
  for (k in 1:K_val) {
    z[i, k] <- rbinom(1, 1, psi[i, k]) #bernoulli trial
    zlam[i, k] <- z[i, k] * lambda[i, k]
  }
}

# L (expected number of encountered individuals)
L <- matrix(NA, nrow=nsite, ncol=nsurv)
for (i in 1:nsite) {
  for (j in 1:nsurv) {
    L[i, j] <- rpois(1, sum(zlam[i, 1:K_val]) ) #sum zlam across each row for number of imperfectly detected individuals
  }
}

# Pi (dirichlet distribution for each encountered individual)
# k (validated species for an individual)
# y (imperfectly observed state for an individual)
Ltot <- sum(L)
site <- sample(1:nsite, Ltot, replace=TRUE) #AIS alternative ways to initialize site?
pi <- matrix(NA, nrow=Ltot, ncol=K_val)
k <- y <- rep(NA, Ltot)
for (l in 1:Ltot) {
  pi[l, 1:K_val] <- zlam[site[l], 1:K_val] / sum(zlam[site[l], 1:K_val])
  if (is.na(pi[l,1])) { #arbitrary column, all in a row would be NA
    pi[l, ] <- rep(1/K_val,K_val)
  } #AIS is this an ok workaround? When pi has a row of NA's, then k is NA for that row, which we don't want. Instead, I give each K_val the same weight
  k[l] <- rcat(1, pi[l, 1:K_val])
  y[l] <- rcat(1, Theta[k[l], 1:K_imp])
}



# Vary the fraction of individuals that are validated ---------------------

# For this simulation, we are changing the fraction of samples that were validated.
# At this point, the object k stores the validated state value for 100% of samples.
# We need to vary the the number of values in k that are NA, thereby varying the fraction
# of individuals that are validated 

k_fracs <- c(1.0, 0.95, 0.7, 0.5, 0.3, 0.15, 0.11, 0.05) 
#AIS do we want to change these values?

k_mat <- matrix(data=k, nrow=Ltot, ncol=length(k_fracs))
colnames(k_mat) <- paste0("k",k_fracs)
for (j in 2:ncol(k_mat)) {
  k_mat[,j] <- k_mat[,j-1]
  alreadyNA <- sum(is.na(k_mat[,j]))
  needtobeNA <- round((1-k_fracs[j])*Ltot)
  ind <- sample(which(!is.na(k_mat[,j])), (needtobeNA - alreadyNA), replace=F)
  k_mat[ind,j]<- NA
}

# Sanity check
# calculate the number of values in each column that are NA
for (j in 1:ncol(k_mat)){
  print(1-sum(is.na(k_mat[ ,j]))/Ltot)
}
# should match k_fracs (in same order)

# Specify which k_frac will be used to initialize z.init, z.dat, and run the model
# manually specify which column index of k_mat
for (current in 1:length(k_fracs)) {
  k_sim <- k_mat[,current]
  #AIS make this line ^ through running the model into a for-loop
  
  # Partially observed occupancy --------------------------------------------
  
  z.dat_df <- data.frame(site,k_sim) %>%
    filter(!is.na(k_sim)) #filter out the NAs so that df contains only validated individuals
  z.dat_df$occ <- 1
  # cast into a nsite x k_val matrix that is full of NAs and 1s
  z.dat_cast <- z.dat_df %>%
    reshape2::acast(site ~ k_sim, fill=-999, drop=F, value.var = "occ")
  z.dat_cast[z.dat_cast == -999] <- NA
  z.dat_cast[z.dat_cast > 0] <- 1
  z.dat <- z.dat_cast
  
  # Initialize Z
  z.init <- z
  z.init[z.dat == 1] <- NA # initialize known values as NA, otherwise model will throw error
           
  # Check that where L>0 for a species, z.init>0 for that species/site/year combo
  for (i in 1:dim(z.init)[1]) {
    for (j in 1:dim(z.init)[2]) {
        if (sum(L[i,], na.rm = TRUE) > 0 ) {
          ifelse(z.init[i,j] == 0, 1, z.init[i,j])
      }
    }
  }
  
  # AIS comapre z.init and z. Is creating a z.init object from z and z.dat the right thing to do?
  
  # Fit model ---------------------------------------------------------------
  jags_d <- list(nsite = dim(L)[1],
                 nsurv = dim(L)[2], 
                 K_val = dim(alpha)[1],
                 K_imp = dim(alpha)[2],
                 L = L, 
                 alpha = alpha,
                 Ltot = sum(L), 
                 site = site,
                 # if the individual was labeled by the expert, true ID is known
                 k = k_mat[ ,current],
                 # for all individuals, we get paratxonomist IDs
                 y = y,
                 z = z.dat,
                 R = R)
  JAGSinits <- function(){ list(z = z.init) }
  nc <- 4
  ni <- 5000
  cl <- makeCluster(nc)
  jm <- jags.parfit(cl = cl,
                    data = jags_d,
                    params = c("psi","lambda", "Theta","eps_site", "eps_spec")
                               #"Tau_spec", "Tau_site"),
                    model = "sim_full_JAGS.txt",
                    n.chains = nc,
                    n.adapt = 2000,
                    n.update = 2000,
                    thin = ni/1000,
                    n.iter = ni) 
      
  saveRDS(jm, paste0("output/sim_k",k_fracs[current],"_jm.rds"))
}

    
#AIS discuss: if we validate 11% of samples, the model runs, but not if we tak 10% or less

# Assess model results ----------------------------------------------------

temp = list.files(path = "output/", pattern = "sim_k*")
for (i in 1:length(temp)) assign(temp[i], readRDS(paste0("output/",temp[i])))

k1_summ <- MCMCsummary(sim_k1_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[1]) %>%
    rownames_to_column()
k0.95_summ <- MCMCsummary(sim_k0.95_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[2]) %>%
    rownames_to_column()
k0.7_summ <- MCMCsummary(sim_k0.7_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[3]) %>%
    rownames_to_column()
k0.5_summ <- MCMCsummary(sim_k0.5_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[4]) %>%
    rownames_to_column()
k0.3_summ <- MCMCsummary(sim_k0.3_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[5]) %>%
    rownames_to_column()
k0.15_summ <- MCMCsummary(sim_k0.15_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[6]) %>%
    rownames_to_column()
k0.11_summ <- MCMCsummary(sim_k0.11_jm.rds, params = 'psi', round=2) %>%
    mutate(fraction = k_fracs[7]) %>%
    rownames_to_column()
           
sims_jm <- rbind(k1_summ, k0.95_summ, k0.7_summ, k0.5_summ, k0.3_summ,
                 k0.15_summ, k0.11_summ) 


# Iterations
  # full vs reduced
  # by species
  # by various fractions of validation

# 1) Plot occupancy estimate with CI's along with true parameter value as a line for each species (see Max's old code). Do for full and reduced model. What I'm imagining is overlapping the posterior distribution for each fraction of validation data simulation. Then plot a straight line over the true occupancy/demographic rate 
psi_1site <- data.frame(true=psi[1,], species=as.factor(1:K_val)) 
  
sims_jm %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species"), sep = ",") %>%
  mutate_at(c('site', 'species'), readr::parse_number) %>%
  mutate(species = as.factor(species),
         fraction = as.factor(fraction)) %>%
  # AIS how to combine all of the sites together for each species? I could take their average (average 2.5%, 50%, 97.5%). For now, I  select just one site.
  filter(site==1) %>%
  ggplot(aes(fraction, `50%`, group=site)) + 
  geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="gray40", alpha=.25, size=1.5) +
  geom_point(size = 3, color="gray20") +
  geom_hline(data=psi_1site, aes(yintercept=true), col="red") +
  facet_wrap(~species) + 
  xlab("Fraction of validated individuals") + 
  ylab("Occupancy probability") + 
  theme_minimal()
  
# 2) Plot encounter rate estimate with CI's along with true parameter value as a line for each species (see same old code). Do for full and reduced model
    
# 3) Report on accuracy, recall, etc. for full and reduced model. Is there evidence that the reduced model is less effective for prioritizing confirmation effort than the proposed model. Address this by calculating accuracy, recall, etc. for full vs reduced model for each simulation. This will speak to reviewer feedback: "it was not clear that this reduced model was markedly worse for the purpose of prioritizing confirmation effort. Were patterns in the theta matrix markedly different for the reduced model (I see Fig 4, but more broadly)?"

# 4) make it more explicitly ecologically meaningful: Discuss whether the credible interval contains the true value. To what degree does partial supervision improve estimation of psi.lambda with varying levels of confirmation effort? Do the same for accuracy, recall, etc. Consider how sensitive the model output is to the fraction validated - Maybe the # IDed by the expert doesn’t affect the occupancy estimate when the occupancy is high, but does when occupancy is low….

# After sending questions to max
    # create R and JAGS scripts for reduced model
    # create add validation metrics script - see validation.R as model






# MISC
# AIS questions:
  # Should we allocate validation randomly across samples or put a weight on samples that are not a validated state? For now, I'll just do random because that's easier
  # Should the same individuals not validated when we take 0.5 of k also not be validated when we take 0.2 of k? That is, as we decrease the fraction of k that informs the model, should we keep the same individuals unvalidated, or randomly select again? To make comparing results across varying k fractions more fair, I think we should take subsets of already unvalidated individuals (e.g if indiviual x is not validated when we take 0.8 of k, then individual x should be not validated for all smaller fractions of k)
  #should we consider varying accuracy of imperfect classifier? probably not. That is set randomly in priors
  # Consider expanding to dynamic model to calculate colonization/extinction. I'm worried just a single season doesn't address reviewer concerns enough
  # This simulation doesn't speak to this: "What constitutes a usefully-sized expert-verified sample, how to allocate expert effort across species/sites/particular samples, the choice of priors for theta, specification of/sensitivity to the assumed count distribution". How do we respond to reviewer?





# SCRAP NOTES 

# See Brett's tutorial: https://github.com/EBIO5460Fall2018/Class-materials/blob/master/12_4_sim_multilevel.md
