# Here we extend the dynamic occupancy model with misclassification in script 13_ to run on actual NEON carabid data from Niwot Ridge. Here, we run the data on all samples from Niwot Ridge

library(jagsUI)
library(dclone)
library(MCMCpack)
library(MCMCvis)
library(ggmcmc)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggplot2)
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
# Each row is a beetle identified by a parataxonomist
ind_dat <- read.csv("occupancy/model_df_by_individual_beetle.csv", header = TRUE) %>%
    mutate(collectDate = as.Date(collectDate)) %>%
    filter(collectDate <= as.Date("2018-12-31")) %>%
    dplyr::select(-c(X))

# Each row is a sample. One sample is a para_morph at a trap on a collection day   
sample_dat  <- read.csv("occupancy/model_df_by_species_in_sample_2015-2018.csv", header = TRUE) %>%
    dplyr::select(-c(X))


# EDA ---------------------------------------------------------------------

# How individuals for each expert species ID?
ind_dat %>%
    group_by(expert_sciname) %>%
    summarize(n=n()) %>% 
    arrange(-n) %>% data.frame()
# 39 unique IDs, including NAs

# How individuals for each parataxonomist species ID?
ind_dat %>%
    group_by(para_sciname) %>%
    summarize(n=n()) %>% 
    arrange(-n) %>% data.frame()
# 45 unique IDs, no NAs

# How individuals for each morphospecies ID?
ind_dat %>%
    group_by(morphospeciesID) %>%
    summarize(n=n()) %>% 
    arrange(-n) %>% data.frame()
# 24 unique IDs, including NAs

# Look at para_scinames associated with morphospecies IDs.
ind_dat %>%
    filter(!is.na(morphospeciesID)) %>%
    group_by(para_sciname, morphospeciesID) %>%
    summarize(n=n()) %>% 
    arrange(-n) %>% data.frame()
# Morphospecies exist when parataxonomist could not ID to species level
# Morphs have unique IDs within a year. What if the same morph is found between 2015 and 2016 - do these have different morph IDs? probably

# How many unique combinations of para_sciname and morphospecies (when a morph exists)?
ind_dat %>%
    group_by(para_morph) %>%
    summarize(n=n()) %>% 
    arrange(-n) %>% data.frame()
# 63 unique IDs

# Are there species that the parataxonomist ID’d, but that the expert never ID’d?
ind_dat %>%
    filter(is.na(expert_sciname),   #we want rows where there is no expert ID
           !para_sciname %in% unique(ind_dat$expert_sciname))  #we want rows where the paratax ID is not an existing expert ID
# No parataxonomist IDs that the expert never ID'd

# Misclassification model -------------------------------------------------
# The parameters below assume both a parataxonomist's and an expert taxonomist's classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n <- ind_dat %>%
    filter(!is.na(expert_sciname)) %>%
    group_by(expert_sciname) %>%
    summarize(n=n()) %>% 
    pull(n) #alphabetical


# Alpha: dirichlet concentration parameter
alpha <- matrix(1, nrow=length(n), 
                ncol=(length(unique(ind_dat$para_morph) ) ) )
diag(alpha) <- 20 

# Theta: misclassification probability matrix [KxK]
# Theta has two parts. 
# Part 1: a square matrix where each row/column represents a 'valid' taxonomic group
# Part 2: matrix appended as new columns on Theta that are “invalid” taxonomic groups (all of the morphospecies)
nc <- 4 #MCMC chains
ni <- 8000 #MCMC iterations
nb <- 2000 #MCMC burnin
theta <- array(NA, dim=c(nc*(ni - nb), dim(alpha)) ) #same dimensions as theta output from model
for (k in 1:nrow(alpha)) {
    theta[ ,k, ] <- MCMCpack::rdirichlet(dim(theta)[1], alpha[k, ])
}


# M_k: M[k,k'] is the number of individuals from species k (according to expert) that were identified as species k' by parataxonomist [KxK]
M <- ind_dat %>%
    filter(is.na(expert_sciname)==F) %>% 
    reshape2::acast(expert_sciname ~ para_morph)


# Combine occupancy and misclassification models to simulate observed data --------
# c_obs: c_obs[i,j,k',l] are the elements of vector C, and represent the number of individuals that were classified as k'. dim: nsite x nsurv x nspec x nyear
c_obs   <- sample_dat %>%
    reshape2::acast(plotID ~ col_index ~ para_morph ~ col_year,
                    fun.aggregate=sum, fill=-999, value.var = "sp_abund")
c_obs[c_obs == -999] <- NA

# Occupancy array. dim: nsite x nspec x nyear
# Create Z data. Use expert identifications to incorporate partly observed presence
Z.dat <- ind_dat %>% #ind_dat df has expert IDs
    filter(!is.na(expert_sciname)) %>%
    mutate(occ = 1) %>%
    reshape2::acast(plotID ~ expert_sciname ~ col_year,
                    fill=-999, drop=F, value.var = "occ")
Z.dat[Z.dat == -999] <- NA
Z.dat[Z.dat > 0] <- 1

# Initialize Z
Z.init <- Z.dat
for (i in 1:dim(Z.init)[1]) {
    for (l in 1:dim(Z.init)[3]) {
        Z.init[i,,l] <- sample(c(0,1), replace=TRUE, size=dim(Z.init)[2])
    }
}
# initialize known values as NA, otherwise model will throw error
Z.init[Z.dat == 1] <- NA

# Check that where c_obs>0 for a species, Z.init>0 for that species/site/year combo
for (i in 1:dim(Z.init)[1]) {
    for (k in 1:dim(Z.init)[2]) {
        for (l in 1:dim(Z.init)[3]) {
            if (sum(c_obs[i,,k,l], na.rm = TRUE) > 0 ) {
                ifelse(Z.init[i,k,l] == 0, 1, Z.init[i,k,l])
            }
        }
    }
}


# JAGS model --------------------------------------------------------------

# Run model in JAGS. 
JAGSdata <- list(nsite = dim(c_obs)[1], 
                 nsurv = dim(c_obs)[2], 
                 nspec_exp = dim(alpha)[1], #species ID'd by expert
                 nspec_para = dim(alpha)[2], #paratax ID's and morphospecies
                 nyear = dim(c_obs)[4], 
                 n = n,
                 alpha = alpha,
                 M = M,
                 c_obs = c_obs,
                 Z = Z.dat) #bundle data
JAGSinits <- function(){ list(Z = Z.init) }
JAGSparams <- c("psi", "lambda", "theta", "Z", "phi", "gamma", "n.occ", "growth", "turnover") 
#nc, ni, nb defined above with theta
nt <- 1  #MCMC thin



# JAGS model
out <- jags.parfit(cl = makeCluster(2),
                data = JAGSdata,
                inits = JAGSinits,
                params = JAGSparams,
                model = "occupancy/15_neon_dynamic_multisp_misclass_JAGS.txt", 
                n.chains = nc,
                n.iter = ni,
                n.burnin = nb,
                n.thin = nt)
stopCluster(cl)
#started Monday at 15:11 then crashed Tuesday night. Restarted Wednesday 8:53
# AIS try parallelizing with R2jags::jags.parallel or dclone

#save(out,file="occupancy/script15_jags_output.Rdata")
load("occupancy/script15_jags_output.Rdata") #output from jags.parfit()
out_15 <- out
load("occupancy/script14_jags_output.Rdata") #output from jags()

# How well does the model estimate specified parameters?
out_mcmcsumm <- MCMCsummary(out_15)
save(out_15, out_mcmcsumm, file="occupancy/script15_jags_output.Rdata")

# Look at raw numbers
# phi - survival probability
print(out$summary[grep("phi", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
MCMCtrace(out, params = 'phi', type = 'density', ind = F, pdf=F)

# gamma - colonization probability
print(out$summary[grep("gamma", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2)
# Good, these should be low. Higher for Amara quenseli(2) and Pterrostichus restrictus(7), the species that popped up after the first year of sampling

# psi - occupancy prob.
print(out$summary[grep("psi", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
# Pretty low occupancy. This makes sense with the sparse distributions at Niwot. reasonably small sd 
# Plot all species' occupancy through seasons
plot(NA, xlim = c(1,dim(c_obs)[4]), ylim=c(0,1), main="Occupancy by species", xlab = "Year", ylab = "Occupancy probability", frame.plot = FALSE)
for (k in 1:dim(c_obs)[3]) {
    lines(1:dim(c_obs)[4], out$mean$psi[k,], type = "l", col=k+7,  lwd = 2, lty = 1, las = 1)
}
# Why do they seem to converge?
# Plot posterior vs prior for each species in each year
psi_df <- data.frame(value = numeric(), 
                            post_prior = character(),
                            index = character())
for (i in 1:dim(out$sims.list$psi)[2]) {
    for (j in 1:dim(out$sims.list$psi)[3]) {
        post_df <- data.frame(out$sims.list$psi[ ,i,j]) %>%
            rename(value = out.sims.list.psi...i..j.) %>%
            mutate(post_prior = "posterior",
                   index = paste0("psi[",i,",",j,"]"))
        prior_df <- data.frame(psi[ ,i,j]) %>%
            rename(value = psi...i..j.) %>%
            mutate(post_prior = "prior",
                   index = paste0("psi[",i,",",j,"]"))
        psi_df <- rbind(psi_df, post_df, prior_df)
    }
}

# Plot matrix of theta prior and poserior densities
ggplot(psi_df, aes(x=value,y=..scaled..)) +
    geom_density(aes(color=post_prior)) + 
    facet_wrap( ~ index, scales="free_x") +
    xlab("Theta value") + scale_y_continuous(breaks=seq(0, 1, 0.5))


# lambda - expected abundance, given occupancy, dim: nsite x nsurv x nspec x nyear
print(out$summary[grep("lambda", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
range(out$mean$lambda)

# theta
theta_df <- data.frame(value = numeric(), 
                            post_prior = character(),
                            index = character())
for (i in 1:dim(out$sims.list$theta)[2]) {
    for (j in 1:dim(out$sims.list$theta)[3]) {
        post_df <- data.frame(out$sims.list$theta[ ,i,j]) %>%
            rename(value = out.sims.list.theta...i..j.) %>%
            mutate(post_prior = "posterior",
                   index = paste0("theta[",i,",",j,"]"))
        prior_df <- data.frame(theta[ ,i,j]) %>%
            rename(value = theta...i..j.) %>%
            mutate(post_prior = "prior",
                   index = paste0("theta[",i,",",j,"]"))
        theta_df <- rbind(theta_df, post_df, prior_df)
    }
}

# Plot matrix of theta prior and poserior densities
ggplot(theta_df, aes(x=value,y=..scaled..)) +
    geom_density(aes(color=post_prior)) + 
    facet_wrap( ~ index, scales="free_x") +
    xlab("Theta value") + scale_y_continuous(breaks=seq(0, 1, 0.5))

# Z
out$mean$Z
plot(Z.init, out$mean$Z)
plot(apply(c_obs, c(1,3,4), max, na.rm = TRUE), out$mean$Z)
# not sure what to make of this

# n.occ
print(out$summary[grep("n.occ", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2)
# Yep, make sense. Low numbers for Amara quenseli in first year and for Pterostichus restrictus in first 2 years.

# growth
print(out$summary[grep("growth", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
out$mean$growth
# Large growth values for species 2 and 7 - good

# turnover
print(out$summary[grep("turnover", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
out$mean$turnover
# numbers don't look crazy, but also not sure what to make of them

# Look at psi, phi, gamma, growth, turnover, n.occ graphically
species <- as.character(unique(sample_dat$para_sciname))
for (k in 1:length(species)) {
    k_spec <- species[k]
        
    p1 <- ggplot(data = sample_dat %>% filter(para_sciname == k_spec)) + 
        geom_col(aes(x=collectDate,y=sp_abund)) +
        ggtitle(k_spec) + xlab("Collection Date") + ylab("Abundance")
    p2 <- ggplot(data = left_join(sample_dat %>% dplyr::select(col_year) %>% distinct(),
                                  sample_dat %>% filter(para_sciname==k_spec, sp_abund > 0) %>% group_by(col_year) %>% summarize(n_plots=n_distinct(plotID)))) + 
        geom_line(aes(x=col_year, y=n_plots),) +
        geom_line(aes(x=col_year,y=out$mean$n.occ[k,]),col="red") +
        annotate("text", x=2017, y=50, label = "Predicted occupied (n.occ)", col="red") +
        annotate("text", x=2017, y=45, label = "Observed") +
        ggtitle("Site occupancy/detection") + xlab("Year") + ylab("Number of sites")
    p3 <- ggplot() + geom_line(aes(x=unique(sample_dat$col_year),
                                   y=out$mean$psi[k,]), col="red") +
        ggtitle("Occupancy (psi)") + xlab("Year") + ylab("Probability") + ylim(c(0,1))
    p4 <- ggplot() + 
        geom_line(aes(x=unique(sample_dat$col_year),y=out$mean$growth[k,]),col="blue") + 
        annotate("text", x=2017, y=0.6, label = "Growth", col="blue") +
        geom_line(aes(x=unique(sample_dat$col_year),y=c(out$mean$turnover[k,],NA)),col="darkgreen") + 
        annotate("text", x=2017, y=0.5, label = "Turnover", col="darkgreen") +
        geom_hline(yintercept=out$mean$phi[k],col="purple") + 
        annotate("text", x=2017, y=0.4, label = "Survival (phi)", col="purple") +
        geom_hline(yintercept=out$mean$gamma[k],col="orange") + 
        annotate("text", x=2017, y=0.3, label = "Colonization (gamma)", col="orange") +
        ggtitle("Demographic rates") + xlab("Year") + ylab("Rate") + ylim(c(0,1))
        
    grid.arrange(p1,p2,p3,p4,nrow=2)
}
par(mfrow=c(1,1)) #reset plotting

# Visualize predictions of species unobserved by expert taxonomist
#par(mfrow=c(1,1))
#plot(NA,xlim=c(0,1),ylim=c(0,1),
#     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
#abline(0,1)
#points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")


# Core and chain time trial -----------------------------------------------
# 
# ###############
# alpha <- matrix(1, nrow=5,ncol=5)
# diag(alpha) <- 5 
# M <- matrix(0, nrow=5,ncol=5)
# diag(M) <- 7
# c_obs <- array(1:8,c(5,3,5,2))
# Z.init <- array(0:1,c(5,5,2))
# JAGSdata <- list(nsite = 5, 
#                  nsurv = 3, 
#                  nspec_exp = 5, #species ID'd by expert
#                  nspec_para = 5, #paratax ID's and morphospecies
#                  nyear = 2, 
#                  n = rep(7,5),
#                  alpha = alpha,
#                  M = M,
#                  c_obs = c_obs) #bundle data
# JAGSinits <- function(){ list(Z = Z.init) }
# JAGSparams <- c("psi", "lambda", "theta", "Z", "phi", "gamma", "n.occ", "growth", "turnover") 
# #nc, ni, nb defined above with theta
# nt <- 1  #MCMC thin
# #######################
# #4 chains: 38/37sec on 4 cores, 50sec on 3 cores, 38sec/38sec on 2 cores, 72sec on 1 core
# #3 chains on 2 and 4 cores
# #3 chains: on 4 cores is 28sec, and on 2 cores is 37/42sec
