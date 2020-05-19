# Here we extend the dynamic occupancy model with misclassification in script 13_ to run on actual NEON carabid data from Niwot Ridge. Here, we first we run the model with the 7 most abundant species. These 7 species are perfectly identified by the parataxonomists. Not every individual beetle has an expert identification. Later we will try the model on all species entirely

#library(MCMCpack) #rdirchlet
library(jagsUI)
library(dplyr)
library(reshape2)
library(gridExtra) #gridarrange

# Load NEON Niwot Ridge carabid data 
sample_dat  <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE) %>%
    mutate(collectDate = as.Date(collectDate))# each row is a sample. All samples are represented. One sample is a species at a trap on a collection day
ind_dat     <- read.csv("data_derived/model_df_by_individual_beetle.csv", header = TRUE) %>%
    mutate(expert_sciname = as.factor(ifelse(is.na(expert_sciname) == T, "zno_exp_ID", as.character(expert_sciname)))) # each row is a beetle identified by a parataxonomist


# Misclassification model -------------------------------------------------
# The parameters below assume both a parataxonomist's and an expert taxonomist's classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n <- ind_dat %>%
    filter(expert_sciname != "zno_exp_ID") %>%
    group_by(expert_sciname) %>%
    summarize(n=n()) %>% 
    pull(n)

# Alpha: dirichlet concentration parameter
alpha <- matrix(1, nrow=length(n), 
                ncol=length(unique(ind_dat$para_sciname)) ) #AIS every individual has a para_sciname ID, so there are no morphospecies. However, when we extend the model to all species, then we will need to add the number of unique morphospecies to the alpha columns
diag(alpha) <- 15 #place higher probability mass on parataxonomist getting the correct classification

# M_k: M[k,k'] is the number of individuals from species k (according to expert) that were identified as species k' by parataxonomist [KxK]
M <- ind_dat %>%
    reshape2::acast(expert_sciname ~ para_sciname)
# AIS there are no morphospecies for these 7 species. How unlikely. That means theta and M stay square


# Combine occupancy and misclassification models to simulate observed data --------

# c_obs: c_obs[i,j,k',l] are the elements of vector C, and represent the number of individuals that were classified as k'
c_obs   <- sample_dat %>%
    reshape2::acast(plot_trap ~ col_index ~ para_sciname ~ col_year, 
                    value.var = "occ")
str(c_obs) # nsite x nsurv x nspec x nyear

# Initialize occupancy array. dim: nsite x nspec x nyear
Z <-  apply(c_obs, c(1,3,4), max, na.rm = TRUE) 
Z[Z == "-Inf"] <- NA #this happens where a plot isn't sampled across years 

# JAGS model --------------------------------------------------------------

# Define function to run JAGS model
jags_misclass_fn <- function(){
    str(JAGSdata <- list(nsite = dim(c_obs)[1], #AIS although there are 44 unique traps, there are only 40 traps surveyed per year
                         nsurv = dim(c_obs)[2], #AIS is it alright that survey effort is unequal across years?
                         nspec = dim(c_obs)[3],
                         nyear = dim(c_obs)[4], 
                         n = n,
                         alpha = alpha,
                         M = M,
                         c_obs = c_obs)) #bundle data
    JAGSinits <- function(){list(Z = Z) }
    JAGSparams <- c("psi", "lambda", "theta", "Z", "p", "phi", "gamma", "n.occ", "growth", "turnover") #params monitored
    nc <- 4 #MCMC chains
    ni <- 4000 #MCMC iterations
    nb <- 1000 #MCMC burnin
    nt <- 1  #MCMC thin
    
    # JAGS model
    jags_out <- jags(data = JAGSdata,
                     inits = JAGSinits,
                     parameters.to.save = JAGSparams,
                     model.file = "occupancy/14_neon_dynamic_multisp_misclass_JAGS.txt", 
                     n.chains = nc,
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt)
    return(jags_out)
}

# Run model in JAGS. 
out <- jags_misclass_fn()
#save(out,file="occupancy/script14_jags_output.Rdata")

# How well does the model estimate specified parameters?
print(out, dig=2)
traceplot(out, parameters="theta")

# Did the model converge?
# Yes. With 4000 iterations, incuding 1000 for burnin, the model converged


# Do the modeled values make sense?

# Look at raw numbers
# phi - survival probability
print(out$summary[grep("phi", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
# Seem reasonable at first glance

# gamma - colonization probability
print(out$summary[grep("gamma", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
# Good, these should be low. Higher for Amara quenseli(2) and Pterrostichus restrictus(7), the species that popped up after the first year of sampling

# p - detection probability
print(out$summary[grep("p", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
# Funny that these are all around 50% - as good as a coin toss

# psi - occupancy prob.
print(out$summary[grep("psi", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) # nspec x nyear
# Pretty low occupancy - none exceed 0.5. This makes sense with the sparse distributions at Niwot. reasonably small sd 
# Plot all species' occupancy through seasons
plot(NA, xlim = c(1,dim(c_obs)[4]), ylim=c(0,1), main="Occupancy by species", xlab = "Year", ylab = "Occupancy probability", frame.plot = FALSE)
for (k in 1:dim(c_obs)[3]) {
    lines(1:dim(c_obs)[4], out$mean$psi[k,], type = "l", col=k+7,  lwd = 2, lty = 1, las = 1)
}

# lambda - expected abundance, given occupancy, dim: nsite x nsurv x nspec x nyear
print(out$summary[grep("lambda", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
range(out$mean$lambda)
# AIS hmm, this doesn't look right. Why are all of the values 0.4?

# theta
print(out$summary[grep("theta", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
out$mean$theta
# Makes sense. There were no misidentifications by the parataxonomist. I wonder why some species are in the 80s% though

# Z
out$mean$Z
plot(apply(c_obs, c(1,3,4), max, na.rm = TRUE), out$mean$Z)
# not sure what to make of this. Not bad, not great

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
    p2 <- ggplot(data = left_join(sample_dat %>% select(col_year) %>% distinct(),
                                  sample_dat %>% filter(para_sciname==k_spec, sp_abund > 0) %>% group_by(col_year) %>% summarize(n=n()))) + 
        geom_line(aes(x=col_year, y=n),) +
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
