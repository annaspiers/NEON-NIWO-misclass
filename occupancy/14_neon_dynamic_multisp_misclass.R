# Here we extend the dynamic occupancy model with misclassification in script 13_ to run on actual NEON carabid data from Niwot Ridge. Here, we first we run the model with the 7 most abundant species. These 7 species are perfectly identified by the parataxonomists. Not every individual beetle has an expert identification. Later we will try the model on all species entirely

library(MCMCpack) #rdirchlet
library(jagsUI)
library(dplyr)
library(reshape2)

# Load NEON Niwot Ridge carabid data 
sample_dat  <- read.csv("data_derived/model_df_by_species_in_sample.csv", header = TRUE)  # each row is a sample. All samples are represented. One sample is a species at a trap on a collection day
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




# AIS - in the morning make graphs to compare modelled data to raw data




# Do the modeled values make sense?
# psi - occupancy prob.
out$mean$psi
# Pretty low occupancy - none exceed 0.5. This makes sense with the sparse distributions at Niwot 

# lambda
#out$mean$lambda

# theta
out$mean$theta
# Makes sense. There were no misidentifications by the parataxonomist. I wonder why some species are in the 80s% though

# Z
#out$mean$Z

# p - detection prob.
out$mean$p
# Why are these all around 50% there must be some statistical reason

# phi - survival prob.
out$mean$phi
# Seem reasonable at first glance

# gamma - colonization prob.
out$mean$gamma
# Good, these should be low. Higher for Amara quenseli(2) and Pterrostichus restrictus(7), the species that popped up after the first year of sampling

# n.occ
out$mean$n.occ
# Yep, make sense. Low numbers for Amara quenseli in first year and for Pterostichus restrictus in first 2 years.

# growth
# turnover



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

# Compare true and recaptured Z 
plot(c(Z), c(out$mean$Z), main="Z")
abline(0,1)

# Visualize predictions of species unobserved by expert taxonomist
#par(mfrow=c(1,1))
#plot(NA,xlim=c(0,1),ylim=c(0,1),
#     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
#abline(0,1)
#points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")

# Plot apparent occupancy 
    psi.app <- apply(apply(y, c(1,3), max), 2, mean) 
    lines(year, psi.app, type = "l", col = "black", lwd = 2) 
    text(0.85*K, 0.06, labels = "red solid – true occupancy\n red dashed – detection probability\n black – observed occupancy")