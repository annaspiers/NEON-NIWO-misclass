


# Here we extend the dynamic occupancy model with misclassification in script
# 14_ to run on actual NEON carabid data from Niwot Ridge. Here, we run the data
# on all samples from Niwot Ridge

library(dclone) #alternative: R2jags::jags.paralle
library(MCMCpack)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyselect)
# library(mailR) #try troubleshooting dependencies later
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
# Each row is a beetle identified by a parataxonomist
ind_dat <- read.csv("occupancy/model_df_by_individual_beetle.csv", header = TRUE) %>%
    mutate(collectDate = as.Date(collectDate)) %>%
    filter(collectDate <= as.Date("2018-12-31")) %>%
    dplyr::select(-c(X))

# Each row is a sample. One sample is a para_morph at a trap on a collection day   
sample_dat  <- read.csv("occupancy/model_df_by_species_in_sample_2015-2018.csv", 
                        header = TRUE) %>%
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
    select(para_sciname) %>%
    filter(!para_sciname %in% unique(ind_dat$expert_sciname)) %>% #we want rows where the paratax ID is not an existing expert ID
    mutate_if(is.factor, as.character) %>%
    distinct()
# 14 parataxonomist IDs that the expert never ID'd

# Are there species that the expert taxonomist ID’d, but that the parataxonomist never ID’d?
ind_dat %>%
    select(expert_sciname) %>%
    filter(!is.na(expert_sciname),
           !expert_sciname %in% unique(ind_dat$para_sciname)) %>% #we want rows where the paratax ID is not an existing expert ID
    mutate_if(is.factor, as.character) %>%
    distinct()
# 7 expert tax IDs that the paratax never ID'd


# Misclassification model ------------------------------------------------- 
# The parameters below assume both a parataxonomist's and an expert taxonomist's
# classification is available

# List the unique expert taxonomist ID's and the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the rownames of theta and M
rownames <- %>% arrange(expert_sciname)

ind_dat %>%
    select(para_sciname) %>%
    filter(!para_sciname %in% unique(ind_dat$expert_sciname)) %>% 
    mutate_if(is.factor, as.character) %>%
    distinct()

ind_dat %>% 
    select(expert_sciname) %>% 
    filter(!is.na(expert_sciname)) %>%
    distinct() 
# List all unique unique taxonomist ID's total
allNames    <- unique(unique(ind_dat$para_morph, na.rm=T), expertnames)
, ) %>%
    distinct()

ind_dat %>%
    filter(!is.na(expert_sciname)) %>%
    select(expert_sciname) %>%
    distinct() %>% 
    mutate_if(is.factor, as.character) %>%
    arrange(expert_sciname) %>%
    pull(expert_sciname) 

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
# Part 1: a square matrix where each row/column represents a 'valid' taxonomic
# group 
# Part 2: matrix appended as new columns on Theta that are “invalid” taxonomic
# groups (all of the morphospecies)
nc <- 4 #MCMC chains
ni <- 4000 #MCMC iterations
nb <- 2000 #MCMC burnin = adapt (1000) + update (1000)
theta <- array(NA, dim=c(nc*(ni - nb), dim(alpha)) ) #same dimensions as theta output from model
for (k in 1:nrow(alpha)) {
    theta[ ,k, ] <- MCMCpack::rdirichlet(dim(theta)[1], alpha[k, ])
}

# M_k: M[k,k'] is the number of individuals from species k (according to expert)
# that were identified as species k' by parataxonomist [KxK]
M_sq <- ind_dat %>%
    filter(!is.na(expert_sciname)) %>%
    reshape2::acast(expert_sciname ~ para_morph) 

colnames(M_sq) %in% expertnames
which(colnames(M_sq) %in% expertnames)

test_sq<- M_sq[,which(colnames(M_sq) %in% expertnames)]

M_sq_tib <- tibble(M_sq)
as.tibble() %>%
    dplyr::select(expertnames)

M_sq <- ind_dat %>%
    filter(!is.na(expert_sciname) |
               para_morph %in% expertnames) %>%
    reshape2::acast(expert_sciname ~ para_morph)

#filter(is.na(expert_sciname)==F) %>% 
#filter(as.character(para_sciname) == as.character(expert_sciname)) %>%
#mutate_if(is.factor, as.character)  %>%
reshape2::acast(expert_sciname ~ para_sciname)

dplyr::select(all_of(expertnames))

vars_select(colnames(M_sq), .include = expertnames)

M <- ind_dat %>%
    filter(is.na(expert_sciname)==F) %>% 
    reshape2::acast(expert_sciname ~ para_morph)



c(ind_dat %>%
      filter(!is.na(expert_sciname)) %>%
      select(expert_sciname) %>%
      distinct() %>%
      mutate_if(is.factor, as.character) %>% pull(expert_sciname) ) )


select(ind_dat %>%
           filter(!is.na(expert_sciname)) %>%
           select(expert_sciname) %>%
           distinct() %>%
           mutate_if(is.factor, as.character) %>% pull(expert_sciname) )

M <- ind_dat %>%
    filter(is.na(expert_sciname)==F) %>% 
    reshape2::acast(expert_sciname ~ para_morph)
# Reorder M so that first para_morph names match expert_scinames

#try creating M as a square matrix, then left_join full rectangular matrix

# Combine occupancy and misclassification models to simulate observed data --------
# c_obs: c_obs[i,j,k',l] are the elements of vector C, and represent the number
# of individuals that were classified as k'. dim: nsite x nsurv x nspec x nyear
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
cl <- makeCluster(4)
out <- jags.parfit(cl = cl, 
                   data = JAGSdata,
                   params = JAGSparams,
                   model = "occupancy/15_neon_dynamic_multisp_misclass_JAGS.txt", 
                   inits = JAGSinits,
                   n.chains = nc,
                   n.adapt = 2000,
                   n.update = 2000,
                   thin = nt,
                   n.iter = 4000)
# ran overnight, less than 11 hours
# saveRDS(out, "occupancy/script15_jags_out.rds")

# How well does the model estimate specified parameters?
# out_mcmcsumm <- MCMCsummary(out)
# saveRDS(out_mcmcsumm, "occupancy/script15_mcmcsumm.rds")


# View JAGS output --------------------------------------------------------

out <- readRDS("occupancy/script15_jags_out.rds")
out_mcmcsumm <- readRDS("occupancy/script15_mcmcsumm.rds")

# Did model converge?
hist(out_mcmcsumm$Rhat)
range(out_mcmcsumm$Rhat, na.rm=TRUE)

# mcmc.list visualization resource: https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Look at raw numbers
# phi - survival probability
MCMCsummary(out, params = 'phi', round=2)
MCMCtrace(out, params = 'phi', type = 'density', ind = F, pdf=F)
#phi[35] Rhat=1.35. 
#phi indices 14,16,21,25,32 are bimodal

# gamma - colonization probability
MCMCsummary(out, params = 'gamma', round=2)
MCMCtrace(out, params = 'gamma', type = 'density', ind = F, pdf=F)
#gamma[35] Rhat=1.26. same index as above for phi. 
#all densities peak close to 0 - makes sense

# psi - occupancy prob.
MCMCsummary(out, params = 'psi', round=2)
MCMCtrace(out, params = 'psi', type = 'density', ind = F, pdf=F)
# a few Rhat values > 1.1
#some have wide posterior densities.

# Plot all species' occupancy through seasons
par(mfrow=c(1,1))
plot(NA, xlim = c(1,dim(c_obs)[4]), ylim=c(0,1), main="Occupancy by species", 
     xlab = "Year", ylab = "Occupancy probability", frame.plot = FALSE)
for (k in 1:dim(Z.dat)[2]) {
    lines(1:dim(Z.dat)[3], MCMCsummary(out, params=paste0('psi\\[',k,',\\d\\]'), ISB=F)$mean,
          type = "l", col=k+7,  lwd = 2, lty = 1, las = 1)
}

# lambda - expected abundance, given occupancy, dim: nsite x nsurv x nspec x nyear
#lambda_summ <- MCMCsummary(out, params = 'lambda', round=2) #takes 5ish min
#saveRDS(lambda_summ, "occupancy/lambda_summ.rds")
lambda_summ <- readRDS("occupancy/lambda_summ.rds")
MCMCtrace(out, params = 'lambda', type = 'density', ind = F, pdf=F)
range(lambda_summ$mean)

# theta
#theta_summ <- MCMCsummary(out, params = 'theta', round=2)
#saveRDS(theta_summ, "occupancy/theta_summ.rds")
theta_summ <- readRDS("occupancy/theta_summ.rds")
MCMCtrace(out, params = 'theta', type = 'density', ind = F, pdf=F)

theta_df <- data.frame(expert_index = rep(1:dim(theta)[2], dim(theta)[3]) ,
                       paramorph_index = rep(1:dim(theta)[3], each = dim(theta)[2]))  %>% 
    left_join(ind_dat %>%
                  filter(!is.na(expert_sciname)) %>%
                  select(expert_sciname) %>%
                  distinct() %>%
                  mutate_if(is.factor, as.character) %>%
                  arrange(expert_sciname) %>%
                  mutate(expert_index = 1:dim(theta)[2]) )%>%
    left_join(ind_dat %>%
                  select(para_morph) %>%
                  distinct() %>%
                  mutate_if(is.factor, as.character) %>%
                  arrange(para_morph) %>%
                  mutate(paramorph_index = 1:dim(theta)[3]) ) %>%
    mutate(theta_mean = theta_summ$mean)

# Plot heatmap of theta values
ggplot(theta_df, aes(para_morph, expert_sciname, fill= theta_mean)) + 
    geom_tile() +
    scale_fill_gradient(low="darkblue", high="white") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))

#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html

# Z
out$mean$Z
plot(Z.init, out$mean$Z)
plot(apply(c_obs, c(1,3,4), max, na.rm = TRUE), out$mean$Z)
# not sure what to make of this

# n.occ
print(out$summary[grep("n.occ", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2)
# Yep, make sense. Low numbers for Amara quenseli in first year and for
# Pterostichus restrictus in first 2 years.

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
                                  sample_dat %>% filter(para_sciname==k_spec, sp_abund > 0) %>% 
                                      group_by(col_year) %>% 
                                      summarize(n_plots=n_distinct(plotID)))) + 
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


# Scratch -----------------------------------------------------------------

# Core and chain time trial -----------------------------------------------

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
# JAGSparams <- c("psi", "lambda", "theta", "Z", "phi", "gamma", "n.occ",
# "growth", "turnover")
# #nc, ni, nb defined above with theta
# nt <- 1  #MCMC thin
# #######################
# #4 chains: 38/37sec on 4 cores, 50sec on 3 cores, 38sec/38sec on 2 cores, 72sec on 1 core
# #3 chains on 2 and 4 cores
# #3 chains: on 4 cores is 28sec, and on 2 cores is 37/42sec

# For a while, I used function arguments for jags() in jags.parfit() rather than
# the proper arguments for jags.parfit(). Here are notes that describe the
# errors I ran into
#after 22 hours of running, this error was spat out: 
#Error in unserialize(node$con) : error reading from connection
#stackexchange says its due to running out of memory
#bumping down to 2 cores now. Kicking off Monday 21:10
# with 2 burnin and 10 iterations, took less than 2 hours - but this was using
# jgas() arguments, not jags.parfit


# Create matrix of mean theta posterior
# theta_mat <- matrix(data=NA, nrow=dim(theta)[2], ncol=dim(theta)[3], byrow=T,
#                     dimnames = list(c(ind_dat %>%
#                                           filter(!is.na(expert_sciname)) %>%
#                                           select(expert_sciname) %>%
#                                           distinct() %>%
#                                           mutate_if(is.factor, as.character) %>%
#                                           arrange(expert_sciname) %>%
#                                           pull(expert_sciname)),
#                                     c(ind_dat %>%
#                                           select(para_morph) %>%
#                                           distinct() %>%
#                                           mutate_if(is.factor, as.character) %>%
#                                           arrange(para_morph) %>%
#                                           pull(para_morph)) ) )
# for (k in 1:dim(theta)[1]) {
#     theta_mat[k,] <- MCMCsummary(out, params=paste0('theta\\[',k,',[0-9]+\\]'), 
#                                  ISB=F)$mean
# }
# #why is it missing column 37?
# theta_mat
# heatmap(theta_mat, Colv = NA, Rowv = NA, scale = "none",
#         col = heat.colors(100),
#         xlab="expert taxonomist", ylab="parataxonomist", 
#         main="Theta Posterior Means")
