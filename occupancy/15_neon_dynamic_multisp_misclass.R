# Here we create a dynamic occupancy model with misclassification using
# actual NEON carabid data from Niwot Ridge. Here, we run the data
# on all samples from Niwot Ridge 2015-2018
# mcmc.list visualization resource: https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

library(dclone) #alternative: R2jags::jags.parallel
library(MCMCpack)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyselect)
library(tibble)
library(assertthat)
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
all_paratax_df <- readRDS("occupancy/all_paratax_df.rds") #each row is a subsample (scientificname per trap per plot per collectdate)
pinned_df <- readRDS("occupancy/pinned_df.rds") #each row is an individual beetle that was pinned
expert_df <- readRDS("occupancy/expert_df.rds") #each row is an individual beetle identified by an expert

# Filter datasets to abundant parataxonomist species/morphospecies IDs
common_morpho <- all_paratax_df %>%
  group_by(scimorph_combo, individualCount) %>%
  summarize(total = sum(individualCount)) %>%
  arrange(-total) %>%
  filter(total > 10)

all_paratax_df <- all_paratax_df %>%
  filter(scimorph_combo %in% common_morpho$scimorph_combo)
pinned_df <- pinned_df %>%
  filter(scimorph_combo %in% common_morpho$scimorph_combo)
expert_df <- expert_df %>%
  left_join(pinned_df %>% select(individualID, scimorph_combo)) %>%
  filter(scimorph_combo %in% common_morpho$scimorph_combo)

# List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the row names of theta and M
rownames <- c(expert_df %>%   #1) the unique expert taxonomist ID's
                distinct(scientificName) %>%
                pull(scientificName),
              all_paratax_df %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                distinct(scientificName) %>%
                filter(!scientificName %in% unique(expert_df$scientificName)) %>% 
                pull(scientificName) ) 
rownames <- sort(rownames)
#AIS rather than using the parataxonomist scientific IDs that the expert taxonomist didn't use, 
# wouldn't we want to instead use paratazonomist scientificname-morphospecies combo that the expert didn't use?

# List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
# parataxonomist didn't use. 
# These identifications will be the column names of theta and M
extra_colnames <- c(all_paratax_df %>%   #1) the unique parataxonomist and morphospecies ID's
                      distinct(scimorph_combo) %>%
                      pull(scimorph_combo),
                    expert_df %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                      distinct(scientificName) %>%
                      filter(!scientificName %in% unique(all_paratax_df$scimorph_combo)) %>% 
                      pull(scientificName) ) 
extra_colnames <- sort(extra_colnames)
extra_indices <- which(!extra_colnames %in% rownames)
colnames <- c(rownames, extra_colnames[extra_indices])
rm(extra_colnames, extra_indices)

assertthat::assert_that(all(rownames == colnames[1:length(rownames)]))

# Misclassification model ------------------------------------------------- 
# The parameters below assume both a parataxonomist's and an expert taxonomist's
# classification is available

# n[k]: expert taxonomist's count of individuals of species k 
n_tib <- left_join(tibble(rownames) %>%
                     rename("scientificName" = rownames),
                   expert_df %>%
                     count(scientificName) )
n <- n_tib %>% pull(n) 
n[is.na(n)] <- 0

assertthat::assert_that(all(n_tib$scientificName == rownames))
rm(n_tib)

# Alpha: dirichlet concentration parameter
alpha <- matrix(3, nrow=length(rownames),
                ncol=length(colnames) ) 
diag(alpha) <- 180

test_alpha <- MCMCpack::rdirichlet(1000, c(200,rep(1.5,15)))
hist(test_alpha[,-1], xlim=c(0,1),col="red")
hist(test_alpha[,1],  add=T)

# Theta: misclassification probability matrix [KxK]
# Theta has two parts. 
# Part 1: a square matrix where each row/column represents a 'valid' taxonomic
# group 
# Part 2: matrix appended as new columns on Theta that are “invalid” taxonomic
# groups (all of the morphospecies)
# theta <- array(NA, dim=c(nc*(ni - nb), length(rownames), length(colnames)) ) #same dimensions as theta output from model
# for (k in 1:nrow(alpha)) {
#     theta[ ,k, ] <- MCMCpack::rdirichlet(dim(theta)[1], alpha[k, ])
# }


# M_k: M[k,k'] is the number of individuals from species k (according to expert)
# that were identified as species k' by parataxonomist [KxK]

# create an empty df with the columns and rownames that I want
M <- matrix(0, nrow = length(rownames), ncol = length(colnames))
M <- data.frame(M, row.names = rownames)
colnames(M) <- colnames

# Grab values from casted M matrix and fill in values in final M matrix
M_cast <- left_join(expert_df %>% select(individualID, scientificName_exp=scientificName),
                    pinned_df %>% select(individualID, scimorph_combo)) %>% 
  reshape2::acast(scientificName_exp ~ scimorph_combo)

for(col in colnames(M_cast)) {
  for(row in row.names(M_cast)) {
    M[row,col] <- M_cast[row,col]
  }
}
rm(col, row, M_cast)
# Convert M into matrix
M <- as.matrix(M)

assertthat::assert_that(all(labels(M)[1][[1]] == rownames))
assertthat::assert_that(all(labels(M)[2][[1]] == colnames))

# Combine occupancy and misclassification models to simulate observed data --------
# c_obs: c_obs[i,j,k',l] are the elements of vector C, and represent the number
# of individuals that were classified as k'. dim: nsite x nsurv x nspec x nyear
c_obs <- array(NA, dim = c(length(unique(all_paratax_df$plotID)),length(unique(all_paratax_df$col_index)),
                           length(colnames),length(unique(all_paratax_df$col_year))),
                 dimnames = list(sort(unique(all_paratax_df$plotID)), #plots
                                 unique(all_paratax_df$col_index), #surveys
                                 colnames, #para_morphs
                                 unique(all_paratax_df$col_year))) #year

# Grab values from casted c_obs array and fill in values in final c_obs array
c_obs_cast   <- all_paratax_df %>%
  reshape2::acast(plotID ~ col_index ~ scimorph_combo ~ col_year,
                  fun.aggregate=sum, fill=-999, value.var = "individualCount")
c_obs_cast[c_obs_cast == -999] <- NA

assertthat::assert_that(sum(c_obs_cast, na.rm=T) == sum(all_paratax_df$individualCount, na.rm=T))

for(plot in dimnames(c_obs_cast)[[1]]) {
  for(surv in dimnames(c_obs_cast)[[2]]) {
    for(morph in dimnames(c_obs_cast)[[3]]) {
      for(year in dimnames(c_obs_cast)[[4]]) {
        c_obs[plot,surv,morph,year] <- c_obs_cast[plot,surv,morph,year]
      }
    }
  }
}
c_obs[is.na(c_obs)] <- 0
rm(c_obs_cast)

assertthat::assert_that(all(labels(c_obs)[3][[1]] == colnames))

# Occupancy array. dim: nsite x nspec x nyear
# Create Z data. Use expert identifications to incorporate partly observed presence
Z.dat <- array(NA, dim = c(length(unique(all_paratax_df$plotID)),
                           length(rownames),length(unique(all_paratax_df$col_year))),
               dimnames = list(sort(unique(all_paratax_df$plotID)), #plots
                               rownames, #expert IDs
                               unique(all_paratax_df$col_year))) #year

# Grab values from casted Z.dat array and fill in values in final Z.dat array
Z.dat_cast <- expert_df %>% 
  mutate(occ = 1) %>%
  reshape2::acast(plotID ~ scientificName ~ col_year,
                  fill=-999, drop=F, value.var = "occ")
Z.dat_cast[Z.dat_cast == -999] <- NA
Z.dat_cast[Z.dat_cast > 0] <- 1

assertthat::assert_that(sum(Z.dat_cast, na.rm=T) == nrow(expert_df %>% distinct(plotID, scientificName, col_year)))

for(plot in dimnames(Z.dat_cast)[[1]]) {
  for(morph in dimnames(Z.dat_cast)[[2]]) {
    for(year in dimnames(Z.dat_cast)[[3]]) {
      Z.dat[plot,morph,year] <- Z.dat_cast[plot,morph,year]
    }
  }
}
rm(Z.dat_cast,plot,morph,year)

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
rm(i,k,l)

# JAGS model --------------------------------------------------------------

# Run model in JAGS. 
JAGSdata <- list(nsite = dim(c_obs)[1], 
                 nsurv = dim(c_obs)[2], 
                 nspec_exp = dim(M)[1], #species ID'd by expert
                 nspec_para = dim(M)[2], #paratax ID's and morphospecies
                 nyear = dim(c_obs)[4], 
                 n = n,
                 alpha = alpha,
                 M = M,
                 c_obs = c_obs,
                 R = diag(rep(1, 5)),
                 Z = Z.dat) #bundle data
JAGSinits <- function(){ list(Z = Z.init) }
nc <- 4 #MCMC chains
ni <- 40000
nt <- ni/2000

# JAGS model
JAGSparams <- c("psi", "lambda", "theta", "Z", "phi", 
                "gamma", "n.occ", "growth", "turnover", 
                "eps_spec", "mu_spec", "Tau_spec", 
                "eps_site", "Tau_site")
cl <- makeCluster(nc)

out <- jags.parfit(cl = cl,
                   data = JAGSdata,
                   params = JAGSparams,
                   model = "occupancy/15_neon_dynamic_multisp_misclass_JAGS.txt",
                   inits = JAGSinits,
                   n.chains = nc,
                   n.adapt = 5000,
                   n.update = 5000,
                   thin = nt,
                   n.iter = ni)
saveRDS(out, "occupancy/script15_jags_out.rds")
out_mcmcsumm <- MCMCsummary(out)
saveRDS(out_mcmcsumm, "occupancy/script15_mcmcsumm_out.rds")

# View JAGS output --------------------------------------------------------

out <- readRDS("occupancy/script15_jags_out.rds")
out_mcmcsumm <- readRDS("occupancy/script15_mcmcsumm_out.rds")

# Did model converge?
hist(out_mcmcsumm$Rhat, breaks=40)
range(out_mcmcsumm$Rhat, na.rm=TRUE)

# Look at high Rhat values
out_mcmcsumm <- rownames_to_column(out_mcmcsumm)
out_mcmcsumm %>% filter(Rhat > 1.4)

MCMCtrace(out, params = paste0('eps_spec\\[2,1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('eps_spec\\[2,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('eps_spec\\[2,3\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('eps_spec\\[2,4\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 

# Look at raw numbers
# theta
theta_summ <- MCMCsummary(out, params = 'theta', round=2)
saveRDS(theta_summ, "occupancy/theta_summ.rds")
theta_summ <- readRDS("occupancy/theta_summ.rds")
hist(theta_summ$Rhat)
range(theta_summ$Rhat, na.rm=TRUE)

theta_summ <- rownames_to_column(theta_summ)
theta_summ %>% filter(Rhat > 1.1) #3774 rows total
theta_summ %>% filter(Rhat > 1.4) 

# Look at some theta posteriors with high Rhat values
MCMCtrace(out, params = paste0('theta\\[5,5\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('theta\\[6,6\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('theta\\[13,13\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
left_join(expert_df %>% select(individualID, expert=scientificName), pinned_df %>% select(individualID, para = scientificName)) %>% filter(scientificName==rownames[10])

MCMCtrace(out, params = paste0('theta\\[5,15\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('theta\\[6,15\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('theta\\[13,15\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(out, params = paste0('theta\\[5,16\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
#MCMCsummary(out, params='theta\\[36,[0-9]+\\]', ISB=F)


theta_df <- data.frame(expert_index = rep(1:dim(M)[1], dim(M)[2]) ,
                       paramorph_index = rep(1:dim(M)[2], each = dim(M)[1]),
                       expert_sciname = rep(rownames, dim(M)[2]),
                       para_morph = rep(colnames, each = dim(M)[1]))  %>% 
    mutate(theta_mean = theta_summ$mean)
#make para_morph a factor to force plotting in order
theta_df$para_morph = factor(theta_df$para_morph, levels=colnames) 

# Plot heatmap of theta values
ggplot(theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_mean)) + 
    geom_tile() +
    scale_fill_gradient(low="darkblue", high="white") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))
#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html


# phi - survival probability
MCMCsummary(out, params = 'phi', round=2)
MCMCtrace(out, params = 'phi', type = 'density', ind = F, pdf=F)
#phi[47] Rhat=1.11. 
#several indices are bimodal

# gamma - colonization probability
MCMCsummary(out, params = 'gamma', round=2)
MCMCtrace(out, params = 'gamma', type = 'density', ind = F, pdf=F)
#gamma[47] Rhat=1.1. same index as above for phi. 
#all densities peak close to 0 - makes sense

# psi - occupancy prob.
MCMCsummary(out, params = 'psi', round=2)
MCMCtrace(out, params = 'psi', type = 'density', ind = F, pdf=F)
# a few Rhat values > 1.1. psi[49,1] has Rhat 3.09, psi[36,1] has Rhat 2.06
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
lambda_summ <- MCMCsummary(out, params = 'lambda', round=2) #takes 5ish min
saveRDS(lambda_summ, "occupancy/lambda_summ.rds")
lambda_summ <- readRDS("occupancy/lambda_summ.rds")
hist(lambda_summ$Rhat)
range(lambda_summ$Rhat)
lambda_summ
range(lambda_summ$mean)
MCMCtrace(out, params = 'lambda', type = 'density', ind = F, pdf=F)
range(lambda_summ$mean)
#a good handful of Rhats > 1.1, as large as 3.32

# Z
Z_summ <- MCMCsummary(out, params = 'Z', round=2)
saveRDS(Z_summ, "occupancy/Z_summ.rds")
#Z_summ <- readRDS("occupancy/theta_summ.rds")
MCMCtrace(out, params = 'Z', type = 'density', ind = F, pdf=F)

Z.dat[is.na(Z.dat)] <- 0
Z.init[is.na(Z.init)] <- 0
Z.prior <- Z.dat+Z.init
plot(Z.prior, Z_summ$mean)
#WORK ON ABOVE PLOT - MAKE SURE THEY PLOT 1-1 ACCURATELY, now it's not accurate
plot(apply(c_obs, c(1,3,4), max, na.rm = TRUE), out$mean$Z)











# n.occ
print(out$summary[grep("n.occ", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2)

# growth
print(out$summary[grep("growth", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
out$mean$growth

# turnover
print(out$summary[grep("turnover", row.names(out$summary)), c(1, 2, 3, 7)], dig = 2) 
out$mean$turnover

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
par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
abline(0,1)
points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")

