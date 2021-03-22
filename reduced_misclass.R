# Here we create a reduced model estimating only misclassification probabilities
# without occupancy parameters using NEON Niwot Ridge carabid data 2015-2018

library(dclone) #jags.parfit()
library(MCMCvis)
library(MCMCpack)
library(dplyr)
library(patchwork)
library(ggmcmc)
library(readr)
library(ggplot2)

# Load jags input
source("source/jags_input.R")
jagsinput <- return_jags_input("reduced")
list2env(jagsinput, .GlobalEnv)
rm(jagsinput)

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
if (!file.exists("output/reduced_jm.rds")) {
  jags_d <- list(K_exp = dim(alpha)[1], 
                 K_para = dim(alpha)[2],
                 alpha = alpha,
                 Ltot = sum(L), 
                 # if the individual was labeled by the expert, true ID is known
                 k = y_df$expertID_idx,
                 # for all individuals, we get paratxonomist IDs
                 y = y_df$parataxID_idx)
  JAGSinits <- function(){}
  nc <- 4
  ni <- 20000
  cl <- makeCluster(nc)
  jm <- jags.parfit(cl = cl,
                    data = jags_d,
                    params = c("Theta"),
                    model = "reduced_misclass_JAGS.txt",
                    n.chains = nc,
                    n.adapt = 2000,
                    n.update = 2000,
                    thin = ni/1000,
                    n.iter = ni) 
  
  jm_summ <- MCMCsummary(jm)
  red_theta_summ <- MCMCsummary(jm, params = 'Theta') %>%
      rownames_to_column()
  
  dir.create("output", showWarnings = FALSE)
  saveRDS(jm, "output/reduced_jm.rds")
  saveRDS(jm_summ, "output/reduced_jmsumm.rds")
  saveRDS(red_theta_summ, "output/reduced_theta_summ.rds")
}

