# Here we create a Dynamic occupancy model with misclassification using
# disaggregated individual-level data.This model estimates misclassification
# probabilities with NEON Niwot Ridge carabid data 2015-2018

library(dclone) #jags.parfit()
library(MCMCvis)
library(dplyr)
library(ggmcmc)
library(readr)
library(patchwork) #wrap_plots
options(digits = 3)

# Load jags input
source("source/jags_input.R")
jagsinput <- return_jags_input("full")
list2env(jagsinput, .GlobalEnv)
rm(jagsinput)

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
if (!file.exists("output/full_jm.rds")) {
    jags_d <- list(nsite = dim(L)[1],
                   nsurv = dim(L)[2], 
                   nyear = dim(L)[3], 
                   K_exp = dim(alpha)[1], 
                   K_para = dim(alpha)[2],
                   L = L, 
                   alpha = alpha,
                   Ltot = sum(L), 
                   site = y_df$plotID_idx, #needs to be numeric
                   year = as.numeric(as.factor(y_df$col_year)), #needs to be numeric starting from 1
                   # if the individual was labeled by the expert, true ID is known
                   k = y_df$expertID_idx,
                   # for all individuals, we get paratxonomist IDs
                   y = y_df$parataxID_idx,
                   z = z.dat,
                   R = diag(rep(1, 4)))
    JAGSinits <- function(){ list(z = z.init) }
    nc <- 4
    ni <- 20000
    cl <- makeCluster(nc)
    jm <- jags.parfit(cl = cl,
                      data = jags_d,
                      params = c("logit_psi1", "logit_phi", "logit_gamma",
                                 "psi",
                                 "lambda", "Theta",
                                 "eps_site", "eps_spec", "Tau_spec", "Tau_site",
                                 "log_growth", "turnover"),
                      model = "full_dyn_occ_misclass_JAGS.txt",
                      n.chains = nc,
                      n.adapt = 2000,
                      n.update = 2000,
                      thin = ni/1000,
                      n.iter = ni) 
    
    jm_summ <- MCMCsummary(jm)
    theta_summ <- MCMCsummary(jm, params = 'Theta') %>%
      rownames_to_column()
    
    dir.create("output", showWarnings = FALSE)
    saveRDS(jm, "output/full_jm.rds")
    saveRDS(jm_summ, "output/full_jmsumm.rds")
    saveRDS(theta_summ, "output/full_theta_summ.rds") 
}

