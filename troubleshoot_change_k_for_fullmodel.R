# Troubleshooting simulations
# Run original full dynamic occupancy model 
# Then run again with less validation data

library(jagsUI)
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

nc <- 4
ni <- 20000

# Full JAGS model with NEON data ------------------------------------------
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
jm <- jags(data = jags_d,
                    inits = function(){ list(z = z.init) },
                    parameters.to.save = c("psi","lambda", "Theta"), 
                    model.file = "full_dyn_occ_misclass_JAGS.txt",
                    n.chains = nc,
                    n.adapt = 2000, 
                    n.iter= ni,
                    n.burnin = 200,
                    n.thin = ni/1000)
jm_summ <- MCMCsummary(jm)
theta_summ <- MCMCsummary(jm, params = 'Theta') %>%
        rownames_to_column()
    
dir.create("output", showWarnings = FALSE)
saveRDS(jm, "output/full_jm_og.rds")
saveRDS(jm_summ, "output/full_jmsumm_og.rds")
saveRDS(theta_summ, "output/full_theta_summ_og.rds") 



# Full JAGS model with NEON data WITH FEWER SAMPLES VALIDATED ------------------------

# Take existing data used in main model above, and keep only some of the validation data
# 4772 individuals total, 1764 (37%) are validated.
# Let's keep 20% of the validated data: Now 353 (7.4% of total) samples will be validated
y_df$expertID_idx_less <- y_df$expertID_idx #copy the column of verified sample IDs over
alreadyNA2 <- sum(is.na(y_df$expertID_idx_less)) #caluclate number of samples without validation already
needtobeNA2 <- round((1-0.074)*sum(L)) #calculate number of samples that you want without validation 
ind2 <- sample(which(!is.na(y_df$expertID_idx_less)), (needtobeNA2-alreadyNA2), replace=F) 
y_df$expertID_idx_less[ind2]<- NA

jags_d_less <- list(nsite = dim(L)[1],
               nsurv = dim(L)[2], 
               nyear = dim(L)[3], 
               K_exp = dim(alpha)[1], 
               K_para = dim(alpha)[2],
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = y_df$plotID_idx, 
               year = as.numeric(as.factor(y_df$col_year)), 
               k = y_df$expertID_idx_less, #this has fewer validated samples than before
               y = y_df$parataxID_idx,
               z = z.dat,
               R = diag(rep(1, 4)))
jm_less <- jags(data = jags_d_less,
           inits = function(){ list(z = z.init) },
           parameters.to.save = c("psi","lambda", "Theta"), 
           model.file = "full_dyn_occ_misclass_JAGS.txt",
           n.chains = nc,
           n.adapt = 2000, 
           n.iter= ni,
           n.burnin = 2000,
           n.thin = ni/1000)
jm_summ_less <- MCMCsummary(jm_less)
theta_summ_less <- MCMCsummary(jm_less, params = 'Theta') %>%
    rownames_to_column()

dir.create("output", showWarnings = FALSE)
saveRDS(jm_less, "output/full_jm_less.rds")
saveRDS(jm_summ_less, "output/full_jmsumm_less.rds")
saveRDS(theta_summ_less, "output/full_theta_summ_less.rds") 



# Code to visualize results and compare models ----------------------------

jm_og <- readRDS(paste0("output/full_jm_og.rds"))
jm_less <- readRDS(paste0("output/full_jm_less.rds"))

jm_og_summ <- MCMCsummary(jm_og, param="Theta") %>%
    mutate(fraction = "original") %>%
    rownames_to_column()
jm_less_summ <- MCMCsummary(jm_less, param="Theta") %>%
    mutate(fraction = "less") %>%
    rownames_to_column()

# Posterior for original model
post_og_chains <- MCMCchains(jm_og, params='Theta') 
post_og <- as.data.frame(post_og_chains) %>%
    pivot_longer(colnames(post_og_chains),names_to="index") %>%
    mutate(model = "original")
# Posterior for models with less validation data
post_less_chains <- MCMCchains(jm_less, params='Theta') 
post_less <- as.data.frame(post_less_chains) %>%
    pivot_longer(colnames(post_less_chains),names_to="index") %>%
    mutate(model = "less")
# Prior
prior_theta <- data.frame(index = rep(colnames(post_less_chains),each=nrow(post_less_chains)),
                          value = NA,
                          model = "prior") %>%
    as_tibble()
value_vec <- c()
for (i in 1:ncol(alpha)) {
    temp_vec <- as.vector(MCMCpack::rdirichlet(dim(post_less_chains)[1], alpha[,i]))
    value_vec <- c(value_vec, temp_vec)
}
prior_theta$value <- value_vec
prior_theta$index <- as.character(prior_theta$index)

# Bind them together                            
theta_df <- rbind(prior_theta, post_og, post_less)

# Plot matrix of theta prior and posterior densities
# select_ind <- c("Theta[1,1]","Theta[1,2]","Theta[1,3]","Theta[1,4]",
#                 "Theta[2,1]","Theta[2,2]","Theta[2,3]","Theta[2,4]",
#                 "Theta[3,1]","Theta[3,2]","Theta[3,3]","Theta[3,4]",
#                 "Theta[4,1]","Theta[4,2]","Theta[4,3]","Theta[4,4]")
select_ind <- c("Theta[32,32]","Theta[32,33]","Theta[32,34]","Theta[32,35]",
                "Theta[33,32]","Theta[33,33]","Theta[33,34]","Theta[33,35]",
                "Theta[34,32]","Theta[34,33]","Theta[34,34]","Theta[34,35]",
                "Theta[35,32]","Theta[35,33]","Theta[35,34]","Theta[35,35]")
theta_df %>% 
    filter(index %in% select_ind ) %>%
    ggplot(aes(x=value, y=..scaled.., fill=model)) +
    geom_density( alpha=0.6 ) + 
    facet_wrap( ~ index, scales="free_x") +
    xlab("Theta") + scale_y_continuous(breaks=seq(0, 1, 0.5)) 

# Heat map of posterior differences between models  
full_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                            paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                            expert_sciname = rep(rownames, dim(alpha)[2]),
                            para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
    mutate(theta_median = jm_og_summ$"50%")
full_theta_df$para_morph = factor(full_theta_df$para_morph, levels=colnames) #make para_morph a factor to force plotting in order

red_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                           paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                           expert_sciname = rep(rownames, dim(alpha)[2]),
                           para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
    mutate(theta_median = jm_less_summ$"50%")
#make para_morph a factor to force plotting in order
red_theta_df$para_morph = factor(red_theta_df$para_morph, levels=colnames) 

theta_med_diff_df <- red_theta_df %>%
    rename(red_theta_med = theta_median) %>%
    left_join(full_theta_df %>% 
                  dplyr::select(expert_sciname, para_morph, full_theta_med=theta_median)) %>%
    mutate(median_diff = full_theta_med - red_theta_med,
           diff_pos = ifelse(median_diff == 0, NA, ifelse(median_diff > 0.01, T, F)))

ggplot(theta_med_diff_df, aes(x=para_morph, y=expert_sciname, fill= diff_pos)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = c("blue", "red"), "Main > Reduced \nTheta Median") +
    xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
    scale_y_discrete(limits = rev(levels(as.factor(theta_med_diff_df$expert_sciname))))

# Which species have a higher classification accuracy estimated for the original model?
theta_med_diff_df %>% filter(diff_pos==T) %>% arrange(expert_index)
# For the cases where the original model beats the reduced one, what is special about these species? Look at the abundance for each expert_sciname and para_morph in the data
y_df %>%
    group_by(parataxID_idx) %>%
    summarize(n=n()) %>%
    arrange(-n)
y_df %>%
    filter(!is.na(expertID_idx_less)) %>%
    group_by(expertID_idx_less) %>%
    summarize(n=n()) %>%
    arrange(-n) 

# - [ ] I ran the original model on the NEON dataset and ran it again with less validation data. Similarly to the simulations, we see similar Theta estimates returned by these two models. 
# - [ ] removing validation data from neon beetle data yields similar results where estimates from full and reduced models look similar
# - [ ] show binary diff graph
# - [ ] Which species have higher accuracy for full model compared to model with reduced validation data? What’s weird is the original model yields a higher classification estimate than the reduced model for species that the reduced model actually has validation for.
# - [ ] show density comaprison graph
# - [ ] What does this imply? Not clear. Any validation data that the reduced model has, the original model has. However, there is something that leads the reduced model to have a lower theta median for species where it has validation data. Then where the reduced model has no validation data, it has a similar estimate to the original model’s estimate. Doesn’t make much sense to me
# - [ ] Do we see this same pattern with the simulations?
