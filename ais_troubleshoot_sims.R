

jm_1.0 <- readRDS(paste0("output/simulations/sim_0_jm_full_1.0.rds"))
jm_0.05 <- readRDS(paste0("output/simulations/sim_0_jm_full_0.05.rds"))

jm1.0summ <- MCMCsummary(jm_1.0, param="Theta") %>%
    mutate(fraction = "1.0") %>%
    rownames_to_column()
jm0.05summ <- MCMCsummary(jm_0.05, param="Theta") %>%
    mutate(fraction = "1.0") %>%
    rownames_to_column()

# I would be curious to take two of the extremes (e.g., 0.05 and 1 validation fractions), and
# compare the prior & posterior densities for some of the classification probabilities in Theta.
# I wonder whether the posteriors are different

# 1.0 posterior
post_1.0_chains <- MCMCchains(jm_1.0, params='Theta') 
post_1.0 <- as.data.frame(post_1.0_chains) %>%
    pivot_longer(colnames(post_1.0_chains),names_to="index") %>%
    mutate(model = "1.0")
# 0.05 posterior
post_0.05_chains <- MCMCchains(jm_0.05, params='Theta') 
post_0.05 <- as.data.frame(post_0.05_chains) %>%
    pivot_longer(colnames(post_0.05_chains),names_to="index") %>%
    mutate(model = "0.05")
# Prior
prior_theta <- data.frame(index = rep(colnames(post_0.05_chains),each=nrow(post_0.05_chains)),
                          value = NA,
                          model = "prior") %>%
    as_tibble()
value_vec <- c()
for (i in 1:ncol(alpha)) {
    temp_vec <- as.vector(MCMCpack::rdirichlet(3200, alpha[,i]))
    value_vec <- c(value_vec, temp_vec)
}
prior_theta$value <- value_vec
prior_theta$index <- as.character(prior_theta$index)
                                  
# Bind them together                            
theta_df <- rbind(prior_theta, post_1.0, post_0.05)

# Plot matrix of theta prior and posterior densities
select_ind <- c("Theta[1,1]","Theta[1,2]","Theta[1,3]","Theta[1,4]",
                "Theta[2,1]","Theta[2,2]","Theta[2,3]","Theta[2,4]",
                "Theta[3,1]","Theta[3,2]","Theta[3,3]","Theta[3,4]",
                "Theta[4,1]","Theta[4,2]","Theta[4,3]","Theta[4,4]")
# select_ind <- c("Theta[17,27]","Theta[17,28]","Theta[17,29]","Theta[17,30]",
#                 "Theta[18,27]","Theta[18,28]","Theta[18,29]","Theta[18,30]",
#                 "Theta[19,27]","Theta[19,28]","Theta[19,29]","Theta[19,30]",
#                 "Theta[20,27]","Theta[20,28]","Theta[20,29]","Theta[20,30]")
theta_df %>% 
    filter(index %in% select_ind ) %>%
    ggplot(aes(x=value, y=..scaled.., fill=model)) +
    geom_density( alpha=0.6 ) + 
    facet_wrap( ~ index, scales="free_x") +
    xlab("Theta") + scale_y_continuous(breaks=seq(0, 1, 0.5)) 

# If not, could it be that the prior is too strong, such that the data aren't influencing the posterior that much?


# Questions for max
# What's a good n.eff minimum? I have many parameters with n.eff below 1000.
