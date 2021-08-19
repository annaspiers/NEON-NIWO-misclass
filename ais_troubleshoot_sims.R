

jm_1.0 <- readRDS(paste0("output/simulations/sim_3_jm_full_1.0.rds"))
jm_0.001 <- readRDS(paste0("output/simulations/sim_3_jm_full_0.001.rds"))

jm1.0summ <- MCMCsummary(jm_1.0, param="Theta") %>%
    mutate(fraction = "1.0") %>%
    rownames_to_column()
jm0.001summ <- MCMCsummary(jm_0.001, param="Theta") %>%
    mutate(fraction = "0.001") %>%
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
post_0.001_chains <- MCMCchains(jm_0.001, params='Theta') 
post_0.001 <- as.data.frame(post_0.001_chains) %>%
    pivot_longer(colnames(post_0.001_chains),names_to="index") %>%
    mutate(model = "0.001")
# Prior
prior_theta <- data.frame(index = rep(colnames(post_0.001_chains),each=nrow(post_0.001_chains)),
                          value = NA,
                          model = "prior") %>%
    as_tibble()
value_vec <- c()
for (i in 1:ncol(alpha)) {
    temp_vec <- as.vector(MCMCpack::rdirichlet(dim(post_0.001_chains)[1], alpha[,i]))
    value_vec <- c(value_vec, temp_vec)
}
prior_theta$value <- value_vec
prior_theta$index <- as.character(prior_theta$index)
                                  
# Bind them together                            
theta_df <- rbind(prior_theta, post_1.0, post_0.001)

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




# Heat map of posterior differences between models  
full_theta_df <- data.frame(expert_index = as.factor(rep(1:dim(alpha)[1], dim(alpha)[2])) ,
                            paramorph_index = as.factor(rep(1:dim(alpha)[2], each = dim(alpha)[1])))  %>% 
    mutate(theta_median = jm1.0summ$"50%")

red_theta_df <- data.frame(expert_index = as.factor(rep(1:dim(alpha)[1], dim(alpha)[2])) ,
                           paramorph_index = as.factor(rep(1:dim(alpha)[2], each = dim(alpha)[1])))  %>% 
    mutate(theta_median = jm0.001summ$"50%")

theta_med_diff_df <- red_theta_df %>%
    rename(red_theta_med = theta_median) %>%
    left_join(full_theta_df %>% 
                  dplyr::select(expert_index, paramorph_index, full_theta_med=theta_median)) %>%
    mutate(median_diff = full_theta_med - red_theta_med,
           diff_pos = ifelse(median_diff == 0, NA, ifelse(median_diff > 0.01, T, F)))

ggplot(theta_med_diff_df, aes(x=paramorph_index, y=expert_index, fill= diff_pos)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = c("blue", "red"), "Main > Reduced \nTheta Median") +
    xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
    scale_y_discrete(limits = rev(levels(as.factor(theta_med_diff_df$expert_index))))
cols <- c("-1" = "blue", "0" = "white", "1" = "red")
ggplot(theta_med_diff_df, aes(x=paramorph_index, y=expert_index, fill= median_diff)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red") +
    xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
    scale_y_discrete(limits = rev(levels(as.factor(theta_med_diff_df$expert_index))))



# Visualize dirichlet dist
library(MCMCpack) #rdirichlet()

# Want to be similar to this from original model:
theta_prior <- rdirichlet(1000, c(200, rep(2, 76)))
plot(density(theta_prior[,-1]), xlim=c(0,1), col="black", main=NA) 
lines(density(theta_prior[,1]), col="black")

theta_prior <- rdirichlet(1000, c(500, rep(100, 2)))
lines(density(theta_prior[,-1]), col="red")
lines(density(theta_prior[,1]),col="red")

