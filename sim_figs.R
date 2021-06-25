
# library(jagsUI)
# library(MCMCvis)
# library(ggmcmc)
library(dplyr)
library(ggplot2)
# library(tibble) #rownames_to_column
# options(digits = 3)


# Assess model results ----------------------------------------------------

temp = list.files(path = "output/simulations/", pattern = "sim_*")
sims_jm <- data.frame()
for (i in 1:length(temp)) {
    mcmcdf_temp <- readRDS(paste0("output/simulations/",sort(temp, decreasing=T)[i]))
    sims_jm <- rbind(sims_jm, mcmcdf_temp)
}

# Iterations
# full vs reduced
# by species
# by various fractions of validation


# Plot results ------------------------------------------------------------

# Create new columns for plotting 
# binary column specifying whether estimate 95% CI captures true value
# numeric column specifying how far estimate median is from true value
sims_jm_plot <- sims_jm %>%
    mutate(capt95 = ifelse(true<`2.5%`,0,ifelse(true>`97.5%`,0,1)),
           diff = `50%` - true ) # positive means estimated median is greater than true value

# 1) Coverage
# % of runs where 95% CI captures the true value
# displayed as a table
# Brett says I may need to do a million reps?
sims_jm_plot %>%
    group_by(param, fraction) %>%
    summarize(capt_by95CI = sum(capt95)/n()) %>%
    pivot_wider(names_from=fraction, values_from=capt_by95CI)


# 2) Distance from true value
# Difference between estimate median and true value
# Brett says this is nice
sims_jm_plot %>%
    mutate(fraction = as.factor(fraction)) %>%
    ggplot(aes(x=fraction, y=diff)) + 
    #scale_y_continuous(trans='log10') +
    geom_violin(scale="area") +
    facet_grid(cols=vars(param))



# Calculate validation metrics --------------------------------------------

# temporary code since i didn'ts ave the simulations correctly
sims_jm_val <- sims_jm_plot %>% 
#AIS save this for later. hopefully don't have to do this bc it would take a long time to run the simulations again with a percentage of paratax IDs withheld








    
    










































# 1) Plot occupancy estimate with CI's along with true parameter value as a line for each species. Do for full and reduced model. What I'm imagining is overlapping the posterior distribution for each fraction of validation data simulation. Then plot a straight line over the true occupancy/demographic rate 

# Psi
psi_true_1site <- data.frame(true=psi[2,], species=as.factor(1:K_val)) 
sims_jm %>%
    as_tibble %>%
    #filter by psi
    filter( grepl('psi', rowname) ) %>%
    separate("rowname", into = c("site", "species"), sep = ",") %>%
    mutate_at(c('site', 'species'), readr::parse_number) %>%
    mutate(species = as.factor(species),
           fraction = as.factor(fraction)) %>%
    # AIS how to combine all of the sites together for each species? I could take their average (average 2.5%, 50%, 97.5%). For now, I  select just one site.
    filter(site==2) %>%
    ggplot(aes(fraction, `50%`, group=site)) + 
    geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="gray40", alpha=.25, size=1.5) +
    geom_point(size = 3, color="gray20") +
    geom_hline(data=psi_true_1site, aes(yintercept=true), col="red") +
    facet_wrap(~species) + 
    xlab("Fraction of validated samples") + 
    ylab("Occupancy probability") + 
    theme_minimal()

# Lambda
lambda_true_1site <- data.frame(true=lambda[2,], species=as.factor(1:K_val)) 
sims_jm %>%
    as_tibble %>%
    #filter by psi
    filter( grepl('lambda', rowname) ) %>%
    separate("rowname", into = c("site", "species"), sep = ",") %>%
    mutate_at(c('site', 'species'), readr::parse_number) %>%
    mutate(species = as.factor(species),
           fraction = as.factor(fraction)) %>%
    filter(site==2) %>%
    ggplot(aes(fraction, `50%`, group=site)) + 
    geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="gray40", alpha=.25, size=1.5) +
    geom_point(size = 3, color="gray20") +
    geom_hline(data=lambda_true_1site, aes(yintercept=true), col="red") +
    facet_wrap(~species) + 
    xlab("Fraction of validated samples") + 
    ylab("Encounter rate") + 
    theme_minimal()

# Theta
theta_true <- data.frame(K_val_col=rep(1:K_val,K_imp), 
                         K_imp_col=rep(1:K_imp,each=K_val), 
                         true = c(Theta)) 
sims_jm %>%
    as_tibble %>%
    filter( grepl('Theta', rowname) ) %>%
    separate("rowname", into = c("K_val_col", "K_imp_col"), sep = ",") %>%
    mutate_at(c('K_val_col', 'K_imp_col'), readr::parse_number) %>%
    mutate(K_val_col = as.factor(K_val_col),
           K_imp_col = as.factor(K_imp_col),
           fraction = as.factor(fraction)) %>%
    ggplot(aes(fraction, `50%`)) + 
    geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="red", alpha=.25, size = 1.5) +
    geom_point(size=.75, color="gray20") +
    geom_hline(data=theta_true, aes(yintercept=true), col="red") +
    facet_grid(rows=vars(K_val_col), cols=vars(K_imp_col)) + 
    xlab("Fraction of validated samples") + 
    ylab("Classification probability") + 
    theme_minimal()


# 2) Plot encounter rate estimate with CI's along with true parameter value as a line for each species. Do for full and reduced model

# 3) Report on accuracy, recall, etc. for full and reduced model. Is there evidence that the reduced model is less effective for prioritizing confirmation effort than the proposed model. Address this by calculating accuracy, recall, etc. for full vs reduced model for each simulation. This will speak to reviewer feedback: "it was not clear that this reduced model was markedly worse for the purpose of prioritizing confirmation effort. Were patterns in the theta matrix markedly different for the reduced model (I see Fig 4, but more broadly)?"

# 4) make it more explicitly ecologically meaningful: Discuss whether the credible interval contains the true value. To what degree does partial supervision improve estimation of psi.lambda with varying levels of confirmation effort? Do the same for accuracy, recall, etc. Consider how sensitive the model output is to the fraction validated - Maybe the # IDed by the expert doesn’t affect the occupancy estimate when the occupancy is high, but does when occupancy is low….

# After sending questions to max
# create R and JAGS scripts for reduced model
# create add validation metrics script - see validation.R as model






# MISC
# AIS questions:
# Should we allocate validation randomly across samples or put a weight on samples that are not a validated state? For now, I'll just do random because that's easier
# Should the same individuals not validated when we take 0.5 of k also not be validated when we take 0.2 of k? That is, as we decrease the fraction of k that informs the model, should we keep the same individuals unvalidated, or randomly select again? To make comparing results across varying k fractions more fair, I think we should take subsets of already unvalidated individuals (e.g if indiviual x is not validated when we take 0.8 of k, then individual x should be not validated for all smaller fractions of k)
#should we consider varying accuracy of imperfect classifier? probably not. That is set randomly in priors
# Consider expanding to dynamic model to calculate colonization/extinction. I'm worried just a single season doesn't address reviewer concerns enough
# This simulation doesn't speak to this: "What constitutes a usefully-sized expert-verified sample, how to allocate expert effort across species/sites/particular samples, the choice of priors for theta, specification of/sensitivity to the assumed count distribution". How do we respond to reviewer?



# SCRAP NOTES 

# See Brett's tutorial: https://github.com/EBIO5460Fall2018/Class-materials/blob/master/12_4_sim_multilevel.md
