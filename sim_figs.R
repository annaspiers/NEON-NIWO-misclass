library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)
library(tidyr) #pivot_wider
theme_set(theme_minimal())

# Assess model results ----------------------------------------------------

temp <- list.files(path = "output/simulations", pattern = "sim_out_*", full.name=TRUE)
sims_jm_list <- lapply(temp, readRDS) #takes 1 min
sims_jm <- rbindlist(sims_jm_list)
rm(temp, sims_jm_list)

# remaining <- setdiff(1:max(sims_jm$dataset),sims_jm$dataset) 

range(sims_jm$dataset) #number of iterations
length(unique(sims_jm$dataset))
sims_jm %>% distinct(model) #full vs reduced
sims_jm %>% distinct(fraction) #by various fractions of validation
sims_jm %>% group_by(param,model) %>% summarize(n=n()) 
sims_jm %>% group_by(param,model, fraction) %>% summarize(n=n()) %>% filter(param=="Theta") #why is n=60 for lambda/psi? 2species x 30 sites

# Plot results ------------------------------------------------------------

# Create new columns for plotting 
# binary column specifying whether estimate 95% CI captures true value
# numeric column specifying how far estimate median is from true value
sims_jm_plot <- sims_jm %>%
    mutate(capt95 = ifelse(true<`2.5%`,0,ifelse(true>`97.5%`,0,1)),
           diff = `50%` - true, # positive means estimated median is greater than true value
           ci95_width = `97.5%` - `2.5%`) # 95% CI width
rm(sims_jm)

# 1) Coverage
# % of runs where 95% CI captures the true value
# AIS Brett says I may need to do a million reps?
# Table
sims_jm_plot %>%
    group_by(param, fraction, model) %>%
    summarize(capt_by95CI = sum(capt95)/n()) %>%
    pivot_wider(names_from=fraction, values_from=capt_by95CI)
sims_jm_plot %>%
    filter(param=="Theta") %>%
    group_by(param, fraction, model) %>%
    summarize(capt_by95CI = sum(capt95)/n()) %>%
    group_by(param, model, capt_by95CI) %>%
    summarize(n=n()) #few capt_by95CI values because there are only 6 %ages possible for coverage
# since Theta has 6 different combos (2real species x 3 imperfect species)
# Plot
sims_jm_plot %>%
    group_by(param, fraction, dataset, model) %>%
    summarize(capt_by95CI = sum(capt95)/n()) %>%
    mutate(param_model = paste0(param,"_",model)) %>%
    #filter(capt_by95CI != 1) %>%
    ggplot(aes(x=fraction, y=capt_by95CI)) + 
    geom_point(size=0.3,alpha=0.2) +
    geom_smooth() +
    facet_wrap(vars(param_model)) +
    labs(y="% of Estimates within 95% CI", x="Proportion of data validated")
ggsave("figures/coverage.png")
# Just Theta plot
sims_jm_plot %>%
    filter(param=="Theta") %>%
    group_by(param, fraction, dataset, model) %>%
    summarize(capt_by95CI = sum(capt95)/n()) %>%
    mutate(param_model = paste0(param,"_",model)) %>%
    filter(capt_by95CI != 1) %>%
    ggplot(aes(x=fraction, y=capt_by95CI)) + 
    geom_point(size=0.1,alpha=0.5) +
    geom_smooth() +
    facet_wrap(vars(param_model)) +
    labs(y="% of Estimates within 95% CI", x="Proportion of data validated")
ggsave("figures/coverage_theta.png")


# 2) Distance from true value
# Difference between estimate median and true value
sims_jm_plot %>%
    filter(param=="Theta") %>%
    ggplot(aes(x=fraction, y=abs(diff))) + 
    geom_point(size=0.1,alpha=0.01) +
    geom_smooth() +
    facet_grid(cols=vars(model)) +
    labs(y="Estimated Median - True Value", x="Proportion of data validated")
ggsave("figures/dist_from_true_value.png")


# 3) Plot proportion of data labeled vs 95% CI width - plot full vs reduced model
# Theta
sims_jm_plot %>%
    filter(param=="Theta") %>% 
    #head(100) %>%
    mutate(param_name = paste0(param,"[",i,",",j,"]")) %>%
    ggplot(aes(x=fraction, y=ci95_width, col=model)) +
    geom_point(size=.3, alpha=0.05, position = "jitter") +
    geom_smooth() +
    facet_wrap(vars(param_name)) +
    labs(x="Proportion of data labeled", y="95% CI width")
ggsave("figures/95CIwidth_theta.png")
# Psi
sims_jm_plot %>%
    filter(param=="psi") %>% #filter(param=="psi" & fraction==pull(sims_jm_plot[1,"fraction"])) %>% 
    ggplot(aes(x=fraction, y=ci95_width)) +
    geom_point(size=.5, alpha=0.05) +
    geom_smooth()  +
    labs(x="Proportion of data labeled", y="95% CI width",title="Psi")
ggsave("figures/95CIwidth_psi.png")
# Lambda
sims_jm_plot %>%
    filter(param=="lambda") %>% 
    ggplot(aes(x=fraction, y=ci95_width)) +
    geom_point(size=.5, alpha=0.5) +
    geom_smooth()  +
    labs(x="Proportion of data labeled", y="95% CI width",title="Lambda")
ggsave("figures/95CIwidth_lambda.png")




# 4) Plot 95% CI width for full vs reduced model with 1-1 line
sims_jm_plot %>%
    filter(param=="Theta") %>%
    dplyr::select(param, i, j, fraction, dataset, model, ci95_width) %>% #select only columns with shared variables
    pivot_wider(names_from=model, values_from=ci95_width) %>%
    ggplot(aes(x=reduced, y=full)) +
    geom_point(aes(col=fraction),size=0.5) +
    scale_color_viridis(alpha=0.8) +
    geom_abline(aes(intercept=0, slope=1)) +
    labs(x="95% CI width: misclassification only", y="95% CI width: occupancy misclassification") 
ggsave("figures/95CI_fullvred.png")




# Validation metrics ------------------------------------------------------

# Assess model results ----------------------------------------------------

temp = list.files(path = "output/simulations/", pattern = "sim_val_*", full.name=TRUE)
sims_val_list <- lapply(temp, readRDS) #takes 1 min
sims_val <- rbindlist(sims_val_list)
rm(temp, sims_val_list)

range(sims_val$dataset) #number of iterations
length(unique(sims_val$dataset))
sims_val %>% distinct(model) #full vs reduced
sims_val %>% distinct(fraction) #by various fractions of validation

# All metrics together
sims_val %>%
    pivot_longer(cols=c(accuracy,precision,recall,f1), names_to="metric") %>%
    ggplot(aes(x=fraction,y=value)) +
    geom_point(size=0.5,alpha=0.02) +
    stat_smooth() +
    facet_grid(cols=vars(metric)) +
    labs(x="Proportion of data labeled", y="Metric Value")
ggsave("figures/val_all.png")
# Accuracy
sims_val %>%
    ggplot(aes(x=fraction, y=accuracy)) +
    geom_point(size=0.5,alpha=0.03) +
    geom_smooth() +
    labs(x="Proportion of data labeled", y="Accuracy")
ggsave("figures/val_acc.png")
# F1
sims_val %>%
    ggplot(aes(x=fraction, y=f1)) +
    geom_point(size=0.5,alpha=0.03) +
    geom_smooth() +
    labs(x="Proportion of data labeled", y="F1 Score")
ggsave("figures/val_f1.png")
# Precision
sims_val %>%
    ggplot(aes(x=fraction, y=precision)) +
    geom_point(size=0.5,alpha=0.03) +
    geom_smooth() +
    labs(x="Proportion of data labeled", y="Precision")
ggsave("figures/val_prec.png")
# Recall
sims_val %>%
    ggplot(aes(x=fraction, y=recall)) +
    geom_point(size=0.5,alpha=0.03) +
    geom_smooth() +
    labs(x="Proportion of data labeled", y="Recall")
ggsave("figures/val_recall.png")


# # 1) Plot occupancy estimate with CI's along with true parameter value as a line for each 
# species. Do for full and reduced model. What I'm imagining is overlapping the posterior 
# distribution for each fraction of validation data simulation. Then plot a straight line over
# the true occupancy/demographic rate 
# 
# # Psi
# psi_true_1site <- data.frame(true=psi[2,], species=as.factor(1:K_val))
# sims_jm %>%
#     as_tibble %>%
#     #filter by psi
#     filter( grepl('psi', rowname) ) %>%
#     separate("rowname", into = c("site", "species"), sep = ",") %>%
#     mutate_at(c('site', 'species'), readr::parse_number) %>%
#     mutate(species = as.factor(species),
#            fraction = as.factor(fraction)) %>%
#     # AIS how to combine all of the sites together for each species? I could take 
#           their average (average 2.5%, 50%, 97.5%). For now, I  select just one site.
#     filter(site==2) %>%
#     ggplot(aes(fraction, `50%`, group=site)) +
#     geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="gray40", alpha=.25, size=1.5) +
#     geom_point(size = 3, color="gray20") +
#     geom_hline(data=psi_true_1site, aes(yintercept=true), col="red") +
#     facet_wrap(~species) +
#     xlab("Fraction of validated samples") +
#     ylab("Occupancy probability") +
#     theme_minimal()
# 
# # Lambda
# lambda_true_1site <- data.frame(true=lambda[2,], species=as.factor(1:K_val)) 
# sims_jm %>%
#     as_tibble %>%
#     #filter by psi
#     filter( grepl('lambda', rowname) ) %>%
#     separate("rowname", into = c("site", "species"), sep = ",") %>%
#     mutate_at(c('site', 'species'), readr::parse_number) %>%
#     mutate(species = as.factor(species),
#            fraction = as.factor(fraction)) %>%
#     filter(site==2) %>%
#     ggplot(aes(fraction, `50%`, group=site)) + 
#     geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="gray40", alpha=.25, size=1.5) +
#     geom_point(size = 3, color="gray20") +
#     geom_hline(data=lambda_true_1site, aes(yintercept=true), col="red") +
#     facet_wrap(~species) + 
#     xlab("Fraction of validated samples") + 
#     ylab("Encounter rate") + 
#     theme_minimal()
# 
# # Theta
# theta_true <- data.frame(K_val_col=rep(1:K_val,K_imp), 
#                          K_imp_col=rep(1:K_imp,each=K_val), 
#                          true = c(Theta)) 
# sims_jm %>%
#     as_tibble %>%
#     filter( grepl('Theta', rowname) ) %>%
#     separate("rowname", into = c("K_val_col", "K_imp_col"), sep = ",") %>%
#     mutate_at(c('K_val_col', 'K_imp_col'), readr::parse_number) %>%
#     mutate(K_val_col = as.factor(K_val_col),
#            K_imp_col = as.factor(K_imp_col),
#            fraction = as.factor(fraction)) %>%
#     ggplot(aes(fraction, `50%`)) + 
#     geom_linerange(aes(ymin=`2.5%`, ymax=`97.5%`), color="red", alpha=.25, size = 1.5) +
#     geom_point(size=.75, color="gray20") +
#     geom_hline(data=theta_true, aes(yintercept=true), col="red") +
#     facet_grid(rows=vars(K_val_col), cols=vars(K_imp_col)) + 
#     xlab("Fraction of validated samples") + 
#     ylab("Classification probability") + 
#     theme_minimal()


# 2) Plot encounter rate estimate with CI's along with true parameter value as a line for each species. Do for full and reduced model

# 3) Report on accuracy, recall, etc. for full and reduced model. Is there evidence that the reduced model is less effective for 
# prioritizing confirmation effort than the proposed model. Address this by calculating accuracy, recall, etc. for full vs reduced 
# model for each simulation. This will speak to reviewer feedback: "it was not clear that this reduced model was markedly worse for 
# the purpose of prioritizing confirmation effort. Were patterns in the theta matrix markedly different for the reduced model (I see 
# Fig 4, but more broadly)?"

# 4) make it more explicitly ecologically meaningful: Discuss whether the credible interval contains the true value. To what 
# degree does partial supervision improve estimation of psi.lambda with varying levels of confirmation effort? Do the same for 
# accuracy, recall, etc. Consider how sensitive the model output is to the fraction validated - Maybe the # IDed by the expert 
# doesn’t affect the occupancy estimate when the occupancy is high, but does when occupancy is low….


# AIS questions:
# Should we allocate validation randomly across samples or put a weight on samples that are not a validated state? 
# For now, I'll just do random because that's easier
# should we consider varying accuracy of imperfect classifier? probably not. That is set randomly in priors
# This simulation doesn't speak to this: "What constitutes a usefully-sized expert-verified sample, how to 
#allocate expert effort across species/sites/particular samples, the choice of priors for theta, specification 
#of/sensitivity to the assumed count distribution". How do we respond to reviewer?
