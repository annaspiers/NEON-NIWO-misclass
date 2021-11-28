# Here we create a Dynamic occupancy model with misclassification using
# disaggregated individual-level data.This model estimates misclassification
# probabilities with NEON Niwot Ridge carabid data 2015-2018

library(tidyr) #uncount()
library(tibble) #rownames_to_column()
library(dclone) #jags.parfit()
library(MCMCvis)
library(dplyr)
library(patchwork)
library(ggmcmc)
library(readr)
library(ggplot2)
library(gridExtra) #grid.arrange()
library(cowplot) #get_legend()

# Load jags input for full model
source("source/jags_input.R")
jagsinput <- return_jags_input("full")
list2env(jagsinput, .GlobalEnv)
L_full <- L
y_df_full <- y_df
rm(jagsinput, L, y_df, return_jags_input)

# To validate models, withold parataxonomist IDs from individuals that have expert IDs
perc_withld <- 0.2 #percent individuals with expert IDs that will have paratax ID withheld
set.seed(49)
# Select rows at random that will have expert ID withheld
rows_withld <- y_df_full %>%
  rownames_to_column() %>%
  #filter to data that has expert ID, as we need the true expert ID for validation after running the model
  filter(!is.na(expertID_idx)) %>% #we need true expert ID's available for verification
  slice_sample(prop=perc_withld) %>%
  pull(rowname) #select vector of row indices that will be withheld

# Assign expertIDx in selected rows to NA
y_df_val <- y_df_full %>%
  rownames_to_column() %>%
  mutate(expertID_idx_val = ifelse(rowname %in% rows_withld, NA, expertID_idx)) 

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
if (!file.exists("output/val_full_jm.rds")) {
  JAGSinits <- function(){ list(z = z.init) }
  nc <- 4
  ni <- 20000
  na <- 2000
  cl <- makeCluster(nc)
  
  jags_full <- list(K_exp = dim(alpha)[1],
                    K_para = dim(alpha)[2],
                    nsite = dim(L_full)[1],
                    nsurv = dim(L_full)[2],
                    nyear = dim(L_full)[3],
                    alpha = alpha,
                    L = L_full,
                    Ltot = sum(L_full),
                    site = y_df_val$plotID_idx,
                    year = as.numeric(as.factor(y_df_val$col_year)),
                    k = y_df_val$expertID_idx_val, #all expert IDs known those assigned NA for validation
                    y = y_df_val$parataxID_idx, #all individuals have known parataxID
                    z = z.dat,
                    R = diag(rep(1, 4)))
  jm_full <- jags.parfit(cl = cl,
                    data = jags_full,
                    params = c("Theta", "k"), #we only need these two outputs for validation
                   model = "full_dyn_occ_misclass_JAGS.txt",
                    inits = JAGSinits,
                    n.chains = nc,
                    n.adapt = na,
                    n.update = na,
                    thin = ni/1000,
                    n.iter = ni)
  #jm_full <- jags(data = jags_full,
  #                       inits = JAGSinits,
  #                       parameters.to.save = c("Theta", "k"), #we only need these two outputs for validation
  #                       model.file = "full_dyn_occ_misclass_JAGS.txt",
  #                       n.chains = nc,
  #                       n.adapt = na,
  #                       #n.update = 2000,
  #                       #thin = ni/1000,
  #                       n.iter = ni,
  #                       n.burnin = na)
  #jm_full <- jags(data = jags_d_full,
  #                inits = function(){ list(z = z.init) },
  #                parameters.to.save = c("psi","lambda", "Theta", "k"),
  #                model.file = "sim_full_JAGS.txt",
  #                n.chains = nc,
  #                n.adapt = na,
  #                n.iter= ni,
  #                n.burnin = na)
  
  dir.create("output", showWarnings = FALSE)
  saveRDS(jm_full, "output/val_full_jm.rds")
  saveRDS(jm_full, "val_full_jm.rds")
}

# View JAGS output --------------------------------------------------------

full_val_jm <- readRDS("output/val_full_jm.rds")

### Compare the validation results
# multiclass classification problem

# Filter to validation data
holdout <- y_df_val %>% 
  rename(full_individual=rowname) %>%
  filter(!is.na(expertID_idx)) %>%
  dplyr::select(full_individual, trueID_idx=expertID_idx, trueID=expertID,
                predID_idx=parataxID_idx, predID=parataxID)

# Posterior draws of predicted classifications
full_y_out <- MCMCchains(full_val_jm, params = 'y') %>%
  as_tibble() %>%
  dplyr::select(matches(paste0("y\\[",holdout$full_individual,"\\]"))) %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="full_individual",values_to="predID_idx") %>%
  mutate_at(c('full_individual'), readr::parse_number) %>%
  mutate(full_individual=as.character(full_individual)) %>%
  left_join(holdout %>% dplyr::select(full_individual, trueID_idx)) %>%
  mutate(match = (predID_idx==trueID_idx))

# Generate confusion matrices for each posterior draw ---------------------
all_combos <- expand.grid(draw = 1:max(full_y_out$draw), 
                          est_impID_idx = 1:max(y_df_full$parataxID_idx), #predID_idx to est_impID_idx
                          trueID_idx = 1:max(y_df_full$parataxID_idx, 
                                             na.rm = TRUE))  %>%
  as_tibble 

full_cm <- full_y_out %>%
  filter(draw < max_draw) %>%
  count(draw, predID_idx, trueID_idx) %>%
  full_join(all_combos) %>% # fill in implicit zeros
  reshape2::acast(draw ~ trueID_idx ~ predID_idx, 
                  value.var = "n", 
                  fill = 0)

assertthat::assert_that(!any(is.na(full_cm)))                # no NA vals
assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices

get_metrics <- function(confusion_matrix) {
  # confusion_matrix is a (true, pred) square matrix
  true_positives <- diag(confusion_matrix)
  false_positives <- colSums(confusion_matrix) - diag(confusion_matrix)
  false_negatives <- rowSums(confusion_matrix) - diag(confusion_matrix)
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1 <- 2 * (precision * recall) / (precision + recall)
  tibble(parataxID_idx = 1:length(f1),
         precision = precision, 
         recall = recall, 
         f1 = f1)
}

full_metrics <- apply(full_cm, 1, get_metrics) %>% # number of species x max_draws (see full_cm dimensions)
  bind_rows(.id = "draw") 

# compare distribution of macro F1 (the average of species-specific F1 scores)
full_metrics %>%
  group_by(draw) %>%
  summarize(macro_f1 = mean(f1, na.rm = TRUE)) %>%
  ggplot(aes(macro_f1, fill = model)) + 
  geom_density(alpha = .5)
  

## Accuracy
# Overall (across all species)
full_y_out %>% 
  count(match) %>%
  filter(match==T) %>%
  pull(n) / nrow(full_y_out)

reduced_y_out %>% 
  count(match) %>%
  filter(match==T) %>%
  pull(n) / nrow(reduced_y_out)

full_y_out %>% 
  group_by(draw) %>%
  summarize(accuracy=mean(match)) %>%
  ggplot(aes(x=accuracy)) +
  geom_density()

# By species 
full_y_out %>% 
  count(trueID_idx, match) %>%
  filter(match==T) %>%
  left_join(holdout %>% count(trueID_idx) %>% rename(num_inds=n)) %>%
  mutate(acc = n / (num_inds * nrow(full_y_out)/nrow(holdout)))

full_y_out %>% 
  group_by(draw, trueID_idx) %>%
  summarize(accuracy=mean(match)) %>%
  mutate(trueID = rownames[trueID_idx]) %>%
  ggplot(aes(x=accuracy)) +
  geom_density() +
  facet_wrap(~trueID)


# Holdout log likelihood (log probabilities from the categorical distribution) 
# (could use dcat() function: https://www.rdocumentation.org/packages/mvc/versions/1.3/topics/dcat)
theta_summ_full_val <- MCMCchains(full_val_jm, params = 'Theta') %>% 
  as_tibble() %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="index",values_to="value") %>%
  separate("index", into = c("trueID_idx", "predID_idx"), sep = ",") %>%
  mutate_at(c('trueID_idx', 'predID_idx'), readr::parse_number) 

full_loglik <- holdout %>%
  count(trueID_idx, predID_idx) %>%
  left_join(theta_summ_full_val) %>%
  group_by(draw) %>%
  summarize(log_lik = sum(n * log(value)))


full_loglik %>%
  ggplot(aes(log_lik)) + 
  geom_density() + 
  xlab("Holdout log-likelihood") + 
  ylab("Density")

# create validation.png: accuracy, f1 score, holdout log-lik

acc <- full_y_out %>% 
  group_by(draw) %>%
  summarize(accuracy=mean(match)) %>%
  ggplot(aes(x=accuracy)) +
  geom_density(alpha = .6) + 
  xlab("Accuracy") +
  ylab("Density") +
  theme(legend.position = "none")     

f1 <- full_metrics %>%
  group_by(draw) %>%
  summarize(macro_f1 = mean(f1, na.rm = TRUE)) %>%
  ggplot(aes(macro_f1)) + 
  geom_density(alpha = .6) +
  xlab("F1 score") +
  ylab(NULL) + 
  theme(legend.position = "none")     

hll <- full_loglik %>%
  ggplot(aes(log_lik)) + 
  geom_density(alpha = .6) + 
  xlab("Holdout log-likelihood") + 
  ylab(NULL) + 
  theme(legend.position = "none")   

legend <- get_legend( full_loglik %>%
  ggplot(aes(log_lik)) + 
  geom_density(alpha=0.6) +
    theme(legend.key.size = unit(0.8, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size=10)) )

grid.arrange(acc, f1, hll, legend, ncol = 4) 
val <- arrangeGrob(acc, f1, hll, legend, ncol = 4)
ggsave("figures/validation.png", val, width = 8)

