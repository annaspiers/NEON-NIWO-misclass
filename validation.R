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

# Load jags input for full model
source("source/jags_input.R")
jagsinput <- return_jags_input("full")
list2env(jagsinput, .GlobalEnv)
L_full <- L
rm(jagsinput, L)

# Load jags input for reduced model
jagsinput <- return_jags_input("reduced")
list2env(jagsinput, .GlobalEnv)
L_reduced <- L
rm(jagsinput, L)

# To validate models, withold parataxonomist IDs from individuals that have expert IDs
perc_withld <- 0.2 #percent individuals with expert IDs that will have paratax ID withheld
# Select rows at random that will have paratax ID withheld
set.seed(49)
rows_withld <- y_df %>%
  rownames_to_column() %>%
  filter(!is.na(expertID_idx)) %>%
  slice_sample(prop=perc_withld) %>%
  pull(rowname) #select vector of row indices that will be withheld
# Assign parataxID_idx in selected rows to NA
y_df <- y_df %>%
  rownames_to_column() %>%
  mutate(parataxID_idx_val = ifelse(rowname %in% rows_withld, NA, parataxID_idx)) 
#saveRDS(y_df,"data/validation_df.rds")

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
JAGSinits <- function(){ list(z = z.init) }
nc <- 4
ni <- 20000
cl <- makeCluster(nc)

jags_full <- list(K_exp = dim(alpha)[1],
                  K_para = dim(alpha)[2],
                  nsite = dim(L_full)[1],
                  nsurv = dim(L_full)[2],
                  nyear = dim(L_full)[3],
                  alpha = alpha,
                  L = L_full,
                  Ltot = sum(L_full),
                  site = y_df$plotID_idx,
                  year = as.numeric(as.factor(y_df$col_year)),
                  k = y_df$expertID_idx,
                  y = y_df$parataxID_idx_val,
                  z = z.dat,
                  R = diag(rep(1, 4)))
jm_full <- jags.parfit(cl = cl,
                  data = jags_full,
                  params = c("Theta", "y"), #we only need these two outputs for validation
                  model = "full_dyn_occ_misclass_JAGS.txt",
                  n.chains = nc,
                  n.adapt = 2000,
                  n.update = 2000,
                  thin = ni/1000,
                  n.iter = ni)
dir.create("output", showWarnings = FALSE)
jm_summ <- MCMCsummary(jm_full)
saveRDS(jm_full, "output/val_full_jm.rds")
saveRDS(jm_summ, "output/val_full_jmsumm.rds")

jags_reduced <- list(K_exp = dim(alpha)[1],
                     K_para = dim(alpha)[2],
                     alpha = alpha,
                     Ltot = sum(L_reduced),
                     k = y_df %>%
                       filter(!is.na(expertID_idx)) %>%
                       pull(expertID_idx),
                     y = y_df %>%
                       filter(!is.na(expertID_idx)) %>%
                       pull(parataxID_idx_val))
jm_reduced <- jags.parfit(cl = cl,
                       data = jags_reduced,
                       params = c("Theta", "y"),
                       model = "reduced_misclass_JAGS.txt",
                       n.chains = nc,
                       n.adapt = 2000,
                       n.update = 2000,
                       thin = ni/1000,
                       n.iter = ni)
jm_summ <- MCMCsummary(jm_reduced)
saveRDS(jm_reduced, "output/val_reduced_jm.rds")
saveRDS(jm_summ, "output/val_reduced_jmsumm.rds")

# View JAGS output --------------------------------------------------------

full_val_jm <- readRDS("output/val_full_jm.rds")
full_val_jmsumm <- readRDS("output/val_full_jmsumm.rds")
reduced_val_jm <- readRDS("output/val_reduced_jm.rds")
reduced_val_jmsumm <- readRDS("output/val_reduced_jmsumm.rds")

# Did models converge?
hist(full_val_jmsumm$Rhat, breaks=40)
hist(reduced_val_jmsumm$Rhat, breaks=40)


### Compare the validation results
# multiclass classification problem

# Filter to validation data
y_df <- readRDS("data/validation_df.rds")
holdout <- y_df %>% 
  rename(full_individual=rowname) %>%
  filter(!is.na(expertID_idx)) %>%
  rownames_to_column() %>%
  rename(reduced_individual=rowname) %>%
  filter(is.na(parataxID_idx_val)) %>%
  select(full_individual, reduced_individual, trueID_idx=parataxID_idx) %>%
  as_tibble()

# Posterior draws of predicted classifications
full_y_out <- MCMCchains(full_val_jm, params = 'y') %>%
  as_tibble() %>%
  select(matches(paste0("y\\[",holdout$full_individual,"\\]"))) %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="full_individual",values_to="predID_idx") %>%
  mutate_at(c('full_individual'), readr::parse_number) %>%
  mutate(full_individual=as.character(full_individual)) %>%
  left_join(holdout %>% select(full_individual, trueID_idx)) %>%
  mutate(match = (predID_idx==trueID_idx))

reduced_y_out <- MCMCchains(reduced_val_jm, params = 'y') %>%
  as_tibble() %>%
  select(matches(paste0("y\\[",holdout$reduced_individual,"\\]"))) %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="reduced_individual",values_to="predID_idx") %>%
  mutate_at(c('reduced_individual'), readr::parse_number) %>%
  mutate(reduced_individual=as.character(reduced_individual)) %>%
  left_join(holdout %>% select(reduced_individual, trueID_idx)) %>%
  mutate(match = (predID_idx==trueID_idx))

## Accuracy
# Overall (across all species)
acc_full <- full_y_out %>% 
  group_by(draw) %>%
  summarize(accuracy=mean(match))
  count(match) %>%
  filter(match==T) %>%
  pull(n) / nrow(full_y_out)

acc_reduced <- reduced_y_out %>% 
  count(match) %>%
  filter(match==T) %>%
  pull(n) / nrow(reduced_y_out)

full_y_out %>% 
  group_by(draw) %>%
  summarize(accuracy=mean(match)) %>%
  mutate(model="full") %>%
  rbind(reduced_y_out %>% 
          group_by(draw) %>%
          summarize(accuracy=mean(match)) %>%
              mutate(model="reduced")) %>%
  ggplot(aes(x=accuracy, col=model)) +
  geom_density()

# By species 
acc_full_spec <- full_y_out %>% 
  count(trueID_idx, match) %>%
  filter(match==T) %>%
  left_join(holdout %>% count(trueID_idx) %>% rename(num_inds=n)) %>%
  mutate(acc = n / (num_inds * nrow(full_y_out)/nrow(holdout)))

acc_reduced_spec <- reduced_y_out %>% 
  count(trueID_idx, match) %>%
  filter(match==T) %>%
  left_join(holdout %>% count(trueID_idx) %>% rename(num_inds=n)) %>%
  mutate(acc = n / (num_inds * nrow(full_y_out)/nrow(holdout)))

# AIS how to break this into by species by draw posterior density?
full_y_out %>% 
  group_by(draw, trueID_idx) %>%
  summarize(accuracy=mean(match)) %>%
  mutate(model="full") %>%
  rbind(reduced_y_out %>% 
          group_by(draw, trueID_idx) %>%
          summarize(accuracy=mean(match)) %>%
          mutate(model="reduced")) %>%
  ggplot(aes(x=accuracy, col=model)) +
  geom_density() +
  facet_wrap(~trueID_idx)
#AIS change from true species index to name

# Calculate classification metrics

## Recall - the number of correctly predicted spA out of the number of actual spA
# Overall (across all species)
# By species 

## Precision - the number of correctly predicted spA out of all predicted spA individuals
# Overall (across all species)
# By species 

## F1 score
# Overall (across all species)
# By species 

# Holdout log likelihood (log probabilities from the categorical distribution) 
# (could use dcat() function: https://www.rdocumentation.org/packages/mvc/versions/1.3/topics/dcat)
library(mvc) #dcat()
#dcat(, theta row for a given species)
  # Compute overall (across all species) or on a species-by-species basis
