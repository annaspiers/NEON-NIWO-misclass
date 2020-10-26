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
y_df_full <- y_df
rm(jagsinput, L, y_df)

# Load jags input for reduced model
jagsinput <- return_jags_input("reduced")
list2env(jagsinput, .GlobalEnv)
L_reduced <- L
y_df_reduced <- y_df
rm(jagsinput, L, y_df)

# To validate models, withold parataxonomist IDs from individuals that have expert IDs
perc_withld <- 0.2 #percent individuals with expert IDs that will have paratax ID withheld
# Select rows at random that will have paratax ID withheld
set.seed(49)
rows_withld <- y_df_full %>%
  rownames_to_column() %>%
  filter(!is.na(expertID_idx)) %>%
  slice_sample(prop=perc_withld) %>%
  pull(rowname) #select vector of row indices that will be withheld
# Assign parataxID_idx in selected rows to NA
y_df_val <- y_df_full %>%
  rownames_to_column() %>%
  mutate(parataxID_idx_val = ifelse(rowname %in% rows_withld, NA, parataxID_idx)) 

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
                  site = y_df_val$plotID_idx,
                  year = as.numeric(as.factor(y_df_val$col_year)),
                  k = y_df_val$expertID_idx,
                  y = y_df_val$parataxID_idx_val,
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
                     k = y_df_val %>%
                       filter(!is.na(expertID_idx)) %>%
                       pull(expertID_idx),
                     y = y_df_val %>%
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
holdout <- y_df_val %>% 
  rename(full_individual=rowname) %>%
  filter(!is.na(expertID_idx)) %>%
  rownames_to_column() %>%
  rename(reduced_individual=rowname) %>%
  filter(is.na(parataxID_idx_val)) %>%
  dplyr::select(full_individual, reduced_individual, trueID_idx=parataxID_idx, trueID=expertID,
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
  mutate(match = (predID_idx==trueID_idx), 
         model = "full")

reduced_y_out <- MCMCchains(reduced_val_jm, params = 'y') %>%
  as_tibble() %>%
  dplyr::select(matches(paste0("y\\[",holdout$reduced_individual,"\\]"))) %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="reduced_individual",values_to="predID_idx") %>%
  mutate_at(c('reduced_individual'), readr::parse_number) %>%
  mutate(reduced_individual=as.character(reduced_individual)) %>%
  left_join(holdout %>% dplyr::select(reduced_individual, trueID_idx)) %>%
  mutate(match = (predID_idx==trueID_idx), 
         model = "reduced")


# Generate confusion matrices for each posterior draw ---------------------
max_draw <- 2000
all_combos <- expand.grid(draw = 1:max(full_y_out$draw), 
                          predID_idx = 1:max(y_df_full$parataxID_idx), 
                          trueID_idx = 1:max(y_df_full$parataxID_idx, 
                                             na.rm = TRUE))  %>%
  as_tibble %>%
  filter(draw < max_draw)

full_cm <- full_y_out %>%
  filter(draw < max_draw) %>%
  count(draw, predID_idx, trueID_idx) %>%
  full_join(all_combos) %>% # fill in implicit zeros
  reshape2::acast(draw ~ trueID_idx ~ predID_idx, 
                  value.var = "n", 
                  fill = 0)

reduced_cm <- reduced_y_out %>%
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

full_metrics <- apply(full_cm, 1, get_metrics) %>%
  bind_rows(.id = "draw") %>%
  mutate(model = "full")
reduced_metrics <- apply(reduced_cm, 1, get_metrics) %>%
  bind_rows(.id = "draw") %>%
  mutate(model = "reduced")

# compare distribution of macro F1 (the average of species-specific F1 scores)
full_join(full_metrics, reduced_metrics) %>%
  group_by(draw, model) %>%
  summarize(macro_f1 = mean(f1, na.rm = TRUE)) %>%
  ggplot(aes(macro_f1, fill = model)) + 
  geom_density(alpha = .5)
  

# compare f1 scores by species for each model
# note that these plots probably aren't worth including, but still they 
# do highlight which species have the biggest differences in performance
# (again, the substantial differences show up for the common species)
full_join(full_metrics, reduced_metrics) %>%
  select(draw, parataxID_idx, f1, model) %>%
  pivot_wider(names_from = "model", values_from = "f1") %>%
  na.omit %>%
  left_join(count(y_df_full, parataxID_idx, parataxID)) %>%
  ggplot(aes(x = full, y = reduced, color = n)) + 
  geom_point(alpha = .1) + 
  facet_wrap(~forcats::fct_reorder(parataxID, reduced - full)) +
  geom_abline(linetype = "dashed") + 
  xlab("Full model F1 score") + 
  ylab("Reduced model F1 score") + 
  scale_color_viridis_c(trans = "log10", "Individuals")


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
  mutate(model="full") %>%
  rbind(reduced_y_out %>% 
          group_by(draw) %>%
          summarize(accuracy=mean(match)) %>%
              mutate(model="reduced")) %>%
  ggplot(aes(x=accuracy, col=model)) +
  geom_density()

# By species 
full_y_out %>% 
  count(trueID_idx, match) %>%
  filter(match==T) %>%
  left_join(holdout %>% count(trueID_idx) %>% rename(num_inds=n)) %>%
  mutate(acc = n / (num_inds * nrow(full_y_out)/nrow(holdout)))

reduced_y_out %>% 
  count(trueID_idx, match) %>%
  filter(match==T) %>%
  left_join(holdout %>% count(trueID_idx) %>% rename(num_inds=n)) %>%
  mutate(acc = n / (num_inds * nrow(full_y_out)/nrow(holdout)))

full_y_out %>% 
  group_by(draw, trueID_idx) %>%
  summarize(accuracy=mean(match), 
            model = 'full') %>%
  full_join(reduced_y_out %>% 
          group_by(draw, trueID_idx) %>%
          summarize(accuracy=mean(match), 
                    model = 'reduced')) %>%
  mutate(trueID = rownames[trueID_idx]) %>%
  ggplot(aes(x=accuracy, col=model)) +
  geom_density() +
  facet_wrap(~trueID)
# the full model is more accurate than the reduced model most pronouncedly for the most abundant species



# Holdout log likelihood (log probabilities from the categorical distribution) 
# (could use dcat() function: https://www.rdocumentation.org/packages/mvc/versions/1.3/topics/dcat)
theta_summ_full_val <- MCMCchains(full_val_jm, params = 'Theta') %>% 
  as_tibble() %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="index",values_to="value") %>%
  separate("index", into = c("trueID_idx", "predID_idx"), sep = ",") %>%
  mutate_at(c('trueID_idx', 'predID_idx'), readr::parse_number) 

theta_summ_reduced_val <- MCMCchains(reduced_val_jm, params = 'Theta') %>% 
  as_tibble() %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="index",values_to="value") %>%
  separate("index", into = c("trueID_idx", "predID_idx"), sep = ",") %>%
  mutate_at(c('trueID_idx', 'predID_idx'), readr::parse_number) 


full_loglik <- holdout %>%
  count(trueID_idx, predID_idx) %>%
  left_join(theta_summ_full_val) %>%
  group_by(draw) %>%
  summarize(log_lik = sum(n * log(value))) %>%
  mutate(model = "full")


reduced_loglik <- holdout %>%
  count(trueID_idx, predID_idx) %>%
  left_join(theta_summ_reduced_val) %>%
  group_by(draw) %>%
  summarize(log_lik = sum(n * log(value))) %>%
  mutate(model = "reduced")


full_join(full_loglik, reduced_loglik) %>%
  ggplot(aes(log_lik, fill = model)) + 
  geom_density() + 
  xlab("Holdout log-likelihood") + 
  ylab("Density")
