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
library(mvc) #dcat()

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
  mutate(match = (predID_idx==trueID_idx))

reduced_y_out <- MCMCchains(reduced_val_jm, params = 'y') %>%
  as_tibble() %>%
  dplyr::select(matches(paste0("y\\[",holdout$reduced_individual,"\\]"))) %>%
  mutate(draw=1:n()) %>%
  pivot_longer(cols=-draw,names_to="reduced_individual",values_to="predID_idx") %>%
  mutate_at(c('reduced_individual'), readr::parse_number) %>%
  mutate(reduced_individual=as.character(reduced_individual)) %>%
  left_join(holdout %>% dplyr::select(reduced_individual, trueID_idx)) %>%
  mutate(match = (predID_idx==trueID_idx))

 # Calculate validation metrics by species
val_array <- array(NA, dim=c(length(unique(full_y_out$trueID_idx)), 4, length(unique(full_y_out$draw))),
                      dimnames=list(full_y_out %>% distinct(trueID_idx) %>% arrange(trueID_idx) %>% pull(trueID_idx),
                                    c("accuracy","recall","precision","F1"),
                                    unique(full_y_out$draw)))

# Terribly inefficient but gets the job done
for (i in 1:length(unique(full_y_out$draw))) {
  tempdf <- full_y_out %>% filter(draw == i)
  cm <- reshape2::dcast(full_y_out,  trueID_idx ~ predID_idx) %>%
    column_to_rownames(var="trueID_idx")
  
  for (true in dimnames(val_array)[[1]]) {
    val_array[true,"accuracy",i] <- tempdf %>% 
      group_by(trueID_idx) %>% 
      summarise(accuracy = mean(match)) %>% 
      filter(trueID_idx==as.numeric(true)) %>%
      pull(accuracy)
    val_array[true,"recall",i] <- cm[true,as.numeric(true)] / sum(cm[true,])
    val_array[true,"precision",i] <- cm[true,as.numeric(true)] / sum(cm[,as.numeric(true)])
    val_array[true,"F1",i] <- 2*val_array[true,"recall",i]*val_array[true,"precision",i]/
      (val_array[true,"recall",i] + val_array[true,"precision",i])
  }
}
validation_metrics <- function(df) { #AIS rewrite to handle one draw_ID df at a time
    cm <- reshape2::dcast(df,  trueID_idx ~ predID_idx) %>%
      column_to_rownames(var="trueID_idx")
    
      acc <- df %>% 
        group_by(trueID_idx) %>% 
        summarise(accuracy = mean(match)) %>% 
        pull(accuracy)
      recall <- cm[as.character(df$trueID_idx),df$trueID_idx] / sum(cm[as.character(df$trueID_idx),])
      precision <- cm[as.character(df$trueID_idx),df$trueID_idx] / sum(cm[,df$trueID_idx])
      F1 <- 2*recall*precision/ (recall + precision)
      val_list <- list(acc, recall, precision, F1)
      return(val_list)
}
full_y_out %>% 
  #mutate(draw_ID = paste0(draw,"_",trueID_idx)) %>%
  mutate(fdraw = as.factor(draw)) %>%
  split(.$fdraw) %>%
  lapply(validation_metrics())

full_y_out %>%
  group_by(draw) %>%
  summarize(recall = sum(match)/(n()))
# Accuracy
# Recall - the number of correctly predicted spA out of the number of spA
# Precision - the number of correctly predicted spA out of all predicted spA individuals
# F1 score


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
            recall=sum(match)/sum(sum(match)+sum(full_y_out %>% filter())),
            prec=sum(match)/sum(sum(match)+sum(predID_idx))) %>%
  mutate(model="full") %>%
  rbind(reduced_y_out %>% 
          group_by(draw, trueID_idx) %>%
          summarize(accuracy=mean(match),
                    recall=sum(match)/sum(sum(match)+sum(trueID_idx)),
                    prec=sum(match)/sum(sum(match)+sum(predID_idx))) %>%
          mutate(model="reduced")) %>%
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
holdout %>%
  slice(rep(1:n(), each = 4000)) %>%
  mutate(draw = rep(1:4000, nrow(holdout))) %>%
  sample_frac(0.00001) %>%
  rowwise() %>%
  #AIS left_join
  #AIS try rowwise
  mutate(loglik = (log(theta_summ_full_val %>% 
           dplyr::filter(draw==draw & trueID_idx==trueID_idx & predID_idx==predID_idx) %>% #AIS this is returng 16,000
           pull(value))) %>%
  group_by(draw) %>%
  summarize(hldt_loglik = sum())
#AIS visualize       

#dcat(, theta row for a given species)
  # Compute overall (across all species) or on a species-by-species basis
