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
library(mltest) #ml_test
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
all_paratax_by_ind <- readRDS("data/all_paratax_df.rds") %>%
  rename(para_morph_combo = scimorph_combo)  %>%
  uncount(individualCount) %>%
  rownames_to_column() 
all_paratax_by_ind$exp_sciname <- NA
expert_df <- readRDS("data/expert_df.rds") 
pinned_df <- readRDS("data/pinned_df.rds") 
expert_pinned_df <- expert_df %>%
  rename(exp_sciname = scientificName) %>%
  left_join(pinned_df %>% 
              dplyr::select(individualID, subsampleID, 
                            para_morph_combo = scimorph_combo))
rm(expert_df, pinned_df)

# List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the row names of theta and M
rownames <- c(expert_pinned_df %>%   #1) the unique expert taxonomist ID's
                distinct(exp_sciname) %>%
                pull(exp_sciname),
              all_paratax_by_ind %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                distinct(scientificName) %>%
                filter(!scientificName %in% unique(expert_pinned_df$exp_sciname)) %>% 
                pull(scientificName) ) 
rownames <- sort(rownames)

# List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
# parataxonomist didn't use. 
# These identifications will be the column names of theta and M
extra_colnames <- c(all_paratax_by_ind %>%   #1) the unique parataxonomist and morphospecies ID's
                      distinct(para_morph_combo) %>%
                      pull(para_morph_combo),
                    expert_pinned_df %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                      distinct(exp_sciname) %>%
                      filter(!exp_sciname %in% unique(all_paratax_by_ind$para_morph_combo)) %>% 
                      pull(exp_sciname) ) 
extra_colnames <- sort(extra_colnames)
extra_indices <- which(!extra_colnames %in% rownames)
colnames <- c(rownames, extra_colnames[extra_indices])
rm(extra_colnames, extra_indices)

assertthat::assert_that(all(rownames == colnames[1:length(rownames)]))

# Join parataxonomist and expert ID tables ---------------------------
# Ex: If parataxonomist counts 5 animals in one subsample & expert IDs two
# individuals, assign known species IDs to 2 of the 5 rows at random (this is
# valid because the 5 individuals are exchangeable)

# initialize counter for for-loop
subsamps <- expert_pinned_df %>%
  distinct(subsampleID) %>% pull(subsampleID)
# para_new will replace the all_paratax_df since dplyr doesn't play nice in for-loops
# initialize para_new with subsamples that the expert ID never looked at
para_new <- all_paratax_by_ind %>% 
  filter(!(subsampleID %in% subsamps)) 
for (id in subsamps) {
  temp_exp_df <- expert_pinned_df %>%
    dplyr::select(subsampleID, exp_sciname) %>%
    filter(!is.na(exp_sciname), subsampleID == id)
  temp_para_df <- all_paratax_by_ind %>% 
    filter(subsampleID == id)
  for (row in 1:nrow(temp_exp_df)) {
    # check whether exp_sciname is empty in first row of paratax-df
    while (!is.na(temp_para_df$exp_sciname[row])) {
      row = row + 1
    } 
    # assign exp_sciname from pinned-expert-combo-df to corresponding row of paratax-df
    temp_para_df$exp_sciname[row] <- temp_exp_df$exp_sciname[row]
  }
    para_new <- rbind(para_new, temp_para_df)
}

## Sanity check
assertthat::assert_that(sum(!is.na(para_new$exp_sciname)) == nrow(expert_pinned_df))
# Choose a subsample that has fewer expert IDs than sorting IDs: ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=
expert_pinned_df %>%
  filter(!is.na(exp_sciname), subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
para_new %>% filter(subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
# We should see 6 out of 7 individuals in para_new have an assigned exp_sciname

# Plot discrepancies for samples identified to the species level
para_new %>%
  filter(!is.na(exp_sciname)) %>% #taxonRank == "species",
  count(para_morph_combo, exp_sciname) %>%
  arrange(n) %>%
  mutate(discrepancy = para_morph_combo != exp_sciname) %>%
  ggplot(aes(x = para_morph_combo, 
             y = exp_sciname, 
             color = discrepancy)) + 
  geom_point(aes(size = n)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Species ID from parataxonomist") + 
  ylab("Species ID from expert taxonomist") + 
  scale_color_manual(values = c("black", "red"))

## Imperfect species classifications 
# Probability vector `y` for with a record for each detection. 
# Noisy classifier with skill that vary by species.
y_df <- para_new %>%
  dplyr::select(plotID, collectDate, parataxID = para_morph_combo, 
                expertID = exp_sciname, col_year)
y_df <- y_df %>%
  left_join(y_df %>%
              group_by(plotID) %>%
              summarize(n=n()) %>%
              mutate(plotID_idx = 1:n()) %>%
              dplyr::select(plotID, plotID_idx)) %>%
  mutate(parataxID_idx = match(y_df$parataxID, colnames),
         expertID_idx = match(y_df$expertID, rownames))

# Define L
L_full <- reshape2::acast(para_new, plotID ~ collectDate ~ col_year)
L_reduced <- reshape2::acast(subset(para_new,!is.na(exp_sciname)), plotID ~ collectDate ~ col_year)
 
# Define alpha
alpha <- matrix(2, nrow = length(rownames), ncol = length(colnames))
diag(alpha) <- 200

## "Ground truth" data
# We have a subset of the data with known species IDs from expert identification
# We partly observe k, the expertID column in y_df
# We partly observe z
z.dat <- array(NA, dim = c(length(unique(para_new$plotID)),
                           length(rownames),
                           length(unique(all_paratax_by_ind$col_year))),
               dimnames = list(sort(unique(para_new$plotID)), #plots
                               rownames, #expert IDs
                               unique(all_paratax_by_ind$col_year)))
# Grab values from casted z.dat array and fill in values in final z.dat array
z.dat_cast <- expert_pinned_df %>% 
  mutate(occ = 1) %>%
  reshape2::acast(plotID ~ exp_sciname ~ col_year,
                  fill=-999, drop=F, value.var = "occ")
z.dat_cast[z.dat_cast == -999] <- NA
z.dat_cast[z.dat_cast > 0] <- 1

assertthat::assert_that(sum(z.dat_cast, na.rm=T) == 
                          nrow(expert_pinned_df %>% distinct(plotID, exp_sciname, col_year)))

for(plot in dimnames(z.dat_cast)[[1]]) {
  for(spec in dimnames(z.dat_cast)[[2]]) {
    for(year in dimnames(z.dat_cast)[[3]]) {
      z.dat[plot,spec,year] <- z.dat_cast[plot,spec,year]
    }
  }
}
rm(z.dat_cast,plot,spec)

# Initialize Z
z.init <- z.dat
for (i in 1:dim(z.init)[1]) {
  for (t in 1:dim(z.init)[3]) {
    z.init[i,,t] <- sample(c(0,1), replace=TRUE, size=dim(z.init)[2])
  }
}
# initialize known values as NA, otherwise model will throw error
z.init[z.dat == 1] <- NA

# Check that where L>0 for a species, z.init>0 for that species/site/year combo
for (i in 1:dim(z.init)[1]) {
  for (k in 1:dim(z.init)[2]) {
    for (t in 1:dim(z.init)[3]) {
      if (sum(L_full[i,,t], na.rm = TRUE) > 0 ) {
        ifelse(z.init[i,k,t] == 0, 1, z.init[i,k,t])
      }
      if (sum(L_reduced[i,,t], na.rm = TRUE) > 0 ) {
        ifelse(z.init[i,k,t] == 0, 1, z.init[i,k,t])
      }
    }
  }
}
rm(i,k,t)

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
