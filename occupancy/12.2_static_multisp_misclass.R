# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a static occupancy 
# (single season) model estimating misclassification probabilities with 
# NEON data.

library(reshape2) 
library(tidyverse) 
library(gtools) 
library(R2jags)
library(dclone)
library(MCMCvis)
library(ggplot2)

# Load data
# Filter for 1 year (static model)
all_paratax_by_ind_2018 <- readRDS("occupancy/all_paratax_df.rds") %>%
  rename(para_morph_combo = scimorph_combo) %>%
  filter(col_year == "2018") %>%
  uncount(individualCount) %>%
  rownames_to_column() 
all_paratax_by_ind_2018$exp_sciname <- NA
pinned_df_2018 <- readRDS("occupancy/pinned_df.rds") %>%
    filter(col_year == "2018")
expert_df_2018 <- readRDS("occupancy/expert_df.rds") %>%
    filter(col_year == "2018")
expert_pinned_df_2018 <- expert_df_2018 %>%
  rename(exp_sciname = scientificName) %>%
  left_join(pinned_df_2018 %>% dplyr::select(individualID, subsampleID, para_morph_combo = scimorph_combo))

# Filter datasets to abundant parataxonomist species/morphospecies IDs
common_morpho <- all_paratax_by_ind_2018 %>%
  group_by(para_morph_combo) %>%
  summarize(total = n()) %>%
  arrange(-total) %>%
  filter(total > 10)

all_paratax_by_ind_2018 <- all_paratax_by_ind_2018 %>%
  filter(para_morph_combo %in% common_morpho$para_morph_combo)
expert_pinned_df_2018 <- expert_pinned_df_2018 %>%
  filter(subsampleID %in% all_paratax_by_ind_2018$subsampleID)

rm(pinned_df_2018, expert_df_2018)

# Generate names for parataxonomist and expert taxonomist that will be rows/columns in theta
# List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the row names of theta and M
rownames <- c(expert_pinned_df_2018 %>%   #1) the unique expert taxonomist ID's
                distinct(exp_sciname) %>%
                pull(exp_sciname),
              all_paratax_by_ind_2018 %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                distinct(scientificName) %>%
                filter(!scientificName %in% unique(expert_pinned_df_2018$scientificName)) %>% 
                pull(scientificName) ) 
rownames <- sort(rownames)
#AIS rather than using the parataxonomist scientific IDs that the expert taxonomist didn't use, 
# wouldn't we want to instead use paratazonomist scientificname-morphospecies combo that the expert didn't use?

# List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
# parataxonomist didn't use. 
# These identifications will be the column names of theta and M
extra_colnames <- c(all_paratax_by_ind_2018 %>%   #1) the unique parataxonomist and morphospecies ID's
                      distinct(para_morph_combo) %>%
                      pull(para_morph_combo),
                    expert_pinned_df_2018 %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                      distinct(exp_sciname) %>%
                      filter(!exp_sciname %in% unique(all_paratax_by_ind_2018$para_morph_combo)) %>% 
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
subsamps <- expert_pinned_df_2018 %>%
    distinct(subsampleID) %>% pull(subsampleID)
# para_new will replace the all_paratax_df since dplyr doesn't play nice in for-loops
# (I couldn't assign to a filtered object)
# initialize para_new with subsamples that the expert ID never looked at
para_new <- all_paratax_by_ind_2018 %>% 
    filter(!(subsampleID %in% subsamps)) 
for (id in subsamps) {
    temp_exp_df <- expert_pinned_df_2018 %>%
        dplyr::select(subsampleID, exp_sciname) %>%
        filter(subsampleID == id)
    temp_para_df <- all_paratax_by_ind_2018 %>% 
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

# Sanity check
assertthat::assert_that(sum(!is.na(para_new$exp_sciname)) == nrow(expert_pinned_df_2018))
# Choose a subsample that has fewer expert IDs than sorting IDs: ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=
expert_pinned_df_2018 %>%
  dplyr::select(subsampleID, exp_sciname) %>%
  filter(subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
para_new %>% filter(subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
# We should see 6 out of 7 individuals in para_new have an assigned exp_sciname

## Imperfect species classifications 
# Probability vector `y` for with a record for each detection. 
# Noisy classifier with skill that vary by species.
y_df <- para_new %>%
    dplyr::select(plotID, collectDate, parataxID = para_morph_combo, 
                  expertID = exp_sciname)
y_df <- y_df %>%
  left_join(y_df %>%
              group_by(plotID) %>%
              summarize(n=n()) %>%
              mutate(plotID_idx = 1:n()) %>%
              dplyr::select(plotID, plotID_idx)) %>%
  mutate(parataxID_idx = match(y_df$parataxID, colnames),
         expertID_idx = match(y_df$expertID, rownames))

# Define L
L <- reshape2::acast(para_new, plotID ~ collectDate)

# Define alpha
alpha <- matrix(2, nrow = length(rownames), ncol = length(colnames))
diag(alpha) <- 80

## "Ground truth" data
# We have a subset of the data with known species IDs from expert identification
# We partly observe k, the expertID column in y_df
# We partly observe z
z.dat <- array(NA, dim = c(length(unique(para_new$plotID)),
                           length(rownames)),
               dimnames = list(sort(unique(para_new$plotID)), #plots
                               rownames)) #expert IDs
# Grab values from casted z.dat array and fill in values in final z.dat array
z.dat_cast <- expert_pinned_df_2018 %>% 
  mutate(occ = 1) %>%
  reshape2::acast(plotID ~ exp_sciname,
                  fill=-999, drop=F, value.var = "occ")
z.dat_cast[z.dat_cast == -999] <- NA
z.dat_cast[z.dat_cast > 0] <- 1

assertthat::assert_that(sum(z.dat_cast, na.rm=T) == nrow(expert_pinned_df_2018 %>% distinct(plotID, exp_sciname)))

for(plot in dimnames(z.dat_cast)[[1]]) {
  for(spec in dimnames(z.dat_cast)[[2]]) {
    z.dat[plot,spec] <- z.dat_cast[plot,spec]
  }
}
rm(z.dat_cast,plot,spec)

# Initialize Z
z.init <- z.dat
for (i in 1:dim(z.init)[1]) {
    z.init[i,] <- sample(c(0,1), replace=TRUE, size=dim(z.init)[2])
}
# initialize known values as NA, otherwise model will throw error
z.init[z.dat == 1] <- NA

# Check that where L>0 for a species, z.init>0 for that species/site/year combo
for (i in 1:dim(z.init)[1]) {
  for (k in 1:dim(z.init)[2]) {
      if (sum(L[i,], na.rm = TRUE) > 0 ) {
        ifelse(z.init[i,k] == 0, 1, z.init[i,k])
    }
  }
}
rm(i,k)


## Model fitting
# Fit the model and check convergence.
jags_d <- list(nsite = length(unique(para_new$plotID)),
               K_para = length(colnames),
               K_exp = length(rownames), 
               noccasion = length(unique(para_new$collectDate)), 
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = as.numeric(y_df$plotID_idx), #needs to be numeric
               # if the individual was labeled by the expert, true ID is known
               k = as.numeric(y_df$expertID_idx),
               # for all individuals, we get paratxonomist IDs
               y = as.numeric(y_df$parataxID_idx),
               z = z.dat,
               R = diag(rep(1, 3)))
JAGSinits <- function(){ list(z = z.init) }
nc=4
cl <- makeCluster(nc)
jm <- jags.parfit(cl = cl,
                   data = jags_d,
                   params = c("logit_psi", "logit_p", "log_lambda", "Theta",
                              "eps_site", "eps_spec"),
                   model = "occupancy/12.2_static_multisp_misclass_JAGS.txt",
                   inits = JAGSinits,
                   n.chains = nc,
                   n.adapt = 1000,
                   n.update = 1000,
                   thin = 1,
                   n.iter = 10000)
#took 22ish min
saveRDS(jm, "occupancy/script12.2_jags_jm.rds")
jm_summ <- MCMCsummary(jm)
saveRDS(jm_summ, "occupancy/script12.2_jm_mcmcsumm.rds")

jm <- readRDS("occupancy/script12.2_jags_jm.rds")
jm_summ <- readRDS("occupancy/script12.2_jm_mcmcsumm.rds")

# Check convergence
hist(jm_summ$Rhat)
range(jm_summ$Rhat, na.rm=TRUE)

# Look at high Rhat values
jm_summ <- rownames_to_column(jm_summ)
jm_summ %>% filter(Rhat > 1.1)

MCMCtrace(jm, params = paste0('eps_spec\\[1,1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('eps_spec\\[1,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('eps_spec\\[6,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('eps_spec\\[14,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('logit_p\\[1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('logit_p\\[6\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('logit_p\\[14\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('logit_psi\\[1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 

# theta
theta_summ <- MCMCsummary(jm, params = 'Theta', round=2)
theta_summ <- rownames_to_column(theta_summ)
hist(theta_summ$Rhat)
range(theta_summ$Rhat, na.rm=TRUE)
theta_summ %>% filter(Rhat > 2)

MCMCtrace(jm, params = paste0('Theta\\[2,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('Theta\\[24,2\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
MCMCtrace(jm, params = paste0('Theta\\[24,22\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
theta_prior <- MCMCpack::rdirichlet(1000, c(100,rep(1.5,dim(alpha)[2]-1) ))
plot(density(theta_prior[, 1]))
# left_join(expert_df %>% select(individualID, expert=scientificName), pinned_df %>% select(individualID, para = scientificName)) %>% filter(scientificName==rownames[10])
# MCMCsummary(jm, params='theta\\[36,[0-9]+\\]', ISB=F)

theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_mean = theta_summ$mean)
#make para_morph a factor to force plotting in order
theta_df$para_morph = factor(theta_df$para_morph, levels=colnames) 

# Plot heatmap of theta values
ggplot(theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_mean)) + 
  geom_tile() +
  #scale_fill_gradient(trans="log10",low="darkblue", high="white") +
  scale_fill_gradient(limits=c(0,1),low="darkblue", high="white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))
#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html

# AIS hmm, why is expert ID Amara quenseli matching with a high theta mean with 
# paratax ID Carabus taedatus? HTere's no crossover...
para_new %>% filter(exp_sciname=="Amara quenseli")
para_new %>% filter(para_morph_combo=="Carabus taedatus" & !is.na(exp_sciname))

# AIS why no high probabilities for the morphospecies?
para_new %>% filter(para_morph_combo==colnames[26])
para_new %>% filter(para_morph_combo==colnames[27])
