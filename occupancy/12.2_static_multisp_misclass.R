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

# Load data
# Filter for 1 year (static model)
all_paratax_df_2018 <- readRDS("occupancy/all_paratax_df.rds") %>%
    rename(para_morph_combo = scimorph_combo) %>%
    filter(col_year == "2018")
pinned_df_2018 <- readRDS("occupancy/pinned_df.rds") %>%
    filter(col_year == "2018")
expert_df_2018 <- readRDS("occupancy/expert_df.rds") %>%
    filter(col_year == "2018")

all_paratax_by_ind_2018 <- all_paratax_df_2018 %>%
  uncount(individualCount) %>%
  rownames_to_column() 
all_paratax_by_ind_2018$exp_sciname <- NA

# Generate names for parataxonomist and expert taxonomist that will be rows/columns in theta
# List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the row names of theta and M
rownames <- c(expert_df_2018 %>%   #1) the unique expert taxonomist ID's
                distinct(scientificName) %>%
                pull(scientificName),
              all_paratax_df_2018 %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                distinct(scientificName) %>%
                filter(!scientificName %in% unique(expert_df_2018$scientificName)) %>% 
                pull(scientificName) ) 
rownames <- sort(rownames)
#AIS rather than using the parataxonomist scientific IDs that the expert taxonomist didn't use, 
# wouldn't we want to instead use paratazonomist scientificname-morphospecies combo that the expert didn't use?

# List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
# parataxonomist didn't use. 
# These identifications will be the column names of theta and M
extra_colnames <- c(all_paratax_df_2018 %>%   #1) the unique parataxonomist and morphospecies ID's
                      distinct(para_morph_combo) %>%
                      pull(para_morph_combo),
                    expert_df_2018 %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                      distinct(scientificName) %>%
                      filter(!scientificName %in% unique(all_paratax_df_2018$para_morph_combo)) %>% 
                      pull(scientificName) ) 
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
subsamps <- pinned_df_2018 %>%
    left_join(expert_df_2018 %>% 
                  dplyr::select(individualID, exp_sciname = scientificName)) %>%
    dplyr::select(subsampleID, exp_sciname) %>%
    filter(!is.na(exp_sciname)) %>%
    distinct(subsampleID) %>% pull(subsampleID)
# para_new will replace the all_paratax_df since dplyr doesn't play nice in for-loops
# (I couldn't assign to a filtered object)
# initialize para_new with subsamples that the expert ID never looked at
para_new <- all_paratax_by_ind_2018 %>% 
    filter(!(subsampleID %in% subsamps)) 
for (id in subsamps) {
    temp_exp_df <- pinned_df_2018 %>%
        left_join(expert_df_2018 %>% 
                      dplyr::select(individualID, exp_sciname = scientificName)) %>%
        dplyr::select(subsampleID, exp_sciname) %>%
        filter(!is.na(exp_sciname),
               subsampleID == id)
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
assertthat::assert_that(sum(!is.na(para_new$exp_sciname)) == nrow(expert_df_2018))
# Choose a subsample that has fewer expert IDs than sorting IDs: ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=
pinned_df_2018 %>%
  left_join(expert_df_2018 %>% 
              dplyr::select(individualID, exp_sciname = scientificName)) %>%
  dplyr::select(subsampleID, exp_sciname) %>%
  filter(!is.na(exp_sciname), subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
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
              dplyr::select(plotID, plotID_idx))

# Define L
L <- reshape2::acast(para_new, plotID ~ collectDate)

# Define alpha
alpha <- matrix(1, nrow = length(rownames), ncol = length(colnames))
diag(alpha) <- 10

## "Ground truth" data
# We have a subset of the data with known species IDs from expert identification
# We partly observe k, the expertID column in y_df
# We partly observe z
z.dat <- array(NA, dim = c(length(unique(para_new$plotID)),
                           length(rownames)),
               dimnames = list(sort(unique(para_new$plotID)), #plots
                               rownames)) #expert IDs
# Grab values from casted z.dat array and fill in values in final z.dat array
z.dat_cast <- expert_df_2018 %>% 
  mutate(occ = 1) %>%
  reshape2::acast(plotID ~ scientificName,
                  fill=-999, drop=F, value.var = "occ")
z.dat_cast[z.dat_cast == -999] <- NA
z.dat_cast[z.dat_cast > 0] <- 1

assertthat::assert_that(sum(z.dat_cast, na.rm=T) == nrow(expert_df_2018 %>% distinct(plotID, scientificName)))

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
               site = y_df$plotID_idx, #needs to be numeric
               # if the individual was labeled by the expert, true ID is known
               k = y_df$expertID,
               # for all individuals, we get paratxonomist IDs
               y = y_df$parataxID,
               z = z.dat)

#ji <- function(){ list(z = z.init) }

# jm <- jags.parallel(
#     data = jags_d, 
#     inits = list(z = z.init), #ji,
#     parameters.to.save = c("psi", "p", "phi", "gamma", "lambda", 
#                            "Theta"), 
#     model.file = "occupancy/12.2_static_multisp_misclass_JAGS.txt", 
#     n.chains = 6, 
#     n.iter = 100000, 
#     n.thin = 1,
#     DIC = FALSE)
JAGSinits <- function(){ list(z = z.init) }
nc=4
cl <- makeCluster(nc)
jm <- jags.parfit(cl = cl,
                   data = jags_d,
                   params = c("psi", "p", "phi", "gamma", "lambda", "Theta"),
                   model = "occupancy/12.2_static_multisp_misclass_JAGS.txt",
                   inits = JAGSinits,
                   n.chains = nc,
                   n.adapt = 1000,
                   n.update = 1000,
                   thin = 1,
                   n.iter = 10000)

jm_summ <- MCMCsummary(jm)

# Check convergence
par(mfrow=c(1,1))
hist(jm_summ$Rhat)
#yusssss
