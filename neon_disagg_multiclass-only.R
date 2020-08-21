# We're revisiting the framing of this model and considering Andy Royal's 
# disaggregated misclassification model. Here we create a dynamic occupancy 
# model estimating misclassification probabilities with NEON Niwot Ridge carabid 
# data 2015-2018

# library(reshape2) 
library(tidyr) #uncount()
library(tibble) #rownames_to_column()
library(dclone) #alternative: R2jags::jags.parallel
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(boot) #inv.logit()
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
# Load data
all_paratax_by_ind <- readRDS("data/all_paratax_df.rds") %>%
  rename(para_morph_combo = scimorph_combo)  %>%
  uncount(individualCount) %>%
  rownames_to_column() 
all_paratax_by_ind$exp_sciname <- NA
expert_df <- readRDS("data/expert_df.rds") 
pinned_df <- readRDS("data/pinned_df.rds") 
expert_pinned_df <- expert_df %>%
  rename(exp_sciname = scientificName) %>%
  left_join(pinned_df %>% dplyr::select(individualID, subsampleID, para_morph_combo = scimorph_combo))

# subsample aRw+8bWvW/wOIENap0etmXiCq6V1RZPfCkRcy1w8BLk= had 
# 6 individuals sorted, but 8 pinned and expertly identified. 
# This is the only subsample with less sorted individuals than expertly ID'ed
# As a workaround, we'll remove two expertly identified inidividuals
removed_inds <- expert_pinned_df %>%
  filter(subsampleID == "aRw+8bWvW/wOIENap0etmXiCq6V1RZPfCkRcy1w8BLk=") %>%
  pull(individualID) %>%
  tail(2)
expert_pinned_df <- expert_pinned_df %>%
  filter(!(individualID %in% removed_inds))

rm(expert_df, pinned_df, common_morpho, removed_inds)

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
                    expert_pinned_df %>%   #2) the expert  IDs that the parataxonomist didn't use. 
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
# (I couldn't assign to a filtered object)
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
# Choose a subsample that has fewer expert IDs than sorting IDs: 
  #ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=
expert_pinned_df %>%
  filter(!is.na(exp_sciname), subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
para_new %>% filter(subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
# We should see 6 jm of 7 individuals in para_new have an assigned exp_sciname

# Keep only expert IDed individuals for reduced model
para_new <- para_new %>%
  filter(!is.na(exp_sciname))

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
L <- reshape2::acast(para_new, plotID ~ collectDate ~ col_year)

# Define alpha
alpha <- matrix(2, nrow = length(rownames), ncol = length(colnames))
diag(alpha) <- 80
# visualize alpha
test_alpha <- MCMCpack::rdirichlet(1000, c(alpha[1,1], rep(alpha[1,2], length(colnames)-1)))
hist(test_alpha[,-1], xlim=c(0,1),col="red")
hist(test_alpha[,1],  add=T)

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
jags_d <- list(K_exp = dim(alpha)[1], 
               K_para = dim(alpha)[2],
               alpha = alpha,
               Ltot = sum(L), 
               # if the individual was labeled by the expert, true ID is known
               k = y_df$expertID_idx,
               # for all individuals, we get paratxonomist IDs
               y = y_df$parataxID_idx)
JAGSinits <- function(){ list(z = z.init) }
nc <- 4
ni <- 20000
cl <- makeCluster(nc)
jm <- jags.parfit(cl = cl,
                  data = jags_d,
                  params = c("Theta"),
                  model = "neon_disagg_multiclass-only_JAGS.txt",
                  n.chains = nc,
                  n.adapt = 1000,
                  n.update = 1000,
                  thin = ni/1000,
                  n.iter = ni) 

jm_summ <- MCMCsummary(jm)

saveRDS(jm, "output/red-jm-allspec.rds")
saveRDS(jm_summ, "output/red-jm_mcmcsumm-allspec.rds")
# saveRDS(jm, "output/red-jm-over5spec.rds")
# saveRDS(jm_summ, "output/red-jm_mcmcsumm-over5spec.rds")

# View JAGS output --------------------------------------------------------

jm <- readRDS("output/jm-allspec.rds")
jm_summ <- readRDS("output/jm_mcmcsumm-allspec.rds")
# Did model converge?
hist(jm_summ$Rhat, breaks=40)
range(jm_summ$Rhat, na.rm=TRUE)
jm_summ <- rownames_to_column(jm_summ)
jm_summ %>% filter(Rhat > 1.1)

### Look at raw numbers
# THETA
theta_summ <- MCMCsummary(jm, params = 'Theta', round=2)
theta_summ <- rownames_to_column(theta_summ)
saveRDS(theta_summ, "output/theta_summ-allspec.rds")
#theta_summ <- readRDS("output/theta_summ-allspec.rds")
hist(theta_summ$Rhat)
hist(theta_summ$mean)
range(theta_summ$Rhat, na.rm=TRUE)

# Visualize theta
theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_mean = theta_summ$mean,
         theta_median = theta_summ$"50%")
#make para_morph a factor to force plotting in order
theta_df$para_morph = factor(theta_df$para_morph, levels=colnames) 
# Plot heatmap of theta values
theta_med <- ggplot(theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_median)) + 
  geom_tile() +
  scale_fill_gradient(limits=c(0,1),low="darkblue", high="white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))
#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html
png("figures/theta_conf_matrix_allspec.png")
print(theta_med)
dev.off()

# Visualize predictions of species unobserved by expert taxonomist
par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
abline(0,1)
points(x=jm$mean$theta[nspec,],y=theta[nspec,],col="red")
