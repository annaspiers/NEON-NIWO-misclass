# We're having trouble getting script 15 to fit: the full dynamic occupancy
# model with misclassification. We have an idea that estimating alpha values
# other than the diagonal of the square component would be helpful in getting
# the full occupancy model to converge. Here we try fitting just the
# misclassificaiton model on its own in two cases:
  # Case 1: estimate the weight along the diagonal of the square portion of
    # alpha (alpha_diag) AND the weights in the morphospecies columns (alpha_morph)
  # Case 2: estimate the weight along the diagonal of the square portion of
    # alpha (alpha_diag),the weights in the morphospecies columns (alpha_morph),
    # AND the weights in the non-diagonal portion of the square matrix

library(dclone) #alternative: R2jags::jags.paralle
library(jagsUI)
library(MCMCpack)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyselect)
library(tibble)
# library(mailR) #try troubleshooting dependencies later
options(digits = 3)

# Load NEON Niwot Ridge carabid data 
# Each row is a beetle identified by a parataxonomist
ind_dat <- read.csv("occupancy/model_df_by_individual_beetle.csv", header = TRUE) %>%
    mutate(collectDate = as.Date(collectDate)) %>%
    filter(collectDate <= as.Date("2018-12-31")) %>%
    dplyr::select(-c(X))

# Each row is a sample. One sample is a para_morph at a trap on a collection day   
sample_dat  <- read.csv("occupancy/model_df_by_species_in_sample_2015-2018.csv", 
                        header = TRUE) %>%
    dplyr::select(-c(X))

# Misclassification model ------------------------------------------------- 
# The parameters below assume both a parataxonomist's and an expert taxonomist's
# classification is available

# List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
# expert taxonomist didn't use. 
# These identifications will be the row names of theta and M
rownames <- c(ind_dat %>%   #1) the unique expert taxonomist ID's
                    select(expert_sciname) %>%
                    filter(!is.na(expert_sciname)) %>%
                    mutate_if(is.factor, as.character) %>%
                    distinct() %>%
                    pull(expert_sciname),
              ind_dat %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                select(para_sciname) %>%
                filter(!para_sciname %in% unique(ind_dat$expert_sciname)) %>% 
                mutate_if(is.factor, as.character) %>%
                distinct() %>%
                pull(para_sciname) ) 
rownames <- sort(rownames)

# List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
# parataxonomist didn't use. 
# These identifications will be the column names of theta and M
extra_colnames <- c(ind_dat %>%   #1) the unique parataxonomist and morphospecies ID's
                       select(para_morph) %>%
                       mutate_if(is.factor, as.character) %>%
                       distinct() %>%
                       pull(para_morph),
                     ind_dat %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                       select(expert_sciname) %>%
                       filter(!is.na(expert_sciname),
                         !expert_sciname %in% unique(ind_dat$para_morph)) %>% 
                       mutate_if(is.factor, as.character) %>%
                       distinct() %>%
                       pull(expert_sciname) ) 
extra_colnames <- sort(extra_colnames)
extra_indices <- which(!extra_colnames %in% rownames)
colnames <- c(rownames, extra_colnames[extra_indices])
rm(extra_colnames, extra_indices)


# n[k]: expert taxonomist's count of individuals of species k 
n <- left_join(data.frame(rownames) %>%
                 rename("expert_sciname" = rownames) %>%
                 mutate_if(is.factor, as.character) ,
               ind_dat %>%
                 filter(!is.na(expert_sciname)) %>%
                 group_by(expert_sciname) %>%
                 summarize(n=n()) %>%
                 mutate_if(is.factor, as.character) ) %>% 
    pull(n) 
n[is.na(n)] <- 0

# Theta: misclassification probability matrix [KxK]
# Theta has two parts. 
# Part 1: a square matrix where each row/column represents a 'valid' taxonomic
# group 
# Part 2: matrix appended as new columns on Theta that are “invalid” taxonomic
# groups (all of the morphospecies)


# M_k: M[k,k'] is the number of individuals from species k (according to expert)
# that were identified as species k' by parataxonomist [KxK]

# create an empty df with the columns and rownames that I want
M <- matrix(0, nrow = length(rownames), ncol = length(colnames))
M <- data.frame(M, row.names = rownames)
colnames(M) <- colnames

# Grab values from casted M matrix and fill in values in final M matrix
M_cast <- ind_dat %>%
  filter(is.na(expert_sciname)==F) %>% 
  reshape2::acast(expert_sciname ~ para_morph)

for(col in colnames(M_cast)) {
  for(row in row.names(M_cast)) {
    M[row,col] <- M_cast[row,col]
  }
}
rm(col, row, M_cast)

# Convert M into matrix
M <- data.matrix(M)


# JAGS model --------------------------------------------------------------


# Run model in JAGS. 
JAGSdata <- list(nspec_exp = dim(M)[1], #species ID'd by expert
                 nspec_para = dim(M)[2], #paratax ID's and morphospecies
                 n = n,
                 M = M) #bundle data
nc <- 4 #MCMC chains
ni <- 4000 #MCMC iterations
nt <- 1  #MCMC thin

# JAGS model for Case 1
# Case 1: estimate the weight along the diagonal of the square portion of
# alpha (alpha_diag) AND the weights in the morphospecies columns (alpha_morph)
cl <- makeCluster(4)
out1 <- jags.parfit(cl = cl, 
                   data = JAGSdata,
                   params = c("theta","alpha_diag", "alpha_morph") ,
                   model = "occupancy/16a_neon_just_multisp_misclass_JAGS.txt", 
                  # inits = JAGSinits,
                   n.chains = nc,
                   n.adapt = 2000,
                   n.update = 2000,
                   thin = nt,
                   n.iter = ni) 
saveRDS(out1, "occupancy/script16_jags_out1.rds")

# How well does the model estimate specified parameters?
library(MCMCvis)
out1_mcmcsumm <- MCMCsummary(out1)
saveRDS(out1_mcmcsumm, "occupancy/script16a_mcmcsumm1.rds")


# JAGS model for Case 2
# Case 2: estimate the weight along the diagonal of the square portion of
# alpha (alpha_diag),the weights in the morphospecies columns (alpha_morph),
# AND the weights in the non-diagonal portion of the square matrix
# AIS could never get this to run because of the error: "Slicer stuck at value with infinite density"
# Choosing to run the dynamic occupancy model with Case 1 as the next step
cl <- makeCluster(4)
out2 <- jags.parfit(cl = cl, 
            data = JAGSdata, 
            params = c("theta","alpha_diag", "alpha_morph", "alpha_nondiag") ,
            model = "occupancy/16b_neon_just_multisp_misclass_JAGS.txt", 
            # inits = JAGSinits,
            n.chains = nc,
            n.adapt = 2000,
            n.update = 2000,
            thin = nt,
            n.iter = ni) 
saveRDS(out2, "occupancy/script16b_jags_out1.rds")


# How well does the model estimate specified parameters?
library(MCMCvis)
out2_mcmcsumm <- MCMCsummary(out1)
saveRDS(out2_mcmcsumm, "occupancy/script16b_mcmcsumm2.rds")

# View JAGS output --------------------------------------------------------

out1 <- readRDS("occupancy/script16_jags_out1.rds")
out1_mcmcsumm 

# Did model converge?
hist(out1_mcmcsumm$Rhat, breaks=40)
range(out1_mcmcsumm$Rhat, na.rm=TRUE)

# Look at high Rhat values
#out_mcmcsumm <- rownames_to_column(out_mcmcsumm)
#out_mcmcsumm %>% filter(Rhat > 1.1)

# mcmc.list visualization resource: https://cran.r-project.org/web/packages/MCMCvis/vignettes/MCMCvis.html

# Look at raw numbers
# alpha_diag
MCMCsummary(out1, params = 'alpha_diag')

# alpha_morph
MCMCsummary(out1, params = 'alpha_morph')

# theta
theta_summ1 <- MCMCsummary(out1, params = 'theta', round=2)
saveRDS(theta_summ1, "occupancy/theta_summ16a.rds")
theta_summ1 <- readRDS("occupancy/theta_summ16a.rds")
hist(theta_summ1$Rhat, breaks=20)
range(theta_summ1$Rhat)
theta_summ1
range(theta_summ1$mean)
MCMCtrace(out1, params = paste0('theta\\[17,17\\]'), type = 'trace', ind = F, pdf=F, ISB=F)
MCMCsummary(out1, params='theta\\[36,[0-9]+\\]', ISB=F)

theta_df <- data.frame(expert_index = rep(1:dim(M)[1], dim(M)[2]) ,
                       paramorph_index = rep(1:dim(M)[2], each = dim(M)[1]),
                       expert_sciname = rep(rownames, dim(M)[2]),
                       para_morph = rep(colnames, each = dim(M)[1]))  %>% 
    mutate(theta_mean = theta_summ$mean)
#make para_morph a factor to force plotting in order
theta_df$para_morph = factor(theta_df$para_morph, levels=colnames) 

# Plot heatmap of theta values
ggplot(theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_mean)) + 
    geom_tile() +
    scale_fill_gradient(low="darkblue", high="white") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))
#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html


# Visualize predictions of species unobserved by expert taxonomist
par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species without expert ID")
abline(0,1)
points(x=out$mean$theta[nspec,],y=theta[nspec,],col="red")


# Scratch -----------------------------------------------------------------

# Core and chain time trial -----------------------------------------------

# alpha <- matrix(1, nrow=5,ncol=5)
# diag(alpha) <- 5 
# M <- matrix(0, nrow=5,ncol=5)
# diag(M) <- 7
# c_obs <- array(1:8,c(5,3,5,2))
# Z.init <- array(0:1,c(5,5,2))
# JAGSdata <- list(nsite = 5, 
#                  nsurv = 3, 
#                  nspec_exp = 5, #species ID'd by expert
#                  nspec_para = 5, #paratax ID's and morphospecies
#                  nyear = 2, 
#                  n = rep(7,5),
#                  alpha = alpha,
#                  M = M,
#                  c_obs = c_obs) #bundle data
# JAGSinits <- function(){ list(Z = Z.init) }
# JAGSparams <- c("psi", "lambda", "theta", "Z", "phi", "gamma", "n.occ",
# "growth", "turnover")
# #nc, ni, nb defined above with theta
# nt <- 1  #MCMC thin
# #######################
# #4 chains: 38/37sec on 4 cores, 50sec on 3 cores, 38sec/38sec on 2 cores, 72sec on 1 core
# #3 chains on 2 and 4 cores
# #3 chains: on 4 cores is 28sec, and on 2 cores is 37/42sec

# For a while, I used function arguments for jags() in jags.parfit() rather than
# the proper arguments for jags.parfit(). Here are notes that describe the
# errors I ran into
#after 22 hours of running, this error was spat out: 
#Error in unserialize(node$con) : error reading from connection
#stackexchange says its due to running out of memory
#bumping down to 2 cores now. Kicking off Monday 21:10
# with 2 burnin and 10 iterations, took less than 2 hours - but this was using
# jgas() arguments, not jags.parfit


# Create matrix of mean theta posterior
# theta_mat <- matrix(data=NA, nrow=dim(theta)[2], ncol=dim(theta)[3], byrow=T,
#                     dimnames = list(c(ind_dat %>%
#                                           filter(!is.na(expert_sciname)) %>%
#                                           select(expert_sciname) %>%
#                                           distinct() %>%
#                                           mutate_if(is.factor, as.character) %>%
#                                           arrange(expert_sciname) %>%
#                                           pull(expert_sciname)),
#                                     c(ind_dat %>%
#                                           select(para_morph) %>%
#                                           distinct() %>%
#                                           mutate_if(is.factor, as.character) %>%
#                                           arrange(para_morph) %>%
#                                           pull(para_morph)) ) )
# for (k in 1:dim(theta)[1]) {
#     theta_mat[k,] <- MCMCsummary(out, params=paste0('theta\\[',k,',[0-9]+\\]'), 
#                                  ISB=F)$mean
# }
# #why is it missing column 37?
# theta_mat
# heatmap(theta_mat, Colv = NA, Rowv = NA, scale = "none",
#         col = heat.colors(100),
#         xlab="expert taxonomist", ylab="parataxonomist", 
#         main="Theta Posterior Means")
