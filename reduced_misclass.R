# Here we create a reduced model estimating only misclassification probabilities
# without occupancy parameters using NEON Niwot Ridge carabid data 2015-2018

library(tidyr) #uncount()
library(tibble) #rownames_to_column()
library(dclone) #jags.parfit()
library(MCMCvis)
library(MCMCpack)
library(dplyr)
library(patchwork)
library(ggmcmc)
library(readr)
library(ggplot2)
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
diag(alpha) <- 200


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
                  model = "reduced_misclass_JAGS.txt",
                  n.chains = nc,
                  n.adapt = 2000,
                  n.update = 2000,
                  thin = ni/1000,
                  n.iter = ni) 

jm_summ <- MCMCsummary(jm)

saveRDS(jm, "output/reduced_jm.rds")
saveRDS(jm_summ, "output/reduced_jmsumm.rds")

# View JAGS output --------------------------------------------------------

jm <- readRDS("output/reduced_jm.rds")
jm_summ <- readRDS("output/reduced_jmsumm.rds")
# Did model converge?
hist(jm_summ$Rhat, breaks=40)
jm_summ <- rownames_to_column(jm_summ)

### Look at raw numbers
# THETA
red_theta_summ <- MCMCsummary(jm, params = 'Theta')
red_theta_summ <- rownames_to_column(red_theta_summ)
saveRDS(red_theta_summ, "output/reduced_theta_summ.rds")
#red_theta_summ <- readRDS("output/reduced-theta_summ.rds")

# Visualize theta
red_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_mean = red_theta_summ$mean,
         theta_median = red_theta_summ$"50%")
#make para_morph a factor to force plotting in order
red_theta_df$para_morph = factor(red_theta_df$para_morph, levels=colnames) 
# Plot heatmap of theta values
theta_med <- ggplot(red_theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_median)) + 
  geom_tile() +
  scale_fill_gradient(limits=c(0,1),low="darkblue", high="white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = rev(levels(as.factor(red_theta_df$expert_sciname))))
#consider on plotly: https://www.r-graph-gallery.com/79-levelplot-with-ggplot2.html
png("figures/reduced_thetaconfusion.png", width=850, height=600)
print(theta_med)
dev.off()

# Compare density plots of theta prior, full model, reduced model 
red_theta_post_mat <- MCMCchains(jm,params='Theta')
reduced_post_df <- as.data.frame(red_theta_post_mat) %>%
  pivot_longer(colnames(red_theta_post_mat),names_to="index") %>%
  mutate(model = "reduced")

theta_prior <- data.frame(index = rep(colnames(red_theta_post_mat),each=nrow(red_theta_post_mat)),
                          value = NA,
                          model = "prior")
value_vec <- c()
for (i in 1:ncol(alpha)) {
  temp_vec <- as.vector(MCMCpack::rdirichlet(4000, alpha[,i]))
  value_vec <- c(value_vec, temp_vec)
}
theta_prior$value <- value_vec
theta_prior$index <- as.character(theta_prior$index)

full_model <- readRDS("output/full_jm.rds")
full_theta_post_mat <- MCMCchains(full_model,params="Theta")
full_post_df <- as.data.frame(full_theta_post_mat) %>%
  pivot_longer(colnames(full_theta_post_mat),names_to="index") %>%
  mutate(model = "full")

theta_df <- rbind(theta_prior, reduced_post_df, full_post_df)

# Plot matrix of theta prior and poserior densities
png("figures/comparedensities.png")
print(ggplot(theta_df %>% filter(index==c(#"Theta[45,45]","Theta[45,46]","Theta[45,47]","Theta[45,48]",
                                          "Theta[46,46]","Theta[46,47]","Theta[46,48]",
                                          "Theta[47,46]","Theta[47,47]","Theta[47,48]",
                                          "Theta[48,46]","Theta[48,47]","Theta[48,48]")), 
             aes(x=value,y=..scaled..,fill=model)) +
        geom_density( alpha=0.4,color=NA) + 
        facet_wrap( ~ index, scales="free_x") +
        #facet_wrap( ~ index) +
        #  xlim(c(0,1)) +
        xlab("Theta") + scale_y_continuous(breaks=seq(0, 1, 0.5)))
dev.off()

# Heat map of posterior differences between two confusion matrices 
# Consider the added value of full occupancy model vs just misclass model  
red_theta_df

full_theta <- readRDS("output/full_theta_summ.rds")
full_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_median = full_theta$"50%")
#make para_morph a factor to force plotting in order
full_theta_df$para_morph = factor(full_theta_df$para_morph, levels=colnames) 

theta_med_diff_df <- red_theta_df %>%
  rename(red_theta_med = theta_median) %>%
  left_join(full_theta_df %>% select(expert_sciname, para_morph, full_theta_med = theta_median)) %>%
  mutate(median_diff = full_theta_med - red_theta_med,
         diff_pos = ifelse(median_diff == 0, NA, ifelse(median_diff > 0, T, F)))

png("figures/thetadifference.png", width=850, height=600)
print(ggplot(theta_med_diff_df, aes(x=para_morph, y=expert_sciname, fill= diff_pos)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = c("blue", "red"), "Full > Reduced \nTheta Median") +
  xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
  scale_y_discrete(limits = rev(levels(as.factor(theta_med_diff_df$expert_sciname)))))
dev.off()



# Compare precision between reduced and full model
# Scatterplot of 95% CI widths (x-axis full, y-axis reduced)  
red_theta_summ <- readRDS("output/reduced_theta_summ.rds") %>%
  rename(top = '97.5%', bottom= '2.5%') %>%
  mutate(redCIwidth = top - bottom )
full_theta <- readRDS("output/full_theta_summ.rds")%>%
  rownames_to_column() %>%
  rename(top = '97.5%', bottom= '2.5%') %>%
  mutate(fullCIwidth = top - bottom ) %>%
  left_join(red_theta_summ %>% dplyr::select(rowname,redCIwidth))

ggplot(full_theta) +
  geom_point(aes(x=fullCIwidth, y=redCIwidth), alpha=0.4) +
  xlab("Full Model CI width") + ylab("Reduced Model CI width") +
  xlim(c(0,0.11)) + ylim(c(0,0.11)) +
  geom_abline(intercept=0,slope=1)
ggsave("figures/CIwidthcomparison.png")



# Median

red_theta_summ <- readRDS("output/reduced_theta_summ.rds") %>%
  mutate(model="reduced")
theta_medians <- readRDS("output/full_theta_summ.rds")%>%
  rownames_to_column() %>%
  mutate(model="full") %>%
  full_join(red_theta_summ) %>%
  as_tibble()

theta_medians %>%
  dplyr::select(rowname, model, "50%") %>%
  pivot_wider(names_from=model, values_from="50%") %>%
  ggplot() +
  geom_point(aes(x=full, y=reduced), alpha=0.4) +
  xlab("Full Model Median") + ylab("Reduced Model Median") +
  #xlim(c(0,0.11)) + ylim(c(0,0.11)) +
  geom_abline(intercept=0,slope=1)
