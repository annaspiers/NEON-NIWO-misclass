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
  left_join(pinned_df %>% dplyr::select(individualID, subsampleID, para_morph_combo = scimorph_combo))

# subsample aRw+8bWvW/wOIENap0etmXiCq6V1RZPfCkRcy1w8BLk= had 
# 6 individuals sorted, but 8 pinned and expertly identified. 
# This is the only subsample with less sorted individuals than expertly ID'ed
# As a workaround, we'll remove two expertly identified inidividuals
# Update, this should be fixed in the NEON data portal week of 8/24
removed_inds <- expert_pinned_df %>%
  filter(subsampleID == "aRw+8bWvW/wOIENap0etmXiCq6V1RZPfCkRcy1w8BLk=") %>%
  pull(individualID) %>%
  tail(2)
expert_pinned_df <- expert_pinned_df %>%
  filter(!(individualID %in% removed_inds))

rm(expert_df, pinned_df, removed_inds)

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
# visualize alpha
test_alpha <- MCMCpack::rdirichlet(1000, c(alpha[1,1], rep(alpha[1,2], length(colnames)-1)))
hist(test_alpha[,-1], xlim=c(0,1),col="red")
hist(test_alpha[,1],  add=T)


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

assertthat::assert_that(sum(z.dat_cast, na.rm=T) == nrow(expert_pinned_df %>% distinct(plotID, exp_sciname, col_year)))

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
      if (sum(L[i,,t], na.rm = TRUE) > 0 ) {
        ifelse(z.init[i,k,t] == 0, 1, z.init[i,k,t])
      }
    }
  }
}
rm(i,k,t)

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
jags_d <- list(nsite = dim(L)[1],
               nsurv = dim(L)[2], 
               nyear = dim(L)[3], 
               K_exp = dim(alpha)[1], 
               K_para = dim(alpha)[2],
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = y_df$plotID_idx, #needs to be numeric
               year = as.numeric(as.factor(y_df$col_year)), #needs to be numeric starting from 1
               # if the individual was labeled by the expert, true ID is known
               k = y_df$expertID_idx,
               # for all individuals, we get paratxonomist IDs
               y = y_df$parataxID_idx,
               z = z.dat,
               R = diag(rep(1, 4)))
JAGSinits <- function(){ list(z = z.init) }
nc <- 4
ni <- 100000
cl <- makeCluster(nc)
jm <- jags.parfit(cl = cl,
                  data = jags_d,
                  params = c("logit_psi1", "logit_phi", "logit_gamma",
                             "psi",
                             "lambda", "Theta",
                             "eps_site", "eps_spec", "Tau_spec", "Tau_site",
                             "log_growth", "turnover"),
                  model = "neon_dynamic_disagg_multisp_JAGS.txt",
                  n.chains = nc,
                  n.adapt = 6000,
                  n.update = 6000,
                  thin = ni/1000,
                  n.iter = ni) 
dir.create("output", showWarnings = FALSE)
saveRDS(jm, "output/jm2.rds")
jm_summ <- MCMCsummary(jm)
saveRDS(jm_summ, "output/jm_summ2.rds")

# View JAGS jmput --------------------------------------------------------

jm <- readRDS("output/jm-allspec.rds")
jm_summ <- readRDS("output/jm_mcmcsumm-allspec.rds")

# Did model converge?
hist(jm_summ$Rhat, breaks=40)
range(jm_summ$Rhat, na.rm=TRUE)

# Look at high Rhat values
jm_summ <- rownames_to_column(jm_summ)
jm_summ %>% filter(Rhat > 1.1)

#MCMCtrace(jm, params = paste0('eps_spec\\[1,4\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 

# Look at raw numbers
# theta
theta_summ <- MCMCsummary(jm, params = 'Theta', round=2)
saveRDS(theta_summ, "output/theta_summ.rds")
theta_summ <- readRDS("output/theta_summ-allspec.rds")
hist(theta_summ$Rhat)
hist(theta_summ$mean)
range(theta_summ$Rhat, na.rm=TRUE)
theta_summ <- rownames_to_column(theta_summ)

# # Look at some theta posteriors with high Rhat values
MCMCtrace(jm, params = paste0('Theta\\[1,1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
theta_prior <- MCMCpack::rdirichlet(1000, c(100,rep(1.5,dim(alpha)[2]-1) ))
plot(density(theta_prior[, 1]))


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
    scale_fill_viridis_c("Posterior\nmean") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname))))

# phi - survival probability
MCMCsummary(jm, params = 'logit_phi', round=2)
MCMCtrace(jm, params = 'logit_phi', type = 'trace', ind = F, pdf=F)

# gamma - colonization probability
MCMCsummary(jm, params = 'logit_gamma', round=2)
MCMCtrace(jm, params = 'logit_gamma', type = "trace", ind = F, pdf=F)

# psi - occupancy prob.
(psi_summary <- MCMCsummary(jm, params = 'psi', round=2))
MCMCtrace(jm, params = 'psi', type = 'trace', ind = F, pdf=F)

# Plot all species' occupancy through seasons
rownames_to_column(psi_summary) %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
  mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
  mutate(site = dimnames(z.init)[[1]][site], 
         species = dimnames(z.init)[[2]][species], 
         year = dimnames(z.init)[[3]][year]) %>%
  ggplot(aes(year, mean, group = site)) + 
  geom_line() + 
  facet_wrap(~species) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), color = NA, alpha = .05) +
  xlab("Year") + 
  ylab("Occupancy probability") + 
  theme_minimal()


# Visualize relationships among site and species level parameters ---------

get_rho_df <- function(level) {
  stopifnot(level %in% c("site", "spec"))
  jm %>%
    lapply(function(x) {
      cols <- paste0("Tau_", level) %>%
        grepl(colnames(x))
      as.data.frame(x[, cols]) %>%
        mutate(iter = 1:n()) %>%
        pivot_longer(cols = -iter) %>%
        separate(name, into = c("row", "col"), sep = ",") %>%
        mutate_at(c("row", "col"), readr::parse_number)
    }) %>%
    bind_rows(.id="chain") %>%
    unite("iter", c("chain", "iter"), sep = "_") %>%
    split(.$iter) %>%
    lapply(FUN = function(d) {
      Tau <- reshape2::acast(data = select(d, -iter), row~col, value.var="value")
      Sigma <- solve(Tau)
      Rho <- cov2cor(Sigma)
      reshape2::melt(Rho, varnames = c("row", "col"), value.name = "rho") %>%
        as_tibble %>%
        filter(row > col)
    }) %>%
    bind_rows(.id = "iter")
}


rho_site <- get_rho_df(level = "site")
rho_spec <- get_rho_df(level = "spec")

rho_site %>%
  ggplot(aes(rho)) + 
  geom_density() + 
  facet_grid(row~col)

rho_spec %>%
  ggplot(aes(rho)) + 
  geom_density() + 
  facet_grid(row~col)

eps_d <- ggmcmc::ggs(jm, family = "eps_") %>%
  mutate(level = ifelse(grepl("site", Parameter), "site", "spec")) %>%
  separate(Parameter, into = c("row", "col"), sep = ",") %>%
  mutate_at(c('row', 'col'), parse_number) %>%
  unite("iter", c("Iteration" , "Chain")) %>%
  group_by(row, col, level) %>%
  summarize(med = median(value)) %>%
  ungroup %>%
  mutate(param = c("e_psi1", "e_phi", "e_gamma", "e_lambda")[col])

diag_hists <- function(lev, column, xmin=-5, xmax=5, binwidth=.25, 
                       xlabel) {
  stopifnot(lev %in% eps_d$level)
  df <- eps_d %>%
    filter(level == lev, 
           col == column)
  df %>%
    ggplot(aes(x=med)) + 
      geom_histogram(binwidth = binwidth, fill = "#1F968BFF") + 
      xlab(xlabel) + 
      ylab("") + 
      xlim(xmin, xmax) + 
      theme(axis.text.y = element_blank())
}

diag_hists(lev = "site", column = 1, 
           xlabel = as.expression(bquote(epsilon^(psi[1]))))
diag_hists(lev = "spec", column = 1, 
           xlabel = as.expression(bquote(epsilon^(psi[1]))))

lower_diag <- function(d1, d2, lev, xmin=-5, xmax=5) {
  df <- eps_d %>%
    select(-col) %>%
    filter(level == lev) %>%
    pivot_wider(names_from = param, values_from = med)
  df$x <- df[[d1]]
  df$y <- df[[d2]]
  
  df %>%
    ggplot(aes(x, y)) + 
    geom_point(color = "#482677FF") + 
    xlab("") + 
    ylab("") + 
    xlim(xmin, xmax) + 
    ylim(xmin, xmax)
}

lower_diag("e_psi1", "e_lambda", lev = "site")
lower_diag("e_lambda", "e_psi1", lev = "site")


upper_diag <- function(d, r, c, xlabel) {
  stopifnot(r > c)
  d %>%
    filter(row == r, col == c) %>%
    ggplot(aes(rho)) + 
      geom_density(fill = "#482677FF", alpha = .7) + 
      xlab(xlabel) + 
      ylab("") + 
      xlim(-1, 1) + 
      theme(axis.text.y = element_blank())
}

upper_diag(rho_spec, r = 4, c = 3, 
           xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi)))))
upper_diag(rho_site, r = 4, c = 3, 
           xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi)))))

# Cobble together all of the panels (site level)
site_xmin <- -1.7
site_xmax <- 1.2

wrap_plots(
  list(
    # Column 1
    diag_hists(lev = "site", column = 1, 
               xlabel = as.expression(bquote(epsilon^(psi[1]))), 
               xmin = site_xmin, xmax = site_xmax), 
    lower_diag("e_psi1", "e_phi", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_psi1", "e_gamma", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_psi1", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    # Column 2
    upper_diag(rho_site, r = 2, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(phi), epsilon^(psi[1]))))),
    diag_hists(lev = "site", column = 2, 
               xlabel = as.expression(bquote(epsilon^(phi))), 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_phi", "e_gamma", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_phi", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    # Column 3
    upper_diag(rho_site, r = 3, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(psi[1]))))),
    upper_diag(rho_site, r = 3, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(phi))))),
    diag_hists(lev = "site", column = 3, 
               xlabel = as.expression(bquote(epsilon^(gamma))), 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_gamma", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax), 
    # Column 4
    upper_diag(rho_site, r = 4, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(psi[1]))))),
    upper_diag(rho_site, r = 4, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi))))),
    upper_diag(rho_site, r = 4, c = 3, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(gamma))))),
    diag_hists(lev = "site", column = 4, 
               xlabel = as.expression(bquote(epsilon^(lambda))), 
               xmin = site_xmin, xmax = site_xmax)
    ), 
  nrow = 4, byrow = FALSE)
dir.create("fig", showWarnings = FALSE)
ggsave("fig/site-level-ranefs.pdf", width = 7, height = 5)

# Cobble together all of the panels (species level)
spec_xmin <- -4.7
spec_xmax <- 1.8

wrap_plots(
  list(
    # Column 1
    diag_hists(lev = "spec", column = 1, 
               xlabel = as.expression(bquote(epsilon^(psi[1]))), 
               xmin = spec_xmin, xmax = spec_xmax), 
    lower_diag("e_psi1", "e_phi", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_psi1", "e_gamma", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_psi1", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    # Column 2
    upper_diag(rho_spec, r = 2, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(phi), epsilon^(psi[1]))))),
    diag_hists(lev = "spec", column = 2, 
               xlabel = as.expression(bquote(epsilon^(phi))), 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_phi", "e_gamma", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_phi", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    # Column 3
    upper_diag(rho_spec, r = 3, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(psi[1]))))),
    upper_diag(rho_spec, r = 3, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(phi))))),
    diag_hists(lev = "spec", column = 3, 
               xlabel = as.expression(bquote(epsilon^(gamma))), 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_gamma", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax), 
    # Column 4
    upper_diag(rho_spec, r = 4, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(psi[1]))))),
    upper_diag(rho_spec, r = 4, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi))))),
    upper_diag(rho_spec, r = 4, c = 3, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(gamma))))),
    diag_hists(lev = "spec", column = 4, 
               xlabel = as.expression(bquote(epsilon^(lambda))), 
               xmin = spec_xmin, xmax = spec_xmax)
  ), 
  nrow = 4, byrow = FALSE)
ggsave("fig/spec-level-ranefs.pdf", width = 7, height = 5)






# lambda - expected abundance, given occupancy, dim: nsite x nsurv x nspec x nyear
lambda_summ <- MCMCsummary(jm, params = 'lambda', round=2) #takes 5ish min
saveRDS(lambda_summ, "output/lambda_summ.rds")
lambda_summ <- readRDS("output/lambda_summ.rds")
hist(lambda_summ$Rhat)
range(lambda_summ$Rhat)
lambda_summ
range(lambda_summ$mean)
MCMCtrace(jm, params = 'lambda', type = 'density', ind = F, pdf=F)
range(lambda_summ$mean)
#a good handful of Rhats > 1.1, as large as 3.32

# Z
Z_summ <- MCMCsummary(jm, params = 'Z', round=2)
saveRDS(Z_summ, "output/Z_summ.rds")
#Z_summ <- readRDS("output/theta_summ.rds")
MCMCtrace(jm, params = 'Z', type = 'density', ind = F, pdf=F)

Z.dat[is.na(Z.dat)] <- 0
Z.init[is.na(Z.init)] <- 0
Z.prior <- Z.dat+Z.init
plot(Z.prior, Z_summ$mean)
#WORK ON ABOVE PLOT - MAKE SURE THEY PLOT 1-1 ACCURATELY, now it's not accurate
plot(apply(c_obs, c(1,3,4), max, na.rm = TRUE), jm$mean$Z)











# n.occ
print(jm$summary[grep("n.occ", row.names(jm$summary)), c(1, 2, 3, 7)], dig = 2)

# growth
print(jm$summary[grep("growth", row.names(jm$summary)), c(1, 2, 3, 7)], dig = 2) 
jm$mean$growth

# turnover
print(jm$summary[grep("turnover", row.names(jm$summary)), c(1, 2, 3, 7)], dig = 2) 
jm$mean$turnover

# Look at psi, phi, gamma, growth, turnover, n.occ graphically
species <- as.character(unique(sample_dat$para_sciname))
for (k in 1:length(species)) {
    k_spec <- species[k]
    
    p1 <- ggplot(data = sample_dat %>% filter(para_sciname == k_spec)) + 
        geom_col(aes(x=collectDate,y=sp_abund)) +
        ggtitle(k_spec) + xlab("Collection Date") + ylab("Abundance")
    p2 <- ggplot(data = left_join(sample_dat %>% dplyr::select(col_year) %>% distinct(),
                                  sample_dat %>% filter(para_sciname==k_spec, sp_abund > 0) %>% 
                                      group_by(col_year) %>% 
                                      summarize(n_plots=n_distinct(plotID)))) + 
        geom_line(aes(x=col_year, y=n_plots),) +
        geom_line(aes(x=col_year,y=jm$mean$n.occ[k,]),col="red") +
        annotate("text", x=2017, y=50, label = "Predicted occupied (n.occ)", col="red") +
        annotate("text", x=2017, y=45, label = "Observed") +
        ggtitle("Site occupancy/detection") + xlab("Year") + ylab("Number of sites")
    p3 <- ggplot() + geom_line(aes(x=unique(sample_dat$col_year),
                                   y=jm$mean$psi[k,]), col="red") +
        ggtitle("Occupancy (psi)") + xlab("Year") + ylab("Probability") + ylim(c(0,1))
    p4 <- ggplot() + 
        geom_line(aes(x=unique(sample_dat$col_year),y=jm$mean$growth[k,]),col="blue") + 
        annotate("text", x=2017, y=0.6, label = "Growth", col="blue") +
        geom_line(aes(x=unique(sample_dat$col_year),y=c(jm$mean$turnover[k,],NA)),col="darkgreen") + 
        annotate("text", x=2017, y=0.5, label = "Turnover", col="darkgreen") +
        geom_hline(yintercept=jm$mean$phi[k],col="purple") + 
        annotate("text", x=2017, y=0.4, label = "Survival (phi)", col="purple") +
        geom_hline(yintercept=jm$mean$gamma[k],col="orange") + 
        annotate("text", x=2017, y=0.3, label = "Colonization (gamma)", col="orange") +
        ggtitle("Demographic rates") + xlab("Year") + ylab("Rate") + ylim(c(0,1))
    
    grid.arrange(p1,p2,p3,p4,nrow=2)
}
par(mfrow=c(1,1)) #reset plotting

# Visualize predictions of species unobserved by expert taxonomist
par(mfrow=c(1,1))
plot(NA,xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted",ylab="Observed",main="Theta comparison for species withjm expert ID")
abline(0,1)
points(x=jm$mean$theta[nspec,],y=theta[nspec,],col="red")

