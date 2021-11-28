# This script creates figures for the full and reduced models and their
# comparison

library(MCMCvis)
library(MCMCpack)
library(dplyr)
library(patchwork)
library(ggmcmc)
library(readr)
library(ggplot2)

source("source/jags_input.R")

# Load data ---------------------------------------------------------------

jagsinput <- return_jags_input("full")
list2env(jagsinput, .GlobalEnv)
L_full <- L
y_df_full <- y_df
rm(jagsinput, L, y_df)

jagsinput <- return_jags_input("reduced")
list2env(jagsinput, .GlobalEnv)
L_red <- L
y_df_red <- y_df
rm(jagsinput, L, y_df)

full_jm <- readRDS("output/full_jm.rds")
full_jm_summ <- readRDS("output/full_jmsumm.rds")
full_theta_summ <- readRDS("output/full_theta_summ.rds") %>%
  rename(top = '97.5%', bottom= '2.5%') %>%
  mutate(fullCIwidth = top - bottom ) 

red_jm <- readRDS("output/reduced_jm.rds")
red_jm_summ <- readRDS("output/reduced_jmsumm.rds")
red_theta_summ <- readRDS("output/reduced_theta_summ.rds")  %>%
  rename(top = '97.5%', bottom= '2.5%') %>%
  mutate(redCIwidth = top - bottom )


# Figures for full model --------------------------------------------------

# theta
theta_prior <- MCMCpack::rdirichlet(1000, c(100,rep(1.5,dim(alpha)[2]-1) ))

theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       para_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_sciname = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_median = full_theta_summ$"50%")
#make para_sciname a factor to force plotting in order
theta_df$para_sciname = factor(theta_df$para_sciname, levels=colnames) 

# Plot heatmap of theta values
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442") #, "#0072B2", "#D55E00", "#CC79A7")
theta_df %>%
  mutate(theta_med_NA=ifelse(theta_median<.04,NA,theta_median)) %>%
  #mutate(theta_med_disc = as.factor(floor(theta_median*10)/10)) %>% 
  ggplot(aes(x=para_sciname, y=expert_sciname, fill= theta_median)) + 
  geom_tile() +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) + #, "#0072B2", "#D55E00", "#CC79A7")) +
  #scale_fill_brewer(palette = "Set1") +
  #scale_fill_fermenter(type="div") +
  #scale_fill_stepsn("Posterior\nmedian", colours=c("#999999", "#E69F00"),
  #                  values=c(0,0.1,0.3,0.5,0.7,0.9,1),
  #                   breaks=c(0,0.1,0.3,0.5,0.7,0.9,1)) +
  scale_fill_binned("Posterior\nmedian",type="viridis", option="plasma",
                    breaks=c(0,0.1,0.3,0.5,0.7,0.9),
                    guide = guide_coloursteps(even.steps = FALSE)) +
  #scale_fill_viridis_c("Posterior\nmedian") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Imperfect classification") + ylab("True classification") +
  scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname)))) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/thetaconfusion.png", width=7.9, height=5.4)


# Plot all species' occupancy through seasons
psi_full_summary <- MCMCsummary(full_jm, params = 'psi', round=2) %>%
  rownames_to_column()

psi_full_summary %>%
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
ggsave("figures/occthroughtime.png", width=11, height=8)

psi_red_summary <- MCMCsummary(red_jm, params = 'psi', round=2) %>%
  rownames_to_column()

psi_red_summary %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
  mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
  mutate(site = dimnames(z.init)[[1]][site], 
         species = dimnames(z.init)[[2]][species], 
         year = dimnames(z.init)[[3]][year]) %>%
  #filter(species=="Amara convexa") %>%
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
  full_jm %>%
    lapply(function(x) {
      cols <- paste0("Tau_", level) %>%
        grepl(colnames(x))
      as.data.frame(x[, cols]) %>%
        mutate(iter = 1:n()) %>%
        tidyr::pivot_longer(cols = -iter) %>%
        separate(name, into = c("row", "col"), sep = ",") %>%
        mutate_at(c("row", "col"), readr::parse_number)
    }) %>%
    bind_rows(.id="chain") %>%
    unite("iter", c("chain", "iter"), sep = "_") %>%
    split(.$iter) %>%
    lapply(FUN = function(d) {
      Tau <- reshape2::acast(data = dplyr::select(d, -iter), row~col, value.var="value")
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

eps_d <- ggmcmc::ggs(full_jm, family = "eps_") %>%
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
    dplyr::select(-col) %>%
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
ggsave("figures/ranefs-site-level.png", width = 7, height = 5)

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
ggsave("figures/ranefs-spec-level.png", width = 7, height = 5)

ggsave("figures/ranefs.png")


# Figures for reduced model and comparing full/reduced models ---------------------------------

# Visualize theta
red_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                           para_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                           expert_sciname = rep(rownames, dim(alpha)[2]),
                           para_sciname = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_mean = red_theta_summ$mean,
         theta_median = red_theta_summ$"50%")
#make para_sciname a factor to force plotting in order
red_theta_df$para_sciname = factor(red_theta_df$para_sciname, levels=colnames) 

# Plot heatmap of theta values
ggplot(red_theta_df, aes(x=para_sciname, y=expert_sciname, fill= theta_median)) + 
  geom_tile() +
  scale_fill_gradient(limits=c(0,1),low="darkblue", high="white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_discrete(limits = rev(levels(as.factor(red_theta_df$expert_sciname)))) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
ggsave("figures/reduced_thetaconfusion.png", width=7, height=5)

# Compare density plots of theta prior, full model, reduced model 
red_theta_post_mat <- MCMCchains(red_jm,params='Theta')
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

full_theta_post_mat <- MCMCchains(full_jm,params="Theta")
full_post_df <- as.data.frame(full_theta_post_mat) %>%
  pivot_longer(colnames(full_theta_post_mat),names_to="index") %>%
  mutate(model = "full")

theta_df <- rbind(theta_prior, reduced_post_df, full_post_df)

# Plot matrix of theta prior and poserior densities
theta_df %>% 
  filter(index==c("Theta[48,48]","Theta[48,49]","Theta[48,50]","Theta[48,51]",
                  "Theta[49,48]","Theta[49,49]","Theta[49,50]","Theta[49,51]",
                  "Theta[50,48]","Theta[50,49]","Theta[50,50]","Theta[50,51]",
                  "Theta[51,48]","Theta[51,49]","Theta[51,50]","Theta[51,51]")) %>%
  mutate(index2 = index) %>%
  separate("index2", into = c("exp_ind", "para_ind"), sep = ",") %>%
  mutate_at(c('exp_ind', 'para_ind'), readr::parse_number) %>%
  mutate(exp = rownames[exp_ind],
         para = rownames[para_ind] ,
         exp_abbr = ifelse(exp=="Pterostichus (Hypherpes) sp.", "P. (Hypherpes) sp.",
                           ifelse(exp=="Pterostichus adstrictus","P. adstrictus",
                                  ifelse(exp=="Pterostichus melanarius","P. melanarius","P.restrictus"))),
         para_abbr = ifelse(para=="Pterostichus (Hypherpes) sp.", "P. (Hypherpes) sp.",
                            ifelse(para=="Pterostichus adstrictus","P. adstrictus",
                                   ifelse(para=="Pterostichus melanarius","P. melanarius","P.restrictus"))),
         names = paste0("[",exp_abbr,", ",para_abbr,"]")) %>%
  ggplot(aes(x=value, y=..scaled.., fill=model)) +
  geom_density( alpha=0.5 ) + 
  #facet_wrap(exp_abbr ~ para_abbr, scales="free_x") +
  #facet_grid(rows=vars(exp), cols=vars(para), scales="free") +
  #facet_wrap( ~ index, scales="free_x") +
  facet_wrap( ~ names, scales="free_x") +
  xlab("Theta") + ylab("Density") + scale_y_continuous(breaks=seq(0, 1, 0.5))
ggsave("figures/comparedensities.png")

# Heat map of posterior differences between two confusion matrices 
# Consider the added value of full occupancy model vs just misclass model  

full_theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                            para_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                            expert_sciname = rep(rownames, dim(alpha)[2]),
                            para_sciname = rep(colnames, each = dim(alpha)[1]))  %>% 
  mutate(theta_median = full_theta_summ$"50%")
#make para_sciname a factor to force plotting in order
full_theta_df$para_sciname = factor(full_theta_df$para_sciname, levels=colnames) 

theta_med_diff_df <- red_theta_df %>%
  rename(red_theta_med = theta_median) %>%
  left_join(full_theta_df %>% 
              select(expert_sciname, para_sciname, full_theta_med = theta_median)) %>%
  mutate(median_diff = full_theta_med - red_theta_med,
         diff_pos = ifelse(median_diff == 0, NA, ifelse(median_diff > 0.01, T, F)))
# change 0.01 to whatever desired threshold

ggplot(theta_med_diff_df, aes(x=para_sciname, y=expert_sciname, fill= diff_pos)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("blue", "red"), "Full > Reduced \nTheta Median") +
  xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
  scale_y_discrete(limits = rev(levels(as.factor(theta_med_diff_df$expert_sciname))))
ggsave("figures/thetadifference.png")



# Compare theta precision between reduced and full model
# Scatterplot of 95% CI widths (x-axis full, y-axis reduced)  

theta_CIwidth <- full_theta_summ %>%
  left_join(red_theta_summ %>% dplyr::select(rowname,redCIwidth))

theta_CIwidth %>%
  separate("rowname", into = c("expertID_idx", "parataxID_idx"), sep = ",", remove=F) %>%
  mutate_at(c('expertID_idx', 'parataxID_idx'), readr::parse_number) %>%
  left_join(y_df_full %>%
              count(parataxID_idx) ) %>%
  as_tibble() %>%
  ggplot(aes(x=fullCIwidth, y=redCIwidth, col=log(n))) +
  geom_point() +
  geom_text(aes(label=ifelse(parataxID_idx==expertID_idx,rowname,""),hjust=0,vjust=0)) +
  scale_color_gradient(low = "orange", high = "darkblue")  +
  xlab("Full Model CI width") + ylab("Reduced Model CI width") +
  xlim(c(0,0.11)) + ylim(c(0,0.11)) +
  geom_abline(intercept=0,slope=1)

theta_CIwidth %>%
  separate("rowname", into = c("expertID_idx", "parataxID_idx"), sep = ",", remove=F) %>%
  mutate_at(c('expertID_idx', 'parataxID_idx'), readr::parse_number) %>%
  left_join(y_df_full %>%
              count(parataxID_idx) ) %>%
  as_tibble() %>%
  ggplot(aes(x=fullCIwidth, y=redCIwidth)) +
  geom_point(alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "darkblue")  +
  xlab("Full Model CI width") + ylab("Reduced Model CI width") +
  xlim(c(0,0.11)) + ylim(c(0,0.11)) +
  geom_abline(intercept=0,slope=1)
ggsave("figures/CIwidthcomparison.png")




### Table in Appendix S2
all_paratax_by_ind <- readRDS("data/all_paratax_df.rds") %>%
  mutate(parataxID = ifelse(is.na(morphospeciesID), scientificName, morphospeciesID)) %>% 
  uncount(individualCount) %>%
  rownames_to_column() 
all_paratax_by_ind$exp_sciname <- NA
expert_df <- readRDS("data/expert_df.rds") 
pinned_df <- readRDS("data/pinned_df.rds") 
expert_pinned_df <- expert_df %>%
  rename(exp_sciname = scientificName) %>%
  left_join(pinned_df %>% 
              mutate(parataxID = ifelse(is.na(morphospeciesID), scientificName, morphospeciesID)) %>% 
              dplyr::select(individualID, subsampleID, parataxID))#para_morph_combo = scimorph_combo))
rm(expert_df, pinned_df)

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
para_new <- para_new %>% as_tibble()

para_new %>%
  #filter to samples that had expert ID
  filter(!is.na(exp_sciname)) %>%
  group_by(parataxID, exp_sciname) %>%
  summarize(n=n()) %>%
  arrange(exp_sciname,parataxID) %>%
  mutate(match=ifelse(parataxID==exp_sciname,1,0)) %>%
  data.frame()



# # Compare theta median between reduced and full model
# red_theta_summ <- readRDS("output/reduced_theta_summ.rds") %>%
#   mutate(model="reduced")
# theta_medians <- readRDS("output/full_theta_summ.rds") %>%
#   mutate(model="full") %>%
#   full_join(red_theta_summ) %>%
#   as_tibble()
# #AIS bring these up top
# 
# theta_medians %>%
#   dplyr::select(rowname, model, "50%") %>%
#   pivot_wider(names_from=model, values_from="50%") %>%
#   separate("rowname", into = c("expertID_idx", "parataxID_idx"), sep = ",", remove=F) %>%
#   mutate_at(c('expertID_idx', 'parataxID_idx'), readr::parse_number) %>%
#   left_join(y_df_full %>%
#             count(parataxID_idx) ) %>%
#   ggplot(aes(x=full, y=reduced, col=log(n))) +
#   geom_point() +
#   geom_text(aes(label=ifelse(log(n)>4,rowname,"")),hjust=0,vjust=0) +
#   scale_color_gradient(